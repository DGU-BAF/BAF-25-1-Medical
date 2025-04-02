#########################################################################
##### Pathway Score Code
#########################################################################
rm(list=ls())
getwd()
setwd("C:/Users/hemoa/OneDrive/바탕 화면/비어플/[25-1] Medical/dataset")

# 필요한 패키지 로드
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)

#########################################################################
##### 형식에 맞게 전처리 진행

# 1. 데이터 불러오기
ccle <- read.csv("./preprocessing/#_filtered_CCLE_gene_expression.csv", row.names = 1, check.names = FALSE)
kegg <- read.csv("./preprocessing/#_filtered_Kegg_info.csv", stringsAsFactors = FALSE)

# 2. gene expression long format으로 melt
ccle_long <- ccle %>%
  tibble::rownames_to_column("SampleID") %>%
  pivot_longer(-SampleID, names_to = "GeneSymbol", values_to = "Expression")

# 3. NodeInfo 파싱해서 gene-level로 explode
kegg_long <- kegg %>%
  mutate(NodeInfo = str_replace_all(NodeInfo, "\\[|\\]|'", "")) %>%  
  separate_rows(NodeInfo, sep = ",\\s*") %>%  # gene 하나씩 분리
  rename(GeneSymbol = NodeInfo)

#########################################################################
##### Pathway score 계산

# 4. Pathway score 계산 함수
stand.path.score <- function(exp.df, path.genes) {
  path.exp <- exp.df %>% filter(GeneSymbol %in% path.genes)
  
  if (nrow(path.exp) < 3) return(NULL)  # 너무 적으면 skip
  
  path.wide <- path.exp %>%
    pivot_wider(names_from = SampleID, values_from = Expression) %>%
    column_to_rownames("GeneSymbol")
  
  path.mat <- as.matrix(path.wide)
  
  # 기준 유전자: 첫 번째 유전자
  ref.gene <- rownames(path.mat)[1]
  # 1. 유전자 상관계수 방향
  gene.corr <- apply(path.mat, 1, function(row) sign(cor(path.mat[ref.gene, ], row)))
  # 2. 랭크 매트릭스
  rank.mat <- apply(abs(path.mat), 1, rank)
  # 3. Pathway Score 계산
  path.score <- colSums(t(rank.mat) * gene.corr) / nrow(path.mat)
  # 4. 정규화
  standardized <- scale(path.score)[,1]
  return(standardized)
}

# 5. 전체 Pathway에 대해 루프
sample_ids <- unique(ccle_long$SampleID)
pathway_list <- unique(kegg_long$PathwayID)
score_df <- data.frame(SampleID = sample_ids)

for (pw in pathway_list) {
  genes <- kegg_long %>% filter(PathwayID == pw) %>% pull(GeneSymbol)
  score <- stand.path.score(ccle_long, genes)
  if (!is.null(score)) {
    score_df[[pw]] <- score
  }
}

# 6. 저장
write.csv(score_df, "#_CCLE_pathway_scores.csv", row.names = FALSE)

