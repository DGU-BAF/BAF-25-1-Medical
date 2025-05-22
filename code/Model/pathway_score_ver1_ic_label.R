# -------------------------------------------
# 0. 환경 초기화 및 라이브러리 로드
# -------------------------------------------

# 작업 중인 모든 객체를 메모리에서 제거하여 깨끗한 환경에서 시작
rm(list = ls())

# -------------------------------------------
# 1. 데이터 불러오기 및 전처리
# -------------------------------------------

# 약물 감수성 라벨 데이터 읽기
ic_label <- read.csv("C:/Users/USER/비어플 의료/1_drug_label.csv")
# 첫 번째 열(샘플 ID)을 행 이름으로 설정하고, 이후 열 제거
rownames(ic_label) <- ic_label[,1]
ic_label <- ic_label[,-1]
# 데이터 확인 (앞 5개 샘플, 앞 5개 약물)
ic_label[1:5,1:5]

# CCLE 유전자 발현 데이터 읽기 (train/valid/test split 이후 사용)
ccle <- read.csv("C:/Users/USER/비어플 의료/#_filtered_CCLE_gene_expression.csv")
# 첫 번째 열을 행 이름(샘플 ID)으로 설정, 이후 제거
rownames(ccle) <- ccle[,1]
ccle <- ccle[,-1]
# 데이터 차원 확인 및 일부 값 확인
dim(ccle)
ccle[1:5,1:5]

# KEGG pathway 정보 읽기
kegg <- read.csv("C:/Users/USER/비어플 의료/#_filtered_Kegg_info.csv")

# TCGA 발현 데이터 읽기 및 전처리
tcga <- read.csv("C:/Users/USER/비어플 의료/final_tcga.csv")
# 두 번째 열 제거
tcga <- tcga[,-2]
# 결측 또는 문제 샘플 제거(인덱스 202,418,1941,2677)
tcga <- tcga[-c(202,418,1941,2677), ]
# 첫 번째 열을 행 이름으로 설정
rownames(tcga) <- tcga$X
# 나머지 데이터만 남기고 transpose하여 샘플×유전자 형태로 변환
tcga_t <- t(tcga)
tcga_final <- as.data.frame(tcga_t[-1, ])
# 모든 값을 숫자형으로 변환
tcga_final <- data.frame(lapply(tcga_final, as.numeric), row.names = rownames(tcga_t)[-1])
# 일부 값 확인
tcga_final[1:5,1:5]

# -------------------------------------------
# 2. 필수 라이브러리 로드
# -------------------------------------------

library(DESeq2)   # 차등 발현 분석
library(caret)    # 데이터 분할
library(dplyr)    # 데이터 조작
library(tidyr)    # 데이터 형태 변환(pivot)
library(tibble)   # rownames_to_column 등
library(stringr)  # 문자열 처리

# -------------------------------------------
# 3. 결과 저장용 리스트 초기화
# -------------------------------------------

deg_summary_list <- list()  # 나중에 약물별 DEG 요약 저장

# -------------------------------------------
# 4. 파이프라인 함수 정의: run_full_pipeline()
# -------------------------------------------

run_full_pipeline <- function(drug_name, ccle, ic_label, kegg, tcga_raw) {
  message("=== Running: ", drug_name, " ===")
  
  # 4-1) 해당 약물 라벨 추출 및 train/validation/test 분할 (7:1:2)
  ic <- ic_label[, drug_name, drop = FALSE]
  set.seed(42)
  sample_ids   <- rownames(ic)
  train_index  <- createDataPartition(ic[[1]], p = 0.7, list = FALSE)
  train_ids    <- sample_ids[train_index]
  remaining    <- setdiff(sample_ids, train_ids)
  valid_index  <- createDataPartition(ic[remaining,1], p = 1/3, list = FALSE)
  valid_ids    <- remaining[valid_index]
  test_ids     <- setdiff(remaining, valid_ids)
  
  # 분할된 ID로 CCLE 및 라벨 데이터 분리
  ccle_train <- ccle[train_ids, ]
  ccle_valid <- ccle[valid_ids, ]
  ccle_test  <- ccle[test_ids, ]
  ic_train   <- ic[train_ids, , drop = FALSE]
  ic_valid   <- ic[valid_ids, , drop = FALSE]
  ic_test    <- ic[test_ids, , drop = FALSE]
  
  # 4-2) NA 샘플 제거: 각 split에서 NA 라벨 샘플 제외
  for(split in c("train","valid","test")) {
    ic_part <- get(paste0("ic_", split))
    na_ids  <- rownames(ic_part)[is.na(ic_part[[1]])]
    if(length(na_ids) > 0) {
      assign(paste0("ic_", split), ic_part[!rownames(ic_part) %in% na_ids, , drop = FALSE])
      assign(paste0("ccle_", split), get(paste0("ccle_", split))[!rownames(get(paste0("ccle_", split))) %in% na_ids, ])
    }
  }
  
  # 4-3) DEG 분석 (DESeq2 사용)
  colData <- data.frame(
    row.names = rownames(ic_train),
    condition = factor(ic_train[[1]], levels = c("resistant","sensitive"))
  )
  dds <- DESeqDataSetFromMatrix(countData = t(ccle_train), colData = colData, design = ~ condition)
  # 저발현 유전자 필터링
  dds <- dds[rowSums(counts(dds)) > 10, ]
  # 모델 학습 및 결과 추출
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition","sensitive","resistant"))
  
  # 4-4) DEG 필터링 및 요약 정보 저장
  deg_filtered <- res[which((res$padj < 0.05) & abs(res$log2FoldChange) >= 1), ]
  deg_filtered <- deg_filtered[order(-abs(deg_filtered$log2FoldChange)), ]
  deg_top3     <- rownames(deg_filtered)[1:3]
  deg_count    <- nrow(deg_filtered)
  # 전역 리스트에 저장
  deg_summary_list[[drug_name]] <<- c(
    drug      = drug_name,
    gene1     = deg_top3[1],
    gene2     = deg_top3[2],
    gene3     = deg_top3[3],
    DEG_Count = deg_count
  )
  deg_genes <- rownames(deg_filtered)  # pathway 기준유전자 후보
  
  # 4-5) long format 변환: pathway score 계산용
  ccle_train_long <- ccle_train %>% rownames_to_column("SampleID") %>% pivot_longer(-SampleID, names_to="GeneSymbol", values_to="Expression")
  ccle_valid_long <- ccle_valid %>% rownames_to_column("SampleID") %>% pivot_longer(-SampleID, names_to="GeneSymbol", values_to="Expression")
  ccle_test_long  <- ccle_test  %>% rownames_to_column("SampleID") %>% pivot_longer(-SampleID, names_to="GeneSymbol", values_to="Expression")
  
  # 4-6) custom_pathway_score 함수 정의
  custom_pathway_score <- function(exp.df, path.genes, deg.genes, pre_ref = NULL) {
    path.exp <- exp.df %>% filter(GeneSymbol %in% path.genes)
    if(nrow(path.exp) < 3) return(NULL)
    path.wide <- path.exp %>% pivot_wider(names_from=SampleID, values_from=Expression) %>% column_to_rownames("GeneSymbol")
    path.mat <- as.matrix(path.wide)
    genes    <- rownames(path.mat)
    # 기준 유전자 선택: train은 DEG 중 분산 최대로, valid/test는 pre_ref 사용
    if(is.null(pre_ref)) {
      deg_in_path <- intersect(genes, deg.genes)
      if(length(deg_in_path) >=1) {
        deg_sd   <- apply(path.mat[deg_in_path,,drop=FALSE],1,sd)
        ref.gene <- names(which.max(deg_sd))
      } else {
        return(NULL)
      }
    } else {
      ref.gene <- pre_ref
    }
    if(!(ref.gene %in% genes)) return(NULL)
    # 상관계수 부호 계산
    safe_cor <- function(x,y){ if(sd(x,na.rm=TRUE)==0||sd(y,na.rm=TRUE)==0) return(0); cor(x,y,use="complete.obs") }
    gene.corr <- apply(path.mat,1,function(row) sign(safe_cor(path.mat[ref.gene,],row)))
    gene.corr[is.na(gene.corr)] <- 0
    # rank 기반 score 산출 및 표준화
    rank.mat   <- apply(abs(path.mat),1,rank)
    path.score <- colSums(t(rank.mat)*gene.corr)/nrow(path.mat)
    standardized <- scale(path.score)[,1]
    names(standardized) <- colnames(path.mat)
    return(list(score=standardized, ref.gene=ref.gene))
  }
  
  # 4-7) pathway별 score 계산 및 저장
  train_score_df <- data.frame(SampleID = unique(ccle_train_long$SampleID))
  valid_score_df <- data.frame(SampleID = unique(ccle_valid_long$SampleID))
  test_score_df  <- data.frame(SampleID = unique(ccle_test_long$SampleID))
  tcga_score_df  <- data.frame(SampleID = unique(tcga_final %>% rownames_to_column("SampleID") %>% pull(SampleID)))
  ref_list <- list()
  # train 기준으로 ref.gene 뽑고 score 저장
  for(pw in unique(kegg$PathwayID)){
    genes <- kegg %>% filter(PathwayID==pw) %>% mutate(NodeInfo=str_remove_all(NodeInfo,"\\[|\\]|'")) %>% separate_rows(NodeInfo,sep=",\\s*") %>% pull(NodeInfo)
    res_train <- custom_pathway_score(ccle_train_long, genes, deg_genes, pre_ref=NULL)
    if(!is.null(res_train)){
      train_score_df[[pw]] <- res_train$score
      ref_list[[pw]]       <- res_train$ref.gene
    }
  }
  # valid, test, TCGA에도 동일 ref.gene으로 score 계산
  for(pw in names(ref_list)){
    genes <- kegg %>% filter(PathwayID==pw) %>% mutate(NodeInfo=str_remove_all(NodeInfo,"\\[|\\]|'")) %>% separate_rows(NodeInfo,sep=",\\s*") %>% pull(NodeInfo)
    res_valid <- custom_pathway_score(ccle_valid_long, genes, NULL, pre_ref=ref_list[[pw]])
    if(!is.null(res_valid)) valid_score_df[[pw]] <- res_valid$score
    res_test  <- custom_pathway_score(ccle_test_long,  genes, NULL, pre_ref=ref_list[[pw]])
    if(!is.null(res_test))  test_score_df[[pw]]  <- res_test$score
    # TCGA 데이터 long format 생성 후 계산
    tcga_long <- tcga_raw %>% t() %>% as.data.frame() %>% rownames_to_column("GeneSymbol") %>% pivot_longer(-GeneSymbol,names_to="SampleID",values_to="Expression")
    res_tcga  <- custom_pathway_score(tcga_long, genes, NULL, pre_ref=ref_list[[pw]])
    if(!is.null(res_tcga)) tcga_score_df[[pw]] <- res_tcga$score[tcga_score_df$SampleID]
  }
  
  # 4-8) 결과 디렉토리 생성 및 파일 저장
  out_dir <- paste0("C:/Users/USER/비어플 의료/약물별 데이터셋/0420/", drug_name, "/")
  dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
  write.csv(train_score_df, paste0(out_dir, "train_pathway_score_", drug_name, ".csv"), row.names=FALSE)
  write.csv(valid_score_df, paste0(out_dir, "valid_pathway_score_", drug_name, ".csv"), row.names=FALSE)
  write.csv(test_score_df,  paste0(out_dir, "test_pathway_score_",  drug_name, ".csv"), row.names=FALSE)
  write.csv(ic_train,        paste0(out_dir, "ic_train_", drug_name, ".csv"), row.names=TRUE)
  write.csv(ic_valid,        paste0(out_dir, "ic_valid_", drug_name, ".csv"), row.names=TRUE)
  write.csv(ic_test,         paste0(out_dir, "ic_test_", drug_name, ".csv"), row.names=TRUE)
  write.csv(tcga_score_df,   paste0(out_dir, "tcga_pathway_score_", drug_name, ".csv"), row.names=FALSE)
}

# -------------------------------------------
# 5. 모든 약물에 대해 파이프라인 실행 및 DEG summary 저장
# -------------------------------------------

drugs <- colnames(ic_label)
for(drug in drugs) {
  run_full_pipeline(drug, ccle, ic_label, kegg, tcga_final)
}

deg_summary_df <- do.call(rbind, deg_summary_list)
rownames(deg_summary_df) <- NULL
write.csv(deg_summary_df, "C:/Users/USER/비어플 의료/약물별 데이터셋/DEG_summary.csv", row.names=FALSE)
# 최종 DEG 요약 출력
print(deg_summary_df)
