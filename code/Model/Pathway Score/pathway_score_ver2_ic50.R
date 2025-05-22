# -------------------------------------------
# 0.라이브러리 로드
# -------------------------------------------

rm(list = ls()) 

library(DESeq2)
library(caret)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)

# -------------------------------------------
# 1. IC50 데이터 전처리 및 라벨 생성
# -------------------------------------------

# IC50 데이터 불러오기
ic50 <- read.csv("C:/Users/User/BAF-의료/data/ic50_reshape.csv")
rownames(ic50) <- ic50$ARXSPAN_ID  # 행 이름을 샘플 ID로 설정
ic50 <- ic50[ , -1]                # ARXSPAN_ID 열 제거

# pIC50 계산: IC50 값에 -log10을 취해 감수성 수치화 (클수록 더 sensitive)
pic50 <- log10(ic50) * (-1)
pic50_group <- pic50  # 그룹 정보를 저장할 데이터프레임 복사

# 각 약물별로 z-score 기반 그룹 분류: 
#   z > 0 → "sensitive", z ≤ 0 → "resistant"
for (drug in colnames(pic50)) {
  values   <- pic50[[drug]]
  z_values <- scale(values)  # 평균 0, 표준편차 1로 정규화
  group    <- ifelse(
    is.na(z_values), NA,
    ifelse(z_values > 0, "sensitive", "resistant")
  )
  pic50_group[[drug]] <- group
}

# DEG 요약 정보 불러오기
# 사전에 모든 약물에 DEG를 수행해 각 약물별 기준 유전자의 기준을 만족하는 유전자들의 수가 담겨있음
# 기준유전자의 기준은 민감군과 저항군 차이의 p-value가 0.05 이하이며 logfoldchange가 2 이상이어야함함
deg_df    <- read.csv("C:/Users/USER/Downloads/DEG_summary_0410.csv")

# DEG 개수와 NA 개수 기준으로 top10 약물 선택
deg_select <- deg_df[
  (deg_df$DEG_Count > 100) &
    (deg_df$NA_Count  < 200), ]
drug_names <- head(
  deg_select[order(deg_select$DEG_Count, decreasing = TRUE), ],
  10
)$drug

# 선택된 약물들의 라벨 추출
pic_label <- pic50_group[, drug_names]
ic_label  <- ic50[         , drug_names]
drug_names
## 선택한 약물들 엑셀로 내보내기
# 예: drug_names 에 필터링할 약물 이름 벡터가 들어 있다고 가정
# drug_names <- c("ABT.199", "AFATINIB", "ACY.1215", ...)

filtered_df <- deg_df[ deg_df$drug %in% drug_names, ]
print(filtered_df)
# 작업 디렉토리에 "filtered_deg.csv" 로 저장
write.csv(filtered_df,
          file = "C:/Users/User/BAF-의료/data/모델링/약물별 데이터셋_final_0420/ic50에서 파생/filtered_deg.csv",
          row.names = FALSE,    # 행 이름(인덱스) 제외
          fileEncoding = "UTF-8")  # 한글 깨짐 방지


# 각 약물별 민감도/저항성 분포 확인
for (i in seq_along(drug_names)) {
  cat("==", drug_names[i], "==\n")
  print(table(pic_label[, i]))
}

# -------------------------------------------
# 2. CCLE, KEGG, TCGA 데이터 불러오기 및 전처리
# -------------------------------------------

# CCLE 유전자 발현 데이터
ccle <- read.csv("C:/Users/USER/비어플 의료/#_filtered_CCLE_gene_expression.csv")
rownames(ccle) <- ccle[, 1]
ccle <- ccle[ , -1]  # 첫 열 제거

# KEGG pathway 정보
kegg <- read.csv("C:/Users/USER/비어플 의료/#_filtered_Kegg_info.csv")

# TCGA 발현 데이터
tcga <- read.csv("C:/Users/USER/비어플 의료/final_tcga.csv")
tcga <- tcga[ , -2]                     # 불필요 열 제거
tcga <- tcga[-c(202,418,1941,2677), ]    # 중복되는 샘플 제거
rownames(tcga) <- tcga$X
tcga_t <- t(tcga)
tcga_final <- as.data.frame(tcga_t[-1, ])
tcga_final <- data.frame(
  lapply(tcga_final, as.numeric),
  row.names = rownames(tcga_t)[-1]
)

# -------------------------------------------
# 3. custom_pathway_score 함수 
# -------------------------------------------

# 특정 pathway 내 유전자들의 발현 데이터와 DEG 결과를 이용해
# 샘플별 pathway score를 계산하는 함수 정의
custom_pathway_score <- function(exp.df, path.genes, deg.genes, pre_ref = NULL) {
  # 1) pathway 내 유전자만 필터링
  path.exp <- exp.df %>% filter(GeneSymbol %in% path.genes)
  if (nrow(path.exp) < 3) return(NULL)  # 너무 적으면 skip
  
  # 2) wide-format (GeneSymbol × SampleID) 행렬 생성
  path.wide <- path.exp %>%
    pivot_wider(names_from = SampleID, values_from = Expression) %>%
    column_to_rownames("GeneSymbol")
  path.mat <- as.matrix(path.wide)
  genes <- rownames(path.mat)
  
  # 3) 기준 유전자(ref.gene) 선택
  if (is.null(pre_ref)) {
    # 3-1) Train 데이터에서는는 사전 지정된 ref.gene이 없으니까까 
    # deg.genes 중 pathway에 포함된 유전자를 후보로 삼음
    deg_in_path <- intersect(genes, deg.genes)
    
    if (length(deg_in_path) >= 1) {
      # 3-2) 후보군 중 발현 분산(sd)이 가장 큰 유전자를 ref.gene으로 선택
      deg_sd   <- apply(path.mat[deg_in_path, , drop = FALSE], 1, sd)
      ref.gene <- names(which.max(deg_sd))
    } else {
      # 3-3) pathway 내에 DEG가 하나도 없으면 더 이상 score 계산 불가 → NULL 반환 →pathway 제외
      return(NULL)
    }
  } else {
    # 3-4) valid, test, TCGA에서는 ref.gene 그대로 사용
    ref.gene <- pre_ref
  }
  
  # 3-5) 선택된 ref.gene이 실제 pathway 유전자 리스트에 없으면 계산 중단
  if (!(ref.gene %in% genes)) return(NULL)
  
  
  # 4) 안전한 상관계수 함수 정의
  safe_cor <- function(x, y) {
    if (length(x) < 2 || length(y) < 2 ||
        sd(x, na.rm = TRUE) == 0 || sd(y, na.rm = TRUE) == 0) {
      return(0)
    }
    cor(x, y, use = "complete.obs")
  }
  
  # 5) 각 유전자별로 ref.gene과의 상관계수 부호만 추출
  gene.corr <- apply(path.mat, 1, function(row) {
    sign(safe_cor(path.mat[ref.gene, ], row))
  })
  gene.corr[is.na(gene.corr)] <- 0
  
  # 6) 절댓값 발현에 대해 rank 매기기
  rank.mat <- apply(abs(path.mat), 1, rank)
  
  # 7) sign * rank 합산 후 pathway size로 나눠 평균화 → score
  path.score <- colSums(t(rank.mat) * gene.corr) / nrow(path.mat)
  
  # 8) z-score 표준화
  standardized <- scale(path.score)[, 1]
  names(standardized) <- colnames(path.mat)
  
  return(list(score = standardized, ref.gene = ref.gene))
}

# -------------------------------------------
# 4. 메인 파이프라인: train/valid/test 분할 → DEG → pathway score → 저장
# -------------------------------------------

run_full_pipeline <- function(drug_name, ccle, ic_label, kegg, tcga_raw) {
  message("=== Running: ", drug_name, " ===")
  ic <- ic_label[, drug_name, drop = FALSE]
  
  # (1) train/valid/test split (7:1:2)
  set.seed(42)
  ids            <- rownames(ic)
  train_idx      <- createDataPartition(ic[[1]], p = 0.7, list = FALSE)
  train_ids      <- ids[train_idx]
  remaining_ids  <- setdiff(ids, train_ids)
  valid_idx      <- createDataPartition(ic[remaining_ids, 1], p = 1/3, list = FALSE)
  valid_ids      <- remaining_ids[valid_idx]
  test_ids       <- setdiff(remaining_ids, valid_ids)
  
  ccle_train <- ccle[train_ids, ]
  ccle_valid <- ccle[valid_ids, ]
  ccle_test  <- ccle[test_ids, ]
  
  ic_train <- ic[train_ids, , drop = FALSE]
  ic_valid <- ic[valid_ids, , drop = FALSE]
  ic_test  <- ic[test_ids, , drop = FALSE]
  
  # (2) NA 샘플 제거
  for (split in c("train", "valid", "test")) {
    ic_part <- get(paste0("ic_", split))
    na_ids  <- rownames(ic_part)[is.na(ic_part[[1]])]
    if (length(na_ids) > 0) {
      assign(paste0("ic_", split),
             ic_part[!rownames(ic_part) %in% na_ids, , drop = FALSE])
      assign(paste0("ccle_", split),
             get(paste0("ccle_", split))[!rownames(get(paste0("ccle_", split))) %in% na_ids, ])
    }
  }
  
  # (3) DEG 분석: resistant vs sensitive 그룹 간 발현 차이 유전자 탐색
  # --------------------------------------------
  # 1) 샘플별 그룹 정보 메타데이터 생성
  colData <- data.frame(
    row.names = rownames(ic_train),                        
    condition  = factor(                                   
      ic_train[[1]],                                        
      levels = c("resistant", "sensitive")                  
    )
  )
  
  # 2) DESeq2용 객체 생성: gene × sample 형태로 transpose
  dds <- DESeqDataSetFromMatrix(
    countData = t(ccle_train),  
    colData   = colData,        
    design    = ~ condition     
  )
  
  # 3) 노이즈 제거: 총 read count가 10 이하인 유전자는 제거
  dds <- dds[rowSums(counts(dds)) > 10, ]
  
  # 4) 모델 학습: dispersion 추정→모델 피팅→통계 검정 일괄 실행
  dds <- DESeq(dds)
  
  # 5) 그룹 간 비교 결과 추출: sensitive 대 resistant
  res <- results(
    dds,
    contrast = c("condition", "sensitive", "resistant")
  )
  
  # 6) FDR(padj) 기준 및 fold change 기준으로 유의미 DEG 필터링
  deg_filtered <- res[
    (res$padj < 0.05) & (abs(res$log2FoldChange) >= 1),
  ]
  
  # 7) 발현 변화량 크기 기준(절댓값 내림차순) 정렬
  deg_filtered <- deg_filtered[order(-abs(deg_filtered$log2FoldChange)), ]
  
  # 8) top3 DEG 유전자 이름 추출 (후속 분석에 사용)
  deg_top3 <- rownames(deg_filtered)[1:3]
  
  # 9) 전체 DEG 개수 저장 (DEG summary 테이블용)
  deg_count <- nrow(deg_filtered)
  
  # 10) 전역 리스트에 요약 결과 기록
  deg_summary_list[[drug_name]] <<- c(
    drug      = drug_name,
    gene1     = deg_top3[1],
    gene2     = deg_top3[2],
    gene3     = deg_top3[3],
    DEG_Count = deg_count
  )
  
  # 11) pathway scoring에 사용할 DEG 목록으로 저장
  deg_genes <- rownames(deg_filtered)
  
  # (4) long format 변환
  ccle_train_long <- ccle_train %>%
    rownames_to_column("SampleID") %>%
    pivot_longer(-SampleID, names_to = "GeneSymbol", values_to = "Expression")
  ccle_valid_long <- ccle_valid %>%
    rownames_to_column("SampleID") %>%
    pivot_longer(-SampleID, names_to = "GeneSymbol", values_to = "Expression")
  ccle_test_long  <- ccle_test %>%
    rownames_to_column("SampleID") %>%
    pivot_longer(-SampleID, names_to = "GeneSymbol", values_to = "Expression")
  
  # (5) pathway score 계산 (train 기준)
  train_score_df <- data.frame(SampleID = unique(ccle_train_long$SampleID))
  valid_score_df <- data.frame(SampleID = unique(ccle_valid_long$SampleID))
  test_score_df  <- data.frame(SampleID = unique(ccle_test_long$SampleID))
  ref_list <- list()
  
  for (pw in unique(kegg$PathwayID)) {
    genes <- kegg %>%
      filter(PathwayID == pw) %>%
      mutate(NodeInfo = str_remove_all(NodeInfo, "\\[|\\]|'")) %>%
      separate_rows(NodeInfo, sep = ",\\s*") %>%
      pull(NodeInfo)
    
    result <- custom_pathway_score(ccle_train_long, genes, deg_genes, pre_ref = NULL)
    if (!is.null(result)) {
      train_score_df[[pw]] <- result$score
      ref_list[[pw]]      <- result$ref.gene
    }
  }
  
  # (6) valid/test/TCGA에도 같은 ref.gene 사용해 score 계산
  tcga_long <- tcga_raw %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("GeneSymbol") %>%
    pivot_longer(-GeneSymbol, names_to = "SampleID", values_to = "Expression")
  tcga_score_df <- data.frame(SampleID = unique(tcga_long$SampleID))
  
  for (pw in names(ref_list)) {
    genes <- kegg %>%
      filter(PathwayID == pw) %>%
      mutate(NodeInfo = str_remove_all(NodeInfo, "\\[|\\]|'")) %>%
      separate_rows(NodeInfo, sep = ",\\s*") %>%
      pull(NodeInfo)
    
    # TCGA
    result_tcga <- custom_pathway_score(tcga_long, genes, NULL, pre_ref = ref_list[[pw]])
    if (!is.null(result_tcga)) {
      scores <- result_tcga$score
      tcga_score_df[[pw]] <- scores[tcga_score_df$SampleID]
    }
    # Validation
    res_valid <- custom_pathway_score(ccle_valid_long, genes, NULL, pre_ref = ref_list[[pw]])
    if (!is.null(res_valid)) valid_score_df[[pw]] <- res_valid$score
    # Test
    res_test  <- custom_pathway_score(ccle_test_long, genes, NULL, pre_ref = ref_list[[pw]])
    if (!is.null(res_test))  test_score_df[[pw]]  <- res_test$score
  }
  
  # (7) 결과 디렉토리 생성 및 저장
  out_dir <- paste0("C:/Users/USER/비어플 의료/0420/", drug_name, "/")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  write.csv(train_score_df,
            paste0(out_dir, "train_pathway_score_", drug_name, ".csv"),
            row.names = FALSE)
  write.csv(valid_score_df,
            paste0(out_dir, "valid_pathway_score_", drug_name, ".csv"),
            row.names = FALSE)
  write.csv(test_score_df,
            paste0(out_dir, "test_pathway_score_", drug_name, ".csv"),
            row.names = FALSE)
  
  write.csv(ic_train,
            paste0(out_dir, "ic_train_", drug_name, ".csv"),
            row.names = TRUE)
  write.csv(ic_valid,
            paste0(out_dir, "ic_valid_", drug_name, ".csv"),
            row.names = TRUE)
  write.csv(ic_test,
            paste0(out_dir, "ic_test_", drug_name, ".csv"),
            row.names = TRUE)
  
  write.csv(tcga_score_df,
            paste0(out_dir, "tcga_pathway_score_", drug_name, ".csv"),
            row.names = FALSE)
}

# -------------------------------------------
# 5. 모든 약물에 대해 파이프라인 실행 및 DEG summary 저장
# -------------------------------------------

deg_summary_list <- list()
for (drug in colnames(pic_label)) {
  run_full_pipeline(drug, ccle, pic_label, kegg, tcga_final)
}

deg_summary_df <- do.call(rbind, deg_summary_list)
rownames(deg_summary_df) <- NULL
write.csv(deg_summary_df,
          "C:/Users/USER/비어플 의료/약물별 데이터셋/DEG_summary.csv",
          row.names = FALSE)

# 최종 결과 확인
print(deg_summary_df)
