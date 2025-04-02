# 환경 정리
rm(list = ls()) 

# 유전자 발현 데이터 로드
df <- read.csv("C:/Users/User/BAF-의료/data/CCLE_gene_correspond.csv", 
               header = TRUE, row.names = 1, check.names = FALSE)
rcm <- df
rcm <- t(rcm)  # 유전자와 세포주 transpose

# 유전자 이름 정리 (괄호 및 괄호 안 숫자 제거)
rownames(rcm) <- gsub(" \\(.*", "", rownames(rcm))
rownames(rcm) <- gsub("\\(", "-", rownames(rcm))

# read count를 정수로 변환
rcm <- round(rcm, digits = 0)

# 데이터 크기 확인
dim(rcm)

# 발현량 분포 확인
hist(rcm[,1], breaks = 100)
hist(log2(rcm[,1] + 1), breaks = 50)

# 무작위로 9개 세포주 선택하여 분포 시각화
ri <- sample(1:ncol(rcm), 9)
par(mfrow = c(3, 3))
for(i in ri){
  hist(log2(rcm[, i]), breaks = 30, main = colnames(rcm)[i], xlab = "Read count")
}

# 유전자별 표준편차 분포 확인
rcm_sd <- apply(log2(rcm), 1, sd)
hist(rcm_sd, breaks = 20, col = "skyblue", xlab = "SD of read count", main = "")

# edgeR 패키지 사용한 유전자 필터링
library(edgeR)

# 모든 샘플을 하나의 그룹으로 설정
group <- factor(rep("all", ncol(rcm)))

# filterByExpr 함수로 low expression 유전자 필터링
keep <- filterByExpr(rcm, group = group)
rcm_filtered <- rcm[keep, ]
dim(rcm_filtered)

# 전치하여 세포주가 행, 유전자가 열이 되도록
rcm_filtered.t <- t(rcm_filtered)

# 결과 저장
write.csv(rcm_filtered.t, "C:/Users/User/BAF-의료/data/CCLE_gene_EdgeR.csv", row.names = TRUE)
