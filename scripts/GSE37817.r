
# GSE37817.R
# Guilherme Moret Miranda - Riogen
# 27/04/2026
# Download: GEOquery
# Anotação: illuminaHumanv2.db

# ============== PACOTES ==============

library(GEOquery)
library(limma)
library(AnnotationDbi)
library(R.utils)
library(illuminaHumanv2.db)

# ============== EXTRAÇÃO DE DADOS ==============

id_projeto <- "GSE37817"

# carregando o master manifesto
metadata <- read.csv(metadata_path)
metadata <- metadata[metadata$study_ID == "GSE37817", c("sample_ID", "sample_type", "characteristics")]
rownames(metadata) <- metadata$sample_ID
metadata$sample_ID <- NULL


# Criação da pasta se não existir
# Download (usando o geo_dir definido no mestre)
if(!dir.exists(file.path(geo_dir, id_projeto))) {
  getGEOSuppFiles(id_projeto, baseDir = geo_dir)
}

# ====== Leitura dos dados brutos illumina ======
# Este projeto não possui dados brutos disponíveis. A única opção é baixar a matriz normalizada do geo

# --- Etapa de Download ---
message("\n", paste(rep("=", 30), collapse = ""))
message("Baixando dados brutos para ", id_projeto, "...")
message(paste(rep("=", 30), collapse = ""))

# Baixa o GSE
gse <- getGEO('GSE37817', GSEMatrix = TRUE)

# --- Etapa de leitura ---
message("\n", paste(rep("=", 30), collapse = ""))
message("Lendo dados para ", id_projeto, "...")
message(paste(rep("=", 30), collapse = ""))

# Normalmente é a primeira plataforma
norm_GSE37817 <- exprs(gse[[1]])

# remove o ".CEL.gz" ou outras interferências da tabela
colnames(norm_GSE37817) <- toupper(sub(".*(GSM[0-9]+).*", "\\1", colnames(norm_GSE37817)))

# removem genes que apresentem NA
norm_corrigido_GSE37817 <- norm_GSE37817[!is.na(rownames(norm_GSE37817)), ]
norm_corrigido_GSE37817 <- norm_corrigido_GSE37817[rowSums(is.na(norm_corrigido_GSE37817)) == 0, ]
norm_corrigido_GSE37817 <- t(norm_corrigido_GSE37817)

# ====== remove amostras fora do manifesto ======

ids <- intersect(rownames(norm_corrigido_GSE37817), rownames(metadata))

norm_corrigido_GSE37817 <- norm_corrigido_GSE37817[ids, ]
metadata <- metadata[ids, ]

# ====== anotação com EntrezID ======

# --- Etapa de anotação ---
message("\n", paste(rep("=", 30), collapse = ""))
message("Anotando dados de ", id_projeto, "...")
message(paste(rep("=", 30), collapse = ""))


entrez_ids <- mapIds(illuminaHumanv2.db,
                     keys = colnames(norm_corrigido_GSE37817),
                     column = "ENTREZID",
                     keytype = "PROBEID",
                     multiVals = "first")

colnames(norm_corrigido_GSE37817) <- entrez_ids

# Remove genes NA
norm_corrigido_GSE37817 <- norm_corrigido_GSE37817[, !is.na(colnames(norm_corrigido_GSE37817))]

# Seleciona genes duplicados que tem maior variância
# Caso haja genes duplicados, é feito uma correlação de pearson.
# se a correlação for alta (0,7), é feita média dos sinais
# se a correlação for baixa, é selecionado o probe com maior variância
# se existirem probes triplicados ou mais, a comparação é realizada em clusters

# identificar grupos de probes (genes duplicados)
genes <- colnames(norm_corrigido_GSE37817)
grupos <- split(seq_along(genes), genes)

idx_final <- unlist(lapply(grupos, function(i) {
  
  # caso não haja duplicata
  if (length(i) == 1) return(i)
  
  submat <- norm_corrigido_GSE37817[, i, drop = FALSE]
  
  # correlação entre probes do mesmo gene
  cor_mat <- cor(submat, use = "pairwise.complete.obs")
  
  # se tudo NA → fallback variância
  if (all(is.na(cor_mat))) {
    vars <- apply(submat, 2, var, na.rm = TRUE)
    return(i[which.max(vars)])
  }
  
  # distância baseada em correlação
  dist_mat <- as.dist(1 - cor_mat)
  hc <- hclust(dist_mat, method = "average")
  
  clusters <- cutree(hc, h = 1 - 0.7)  # threshold 0.7
  
  # maior cluster
  tab <- table(clusters)
  main_cluster <- as.numeric(names(tab)[which.max(tab)])
  idx_cluster <- i[clusters == main_cluster]
  
  # se cluster confiável, usa média; senão variância
  if (length(idx_cluster) >= 2) {
    return(idx_cluster[1])  # ou poderia usar média depois
  } else {
    vars <- apply(submat, 2, var, na.rm = TRUE)
    return(i[which.max(vars)])
  }
}))

norm_corrigido_GSE37817 <- norm_corrigido_GSE37817[, idx_final]


# ====== Análise Limma ======

message("\n", paste(rep("=", 30), collapse = ""))
message("Executando análise estatística (limma) para ", id_projeto, "...")
message(paste(rep("=", 30), collapse = ""))

# cria a matriz de modelo para o limma 
fator_GSE37817 <- factor(metadata$sample_type,levels = c("non_tumor", "tumor"))
matriz_modelo_GSE37817 <- as.matrix(model.matrix(~0 + fator_GSE37817))
colnames(matriz_modelo_GSE37817) <- c('non_tumor','tumor')

# aplicação do Limma
fit <- lmFit(t(norm_corrigido_GSE37817), matriz_modelo_GSE37817)
contrast.matrix <- makeContrasts(
  Tumor_vs_NonTumor = tumor - non_tumor,
  levels = matriz_modelo_GSE37817)
fit2 <- eBayes(contrasts.fit(fit, contrast.matrix), trend = FALSE)

## calculo de SMD (Corrigido para Metodologia Clássica)

# contando o número de amostras usando a matriz de modelo
n_tumor <- sum(matriz_modelo_GSE37817[, "tumor"])
n_nontumor <- sum(matriz_modelo_GSE37817[, "non_tumor"])

# 1. Extraindo a diferença de médias (logFC)
logFC <- fit2$coefficients[, "Tumor_vs_NonTumor"]

# 2. Reconstruindo o Erro Padrão Ordinário (sem Bayes empírico)
# fit2$sigma é o desvio padrão residual clássico do modelo para cada gene
# fit2$stdev.unscaled é o fator de escala baseado na matriz de design
se_ordinary <- fit2$sigma * fit2$stdev.unscaled[, "Tumor_vs_NonTumor"]

# 3. Calculando a estatística 't' ordinária (Student's t clássico)
t_stat_ordinary <- logFC / se_ordinary

# 4. Calculando o Hedges' g (SMD) e a variância rigorosamente
es_results <- escalc(
  measure = "SMD", 
  ti = t_stat_ordinary, 
  n1i = rep(n_tumor, length(t_stat_ordinary)), 
  n2i = rep(n_nontumor, length(t_stat_ordinary))
)

# montando informações para meta-analise
metafor_GSE37817 <- data.frame(
  logFC    = logFC,
  SE_logFC = se_ordinary,
  SMD      = es_results$yi,
  SE_SMD   = sqrt(es_results$vi)
)
rownames(metafor_GSE37817) <- rownames(fit2$coefficients)

# ====== Salvar arquivo ======
## input metafor

# confere se o pData e matriz estão realmente alinhados
all(rownames(norm_corrigido_GSE37817) == rownames(metadata))

#confere rapidamente se a pasta de salvamento está pronta
out_dir <- file.path(processed_dir, id_projeto)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# arquivo para metafor
saveRDS(
  metafor_GSE37817,
  file = file.path(out_dir, "metafor_GSE37817.rds")
)

# arquivo de matriz de expressão
saveRDS(norm_corrigido_GSE37817,
        file = file.path(processed_dir,
                         id_projeto,
                         "exprs_GSE37817.rds"))

## DEGs identificados nesse projeto

# extração dos resultados do limma
logFC <- fit2$coefficients[, "Tumor_vs_NonTumor"]

# p-value bruto
p_value <- fit2$p.value[, "Tumor_vs_NonTumor"]

# ajuste FDR (Benjamini-Hochberg)
FDR <- p.adjust(p_value, method = "BH")

# classificação de significância
significance <- ifelse(
  abs(logFC) > 1 & FDR < 0.05,
  "significant",
  "not significant"
)

# dataframe final
deg_GSE37817 <- data.frame(
  ENTREZID = rownames(fit2$coefficients),
  logFC = logFC,
  p.value = p_value,
  FDR = FDR,
  significance = significance,
  stringsAsFactors = FALSE
)

# remove genes sem ENTREZID
deg_GSE37817 <- deg_GSE37817[!is.na(deg_GSE37817$ENTREZID), ]

# opcional: ordenar por FDR
deg_GSE37817 <- deg_GSE37817[order(deg_GSE37817$FDR), ]

# salvar
write.csv(
  deg_GSE37817,
  file = file.path(results_dir, "DEGs_tables", "DEGs_GSE37817.csv"),
  row.names = FALSE
)
