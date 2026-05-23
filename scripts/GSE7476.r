
# Projeto GSE7476.R
# Guilherme Moret Miranda - Riogen
# 27/04/2026
# Download: GEOquery
# Normalização RMA com pacote affy
# Anotação: hgu133plus2.db


# ============== PACOTES ==============

library(GEOquery)
library(affy)
library(limma)
library(AnnotationDbi)
library(hgu133plus2.db)
library(hgu133plus2cdf)

# ============== EXTRAÇÃO DE DADOS ==============

id_projeto <- "GSE7476"

# carregando o master manifesto
metadata <- read.csv(metadata_path)
metadata <- metadata[metadata$study_ID == "GSE7476", c("sample_ID", "sample_type", "characteristics")]
rownames(metadata) <- metadata$sample_ID
metadata$sample_ID <- NULL

# --- Etapa de Download ---
message("\n", paste(rep("=", 30), collapse = ""))
message("Baixando dados brutos para ", id_projeto, "...")
message(paste(rep("=", 30), collapse = ""))

# Download (usando o geo_dir definido no mestre)
if(!dir.exists(file.path(geo_dir, id_projeto))) {
  getGEOSuppFiles(id_projeto, baseDir = geo_dir)
}

# definindo pasta do projeto
projeto_dir <- file.path(geo_dir, id_projeto)

# Descompactando
arquivos_tar <- list.files(
  path = projeto_dir,
  pattern = "\\.tar$",
  full.names = TRUE
)
untar(arquivos_tar, exdir = projeto_dir)

# ====== Leitura dos dados brutos ======

# --- Etapa de leitura ---
message("\n", paste(rep("=", 30), collapse = ""))
message("Lendo dados para ", id_projeto, "...")
message(paste(rep("=", 30), collapse = ""))

# Cria uma variável com arquivos que possuem "CEL.gz" e os imprime (exclue os CHP.gz)
cels.GSE7476 <- list.files(
  path = projeto_dir,
  pattern = "[cC][eE][lL]\\.gz$",
  full.names = TRUE
)

# Lendo os dados brutos (raw data)
dados_brutos_GSE7476 <- ReadAffy(filenames=cels.GSE7476)

# ====== Normalização dos dados ======

# --- Etapa de Normalização ---
message("\n", paste(rep("=", 30), collapse = ""))
message("Normalizando dados de ", id_projeto, "...")
message(paste(rep("=", 30), collapse = ""))

# normalização
dados_norm_GSE7476 <- affy::rma(dados_brutos_GSE7476)

# extração de matriz de expressão normalizada
norm_GSE7476 <- exprs(dados_norm_GSE7476)

# remove o ".CEL.gz" ou outras interferências da tabela
colnames(norm_GSE7476) <- toupper(sub(".*(GSM[0-9]+).*", "\\1", colnames(norm_GSE7476)))

# removem genes que apresentem NA
norm_corrigido_GSE7476 <- norm_GSE7476[!is.na(rownames(norm_GSE7476)), ]
norm_corrigido_GSE7476 <- norm_corrigido_GSE7476[rowSums(is.na(norm_corrigido_GSE7476)) == 0, ]
norm_corrigido_GSE7476 <- t(norm_corrigido_GSE7476)
norm_corrigido_GSE7476 <- norm_corrigido_GSE7476[,!grepl("^AFFX", colnames(norm_corrigido_GSE7476))] #remove probes de controle

# ====== remove amostras fora do manifesto ======

ids <- intersect(rownames(norm_corrigido_GSE7476), rownames(metadata))

norm_corrigido_GSE7476 <- norm_corrigido_GSE7476[ids, ]
metadata <- metadata[ids, ]

# ====== anotação com EntrezID ======

# --- Etapa de anotação ---
message("\n", paste(rep("=", 30), collapse = ""))
message("Anotando dados de ", id_projeto, "...")
message(paste(rep("=", 30), collapse = ""))


entrez_ids <- mapIds(hgu133plus2.db,
                     keys = colnames(norm_corrigido_GSE7476),
                     column = "ENTREZID",
                     keytype = "PROBEID",
                     multiVals = "first")

colnames(norm_corrigido_GSE7476) <- entrez_ids

# Remove genes NA
norm_corrigido_GSE7476 <- norm_corrigido_GSE7476[, !is.na(colnames(norm_corrigido_GSE7476))]



# Caso haja genes duplicados, é feito uma correlação de pearson.
# se a correlação for alta (0,7), é selecionado o primeiro probe.
# se a correlação for baixa, é selecionado o probe com maior variância
# se existirem probes triplicados ou mais, a comparação é realizada em clusters

# identificar grupos de probes (genes duplicados)
genes <- colnames(norm_corrigido_GSE7476)
grupos <- split(seq_along(genes), genes)

idx_final <- unlist(lapply(grupos, function(i) {
  
  # caso não haja duplicata
  if (length(i) == 1) return(i)
  
  submat <- norm_corrigido_GSE7476[, i, drop = FALSE]
  
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

# subset final
norm_corrigido_GSE7476 <- norm_corrigido_GSE7476[, idx_final]


# ====== Análise Limma ======

# --- Etapa de análise estatística ---
message("\n", paste(rep("=", 30), collapse = ""))
message("Executando análise estatística (limma) para ", id_projeto, "...")
message(paste(rep("=", 30), collapse = ""))

# cria a matriz de modelo para o limma 
fator_GSE7476 <- factor(metadata$sample_type,levels = c("non_tumor", "tumor"))
matriz_modelo_GSE7476 <- as.matrix(model.matrix(~0 + fator_GSE7476))
colnames(matriz_modelo_GSE7476) <- c('non_tumor','tumor')

# aplicação do Limma
fit <- lmFit(t(norm_corrigido_GSE7476), matriz_modelo_GSE7476)
contrast.matrix <- makeContrasts(
  Tumor_vs_NonTumor = tumor - non_tumor,
  levels = matriz_modelo_GSE7476)
fit2 <- eBayes(contrasts.fit(fit, contrast.matrix), trend = FALSE)

### extração de log2FC e Erro padrão para metafor
logFC <- fit2$coefficients[, "Tumor_vs_NonTumor"]
SE <- fit2$stdev.unscaled[, "Tumor_vs_NonTumor"] * fit2$sigma

# padronizando os LogFC para serem comparáveis entre estudos
sd_study <- sd(logFC, na.rm = TRUE)

logFC_scaled <- logFC / sd_study
logFC_scaled <- logFC_scaled - mean(logFC_scaled, na.rm = TRUE) #centra o logFC, removendo viés global do estudo

SE_scaled    <- SE / sd_study

metafor_GSE7476 <- data.frame(
  GSE7476_logFC = logFC,
  GSE7476_SE = SE,
  GSE7476_logFC_scaled = logFC_scaled,
  GSE7476_SE_scaled = SE_scaled
)

# ====== Salvar arquivo ======
## input metafor

# confere se o pData e matriz estão realmente alinhados
all(rownames(norm_corrigido_GSE7476) == rownames(metadata))

#confere rapidamente se a pasta de salvamento está pronta
out_dir <- file.path(processed_dir, id_projeto)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
rm(out_dir)

# arquivo para metafor
write.csv(
  metafor_GSE7476,
  file = file.path(processed_dir, id_projeto, "metafor_GSE7476.csv"),
  row.names = TRUE
)

# arquivo de matriz de expressão
saveRDS(norm_corrigido_GSE7476,
        file = file.path(processed_dir,
                         id_projeto,
                         "exprs_GSE7476.rds"))

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
deg_GSE7476 <- data.frame(
  ENTREZID = rownames(fit2$coefficients),
  logFC = logFC,
  p.value = p_value,
  FDR = FDR,
  significance = significance,
  stringsAsFactors = FALSE
)

# remove genes sem ENTREZID
deg_GSE7476 <- deg_GSE7476[!is.na(deg_GSE7476$ENTREZID), ]

# opcional: ordenar por FDR
deg_GSE7476 <- deg_GSE7476[order(deg_GSE7476$FDR), ]

# salvar
write.csv(
  deg_GSE7476,
  file = file.path(results_dir, "DEGs_tables", "DEGs_GSE7476.csv"),
  row.names = FALSE
)