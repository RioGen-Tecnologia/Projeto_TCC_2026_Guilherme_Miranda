# GSE13507.R
# Guilherme Moret Miranda - Riogen
# 27/04/2026
# Download: GEOquery
# Normalização normalizeBetweenArrays com pacote limma
# Anotação: anotação fornecida nos dados

# ============== PACOTES ==============

library(GEOquery)
library(limma)
library(AnnotationDbi)
library(R.utils)
library(illuminaHumanv2.db)
library(dplyr)

# ============== EXTRAÇÃO DE DADOS ==============

id_projeto <- "GSE13507"

# carregando o master manifesto
metadata <- read.csv(metadata_path)
metadata <- metadata[metadata$study_ID == "GSE13507", c("sample_ID", "sample_type", "characteristics")]
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

# ====== Leitura dos dados brutos illumina ======

# --- Etapa de leitura ---
message("\n", paste(rep("=", 30), collapse = ""))
message("Lendo dados para ", id_projeto, "...")
message(paste(rep("=", 30), collapse = ""))

arquivos <- list.files(
  path = projeto_dir,
  pattern = ".txt$",
  full.names = TRUE
)
cels.GSE13507 <- read.delim(arquivos, stringsAsFactors = FALSE)

# extrai apenas os valores como numeric
dados_brutos_GSE13507 <- matrix(as.numeric(as.matrix(cels.GSE13507[-1, -1])),
                       nrow = nrow(cels.GSE13507)-1,
                       ncol = ncol(cels.GSE13507)-1)

# define os rownames novamente
rownames(dados_brutos_GSE13507) <- cels.GSE13507[-1, 1]
colnames(dados_brutos_GSE13507) <- toupper(sub(".*(GSM[0-9]+).*", "\\1", colnames(cels.GSE13507[-1, -1])))

# ====== Normalização dos dados ======

# --- Etapa de Normalização ---
message("\n", paste(rep("=", 30), collapse = ""))
message("Normalizando dados de ", id_projeto, "...")
message(paste(rep("=", 30), collapse = ""))

# Log2 transform
exprs_log <- log2(pmax(dados_brutos_GSE13507, 1))

# Quantile normalization
exprs_norm <- normalizeBetweenArrays(exprs_log, method = "quantile")

# removem genes que apresentem NA
norm_corrigido_GSE13507 <- exprs_norm[!is.na(rownames(exprs_norm)), ]
norm_corrigido_GSE13507 <- norm_corrigido_GSE13507[rowSums(is.na(norm_corrigido_GSE13507)) == 0, ]
norm_corrigido_GSE13507 <- t(norm_corrigido_GSE13507)

# ====== remove amostras fora do manifesto ======

ids <- intersect(rownames(norm_corrigido_GSE13507), rownames(metadata))

norm_corrigido_GSE13507 <- norm_corrigido_GSE13507[ids, ]
metadata <- metadata[ids, ]


# ====== anotação com EntrezID ======

# --- Etapa de anotação ---
message("\n", paste(rep("=", 30), collapse = ""))
message("Anotando dados de ", id_projeto, "...")
message(paste(rep("=", 30), collapse = ""))

arquivos_gunzip <- list.files(path = projeto_dir,pattern = "\\.bgx\\.gz$",full.names = TRUE)
gunzip(arquivos_gunzip,remove = FALSE,overwrite = TRUE)
arquivos_gunzip <- list.files(path = projeto_dir,pattern = "\\.bgx$",full.names = TRUE)

# Ler as primeiras linhas para inspecionar
linhas <- readLines(arquivos_gunzip, n = 100)
head(linhas, 20)  # ver onde começa "ProbeID" ou "ID"

# Encontrar a linha que contém o cabeçalho
linha_cabecalho <- grep("^Species", linhas)
linha_cabecalho

anotação_ref <- read.delim(arquivos_gunzip,
                           skip = linha_cabecalho - 1,
                           stringsAsFactors = FALSE)


# 1. criar tabela de anotação com apenas Probe_Id e EntrezID
anot <- anotação_ref %>%
  select(Probe_Id, Entrez_Gene_ID) %>%
  filter(!is.na(Entrez_Gene_ID))  # remove probes sem EntrezID

# 2. filtrar a matriz para conter apenas probes com EntrezID
expr_filtrada <- norm_corrigido_GSE13507[, colnames(norm_corrigido_GSE13507) %in% anot$Probe_Id]

# 3. garantir que a ordem dos probes da anotação bate com a matriz
anot <- anot[match(colnames(expr_filtrada), anot$Probe_Id), ]

# 4. substituir colnames pelos EntrezIDs
colnames(expr_filtrada) <- anot$Entrez_Gene_ID

# Remove genes NA
norm_corrigido_GSE13507 <- expr_filtrada[, !is.na(colnames(expr_filtrada))]
rm(expr_filtrada)

# Seleciona genes duplicados que tem maior variância
# Caso haja genes duplicados, é feito uma correlação de pearson.
# se a correlação for alta (0,7), é feita média dos sinais
# se a correlação for baixa, é selecionado o probe com maior variância
# se existirem probes triplicados ou mais, a comparação é realizada em clusters

# identificar grupos de probes (genes duplicados)
genes <- colnames(norm_corrigido_GSE13507)
grupos <- split(seq_along(genes), genes)

idx_final <- unlist(lapply(grupos, function(i) {
  
  # caso não haja duplicata
  if (length(i) == 1) return(i)
  
  submat <- norm_corrigido_GSE13507[, i, drop = FALSE]
  
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

norm_corrigido_GSE13507 <- norm_corrigido_GSE13507[, idx_final]


# ====== Análise Limma ======

# --- Etapa de análise estatística ---
message("\n", paste(rep("=", 30), collapse = ""))
message("Executando análise estatística (limma) para ", id_projeto, "...")
message(paste(rep("=", 30), collapse = ""))

# cria a matriz de modelo para o limma 
fator_GSE13507 <- factor(metadata$sample_type,levels = c("non_tumor", "tumor"))
matriz_modelo_GSE13507 <- as.matrix(model.matrix(~0 + fator_GSE13507))
colnames(matriz_modelo_GSE13507) <- c('non_tumor','tumor')

# aplicação do Limma
fit <- lmFit(t(norm_corrigido_GSE13507), matriz_modelo_GSE13507)
contrast.matrix <- makeContrasts(
  Tumor_vs_NonTumor = tumor - non_tumor,
  levels = matriz_modelo_GSE13507)
fit2 <- eBayes(contrasts.fit(fit, contrast.matrix), trend = FALSE)

### extração de log2FC e Erro padrão para metafor
logFC <- fit2$coefficients[, "Tumor_vs_NonTumor"]
SE <- fit2$stdev.unscaled[, "Tumor_vs_NonTumor"] * fit2$sigma

# padronizando os LogFC para serem comparáveis entre estudos
sd_study <- sd(logFC, na.rm = TRUE)

logFC_scaled <- logFC / sd_study
logFC_scaled <- logFC_scaled - mean(logFC_scaled, na.rm = TRUE) #centra o logFC, removendo viés global do estudo

SE_scaled    <- SE / sd_study

metafor_GSE13507 <- data.frame(
  GSE13507_logFC = logFC,
  GSE13507_SE = SE,
  GSE13507_logFC_scaled = logFC_scaled,
  GSE13507_SE_scaled = SE_scaled
)


# ====== Salvar arquivo ======

# confere se o pData e matriz estão realmente alinhados
all(rownames(norm_corrigido_GSE13507) == rownames(metadata))

#confere rapidamente se a pasta de salvamento está pronta
out_dir <- file.path(processed_dir, id_projeto)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
rm(out_dir)

# arquivo para metafor
write.csv(
  metafor_GSE13507,
  file = file.path(processed_dir, id_projeto, "metafor_GSE13507.csv"),
  row.names = TRUE
)