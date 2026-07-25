# Avaliação de controle de qualidade por integração por comBat e gráficos de PCA
# Guilherme Moret Miranda - Riogen
# 14/07/2026
# este script carrega as matrizes de expressão analisadas, as integra por
# comBat (sva) e disponibiliza gráficos PCA com o objetivo de análise de controle
# de qualidade.

# ============== CARREGANDO PACOTES ==============

message("================================")
message("CARREGANDO PACOTES")
message("================================")

library(here)
library(sva)
library(ggplot2)
library(plotly)

# ============== DEFINIÇÃO DE DIRETÓRIO E ARQUIVOS ==============
# Definir a pasta de trabalho e o master manifesto, com anotações das amostras

geo_dir      <- here("data", "raw", "GEO")
processed_dir <- here("data", "processed")
results_dir  <- here("results")
figures_dir  <- here("figures")
scripts_dir  <- here("scripts")
metadata_dir <- here("metadata")

dirs <- c(
  geo_dir,
  processed_dir,
  results_dir,
  figures_dir,
  scripts_dir
)

# cria os diretórios se eles não existirem
for (d in dirs) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
  }
}

# Master Manifesto com anotações das amostras
metadata <- read.csv(file.path(metadata_dir,"batch_corrected_samples_metadata_bladder_cancer_TCC-2026.csv"))

# limpando
rm(d,dirs)
gc()

# ============== CARREGAMENTO PROJETOS ==============
# Carrega as matrizes de expressão geradas na análise principal na forma de
# arquivos rds e os integra para integração combat

message("================================")
message("CARREGANDO DATASETS")
message("================================")

# Lista das pastas GSE
gse_dirs <- list.dirs(processed_dir, recursive = FALSE)

# Nome dos datasets
gse_names <- basename(gse_dirs)

# Carregar todas as matrizes .rds
matrix_list <- lapply(gse_dirs, function(dir) {
  rds_file <- list.files(
    dir,
    pattern = "exprs_.*\\.rds$",
    full.names = TRUE
  )
  readRDS(rds_file)
})

# Nomear os elementos da lista
names(matrix_list) <- gse_names



# ============== INTEGRAÇÃO DE MATRIZES ==============
# Criando uma matriz unificada com genes comuns para correção de batch

message("================================")
message("INTEGRANDO DATASETS")
message("================================")

# genes comuns entre todos os estudos
genes_comuns <- Reduce(
  intersect,
  lapply(matrix_list, colnames)
)

message(
  "Número de genes comuns: ",
  length(genes_comuns)
)

# padronizando genes
matrix_list_common <- lapply(
  matrix_list,
  function(x) x[, genes_comuns, drop = FALSE]
)

# integrando amostras
expr <- do.call(
  rbind,
  matrix_list_common
)

# manter somente amostras presentes no metadata
common_samples <- intersect(
  rownames(expr),
  metadata$sample_ID
)

expr <- expr[common_samples, ]

metadata <- metadata[
  match(common_samples, metadata$sample_ID),
]


# ============== CORREÇÃO DE BATCHES ==============

message("================================")
message("CORRIGINDO BATCH EFFECT")
message("================================")


# batch = estudo GEO
batch <- metadata$study_ID


# grupos biológicos
group <- metadata$sample_type


# características clínicas
characteristics <- metadata$characteristics

stopifnot(
  all(rownames(expr) == metadata$sample_ID)
)

# variável para visualização
sample_class <- characteristics

sample_class <- factor(
  sample_class,
  levels = c(
    "MIBC",
    "NMIBC",
    "non_cancer_individual"
  )
)


# modelo preservando grupos biológicos
mod <- model.matrix(~ group)


# ComBat
expr_combat <- ComBat(
  dat   = t(expr),
  batch = batch,
  mod   = mod
)


# voltar para formato amostras x genes
expr_combat <- t(expr_combat)


# ============== GRÁFICOS PCA ==============
# Analisando o efeito e eficiência de correção de batches


## PCA antes

pca_before <- prcomp(expr, center = TRUE, scale. = TRUE)

# Data frame para o gráfico
pca_before_df <- data.frame(
  PC1 = pca_before$x[,1],
  PC2 = pca_before$x[,2],
  Batch = batch,
  Sample_Class = sample_class,
  Sample = rownames(expr)
)

# Variância explicada
var_before <- summary(pca_before)$importance[2, ] * 100

png(file.path(figures_dir,"PCA_pré_combat.png"),height = 2500,width = 2400,res = 300)
ggplot(
  pca_before_df,
  aes(
    PC1,
    PC2,
    color = Batch,
    shape = Sample_Class
  )
) +
  geom_point(size = 3.5, alpha = 0.9) +
  scale_shape_manual(
    values = c(
      MIBC = 1,
      NMIBC = 2,
      non_cancer_individual = 3
    )
  ) +
  
  labs(
    title = "PCA antes do ComBat",
    x = paste0("PC1 (", round(var_before[1],1), "%)"),
    y = paste0("PC2 (", round(var_before[2],1), "%)"),
    color = "Batch",
    shape = "Tipo de amostra"
  ) +
  
  theme_classic(base_size = 13)
dev.off()

## PCA depois

pca_after <- prcomp(
  expr_combat,
  center = TRUE,
  scale. = TRUE
)

pca_after_df <- data.frame(
  PC1 = pca_after$x[,1],
  PC2 = pca_after$x[,2],
  Batch = batch,
  Sample_Class = sample_class,
  Sample = rownames(expr_combat)
)

# grupo para elipse
pca_after_df$Tumor_Status <- ifelse(
  pca_after_df$Sample_Class %in% c("MIBC", "NMIBC"),
  "Tumor",
  "Non-tumor"
)


var_after <- summary(pca_after)$importance[2, ] * 100

png(file.path(figures_dir,"PCA_pós_combat.png"),height = 2500,width = 2400,res = 300)
ggplot(
  pca_after_df,
  aes(
    PC1,
    PC2,
    color = Batch,
    shape = Sample_Class
  )
) +
  stat_ellipse(
    data = pca_after_df,
    aes(
      x = PC1,
      y = PC2,
      group = Tumor_Status
    ),
    level = 0.60,
    linewidth = 1,
    color = "black"
  ) +
  geom_point(
    size = 3.5,
    alpha = 0.9
  ) +
  scale_shape_manual(
    values = c(
      MIBC = 1,
      NMIBC = 2,
      non_cancer_individual = 3
    )
  ) +
  labs(
    title = "PCA após ComBat",
    x = paste0("PC1 (", round(var_after[1],1), "%)"),
    y = paste0("PC2 (", round(var_after[2],1), "%)"),
    color = "Batch",
    shape = "Tipo de amostra"
  ) +
  theme_classic(base_size = 13)
dev.off()

## PCA pós combat 3D

## PCA 3D


# Data frame para o PCA 3D
pca_3d_df <- data.frame(
  PC1 = pca_after$x[, 1],
  PC2 = pca_after$x[, 2],
  PC3 = pca_after$x[, 3],
  Batch = batch,
  Grupo = group,
  GSM = rownames(expr_combat)
)

plot_ly(
  pca_3d_df,
  x = ~PC1,
  y = ~PC2,
  z = ~PC3,
  color = ~Batch,
  symbol = ~Grupo,
  symbols = c("circle-open", "cross"),
  text = ~paste(
    "<b>Amostra:</b>", GSM,
    "<br><b>Dataset:</b>", Batch,
    "<br><b>Grupo:</b>", Grupo
  ),
  hoverinfo = "text",
  type = "scatter3d",
  mode = "markers",
  marker = list(
    size = 5,
    opacity = 0.8
  )
) %>%
  layout(
    title = "PCA 3D após correção por ComBat",
    scene = list(
      xaxis = list(title = paste0("PC1 (", round(var_after[1], 1), "%)")),
      yaxis = list(title = paste0("PC2 (", round(var_after[2], 1), "%)")),
      zaxis = list(title = paste0("PC3 (", round(var_after[3], 1), "%)"))
    ),
    legend = list(
      title = list(text = "Dataset GEO")
    )
  )
