# 
# Análise de correção de efeito de batch
# Guilherme Moret Miranda - Riogen
# 23/05/2026
# Esse script carrega dados de expressão de cada dataset obtido na análise principal
# para integrar-los e corrigir seus efeitos de batch pelo pacote sva (comBat)
# para visualizar a capacidade que as amostras tem de separar e os grupos tumor
# e não-tumor.

library(here)
library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)
library(sva)
library(readr)

# ------------------------------------------------------------
# 1. Carregar metadata
# ------------------------------------------------------------

metadata <- read_csv(
  here("metadata", "samples_metadata_bladder_cancer_TCC-2026.csv")
)

head(metadata)

# ------------------------------------------------------------
# 2. Encontrar arquivos .rds
# ------------------------------------------------------------

exprs_files <- list.files(
  path = here("data", "processed"),
  pattern = "^exprs_.*\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

# ------------------------------------------------------------
# 3. Carregar matrizes e transpor
# ------------------------------------------------------------
# Resultado:
# linhas = genes
# colunas = amostras

exprs_list <- exprs_files %>%
  set_names(
    basename(dirname(.))
  ) %>%
  map(~ {
    
    mat <- readRDS(.x)
    
    mat <- t(mat)
    
    as.matrix(mat)
  })

# ------------------------------------------------------------
# 4. Genes compartilhados
# ------------------------------------------------------------

common_genes <- reduce(
  map(exprs_list, rownames),
  intersect
)

length(common_genes)

# ------------------------------------------------------------
# 5. Filtrar genes compartilhados
# ------------------------------------------------------------

exprs_filtered <- map(
  exprs_list,
  ~ .x[common_genes, , drop = FALSE]
)

# ------------------------------------------------------------
# 6. Integrar datasets
# ------------------------------------------------------------

exprs_merged <- do.call(
  cbind,
  exprs_filtered
)

dim(exprs_merged)

# ------------------------------------------------------------
# 7. Criar vetor batch
# ------------------------------------------------------------

batch <- map2(
  exprs_filtered,
  names(exprs_filtered),
  ~ rep(.y, ncol(.x))
) %>%
  unlist()

table(batch)

# ------------------------------------------------------------
# 8. Ordenar metadata na MESMA ordem das amostras
# ------------------------------------------------------------

metadata_ordered <- metadata %>%
  filter(sample_ID %in% colnames(exprs_merged)) %>%
  slice(match(colnames(exprs_merged), sample_ID))

# Conferência crítica
all(metadata_ordered$sample_ID == colnames(exprs_merged))

# Deve retornar TRUE

# ------------------------------------------------------------
# 9. Criar variável biológica
# ------------------------------------------------------------

group <- factor(
  metadata_ordered$sample_type,
  levels = c("non_tumor", "tumor")
)

table(group)

# ------------------------------------------------------------
# 10. Modelo biológico
# ------------------------------------------------------------
# Preserva sinal tumor vs normal

mod <- model.matrix(~ group)

# ------------------------------------------------------------
# 11. ComBat
# ------------------------------------------------------------

exprs_combat <- ComBat(
  dat = exprs_merged,
  batch = batch,
  mod = mod,
  par.prior = TRUE,
  prior.plots = FALSE
)

# ------------------------------------------------------------
# 12. PCA antes
# ------------------------------------------------------------

pca_before <- prcomp(
  t(exprs_merged),
  scale. = TRUE
)

pca_before_df <- data.frame(
  PC1 = pca_before$x[,1],
  PC2 = pca_before$x[,2],
  Batch = batch,
  Group = group
)

ggplot(
  pca_before_df,
  aes(PC1, PC2, color = Group, shape = Batch)
) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA antes do ComBat"
  )

# ------------------------------------------------------------
# 13. PCA depois
# ------------------------------------------------------------

pca_after <- prcomp(
  t(exprs_combat),
  scale. = TRUE
)

pca_after_df <- data.frame(
  PC1 = pca_after$x[,1],
  PC2 = pca_after$x[,2],
  Batch = batch,
  Group = group
)

ggplot(
  pca_after_df,
  aes(PC1, PC2, color = Group, shape = Batch)
) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA após ComBat"
  )

# ------------------------------------------------------------
# 14. Salvar matriz integrada
# ------------------------------------------------------------

dir.create(
  here("data", "integrated"),
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  exprs_combat,
  here("data", "integrated", "exprs_combat_all.rds")
)