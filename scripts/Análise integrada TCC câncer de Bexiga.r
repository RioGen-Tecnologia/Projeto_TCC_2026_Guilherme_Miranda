
# Análise integrada de projetos de microarray obtidas do Gene Expression Omnibus
# Guilherme Moret Miranda - Riogen
# 25/04/2026
# Esse script obtém, normaliza e análisa cada um dos projetos buscados previamente,
# combina os dados de cada projeto por meta-análise e realiza posteriores análises
# para determinar biomarcadores gênicos de câncer de bexiga.

# ============== CARREGANDO PACOTES ==============

source(here("scripts", "setup.R"))

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
metadata_path <- list.files(
  here("metadata"),
  pattern = "\\.csv$",
  full.names = TRUE
)

# ============== ANÁLISE DE PROJETOS ==============
# Esta sessão executa os scripts secundários que extrem os dados do GEO, normalizam,
# anotam e análisam estatísticamente por limma cada um dos projetos.


# Lista scripts que começam com "GSE"
scripts_projetos <- list.files(
  path = scripts_dir,
  pattern = "^GSE.*\\.[Rr]$",
  full.names = TRUE
)

# Loop de execução dos scripts secundários
arquivos_metafor <- list()

for (script in scripts_projetos) {
  message("Processando: ", basename(script))
  env <- new.env()
  source(script, local = env)
  id_atual <- gsub(".*(GSE[0-9]+).*", "\\1", basename(script))
  nome_objeto <- paste0("metafor_", id_atual)
  if (exists(nome_objeto, envir = env)) {
    arquivos_metafor[[id_atual]] <- get(nome_objeto, envir = env)
  } else {
    warning(nome_objeto, " não encontrado.")
  }
  rm(env)
  gc()
}


# ============== PREPARAÇÃO PARA META-ANÁLISE ==============
# Nesta seção os dados dos diferentes projetos são integrados.
# O input da meta-análise gerado foi os valores de Log2FC e Erro padrão para cada,
# gene identificado em cada projeto.
# O erro padrão das estimativas foi calculado como o produto entre o desvio padrão
# residual moderado e um fator dependente da matriz de design:
# SE = fit2$stdev.unscaled * fit2$sigma
# SE = estrutura do modelo * ruído dos dados

# 1. Preparar a lista convertendo rownames em uma coluna temporária
# Isso evita que o R perca os IDs durante o merge
entrez <- lapply(arquivos_metafor, function(df) {
  df$EntrezID <- rownames(df)
  return(df)
})

# 2. Unificar todos os dataframes pelo EntrezID
input_metafor <- Reduce(function(x, y) merge(x, y, by = "EntrezID", all = TRUE), entrez)

# 3. Transformar o EntrezID de volta em rownames e remover a coluna
rownames(input_metafor) <- input_metafor$EntrezID
input_metafor$EntrezID <- NULL


# ============== META-ANÁLISE ==============
# Salvamos o input da meta-análise e executamos a meta-análise devidamente utilizando
# o pacote metafor.
# Foi utilizado como critério genes que aparecem em no mínimo 3 projetos simultaneamente.

# confeindo rapidamente se a pasta de salvamento está pronta
out_dir <- file.path(results_dir, "meta_analysis")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
rm(out_dir)

#  salvando tabela de input do metafor
write.csv(input_metafor,
          file = file.path(results_dir,
                           "meta_analysis",
                           "metafor_GSE3167.csv"),
          row.names = TRUE)


# obtenção das colunas de efeito e erro padrão
logFC_cols <- grep("_logFC_scaled$", colnames(input_metafor))
SE_cols    <- grep("_SE_scaled$", colnames(input_metafor))

logFC_orig_cols <- grep("_logFC$", colnames(input_metafor))
SE_orig_cols    <- grep("_SE$", colnames(input_metafor))

# Filtrar genes presentes em ≥3 estudos
n_studies <- rowSums(!is.na(input_metafor[, logFC_cols]))
metafor_filtered <- input_metafor[n_studies >= 3, ]

# Função que roda meta-análise para um gene
meta_gene <- function(i){
  # dados escalados (para meta-análise)
  yi  <- as.numeric(metafor_filtered[i, logFC_cols])
  sei <- as.numeric(metafor_filtered[i, SE_cols])
  keep <- !is.na(yi) & !is.na(sei)
  # exigir pelo menos 3 estudos (mais robusto)
  if (sum(keep) < 2) {
    return(c(NA, NA, NA, NA, NA))
  }
  out <- tryCatch({
    # meta-análise no espaço padronizado
    fit <- rma(yi = yi[keep],
               sei = sei[keep],
               method = "REML",
               test = "knha")
    # reconstrução do logFC original
    yi_orig  <- as.numeric(metafor_filtered[i, logFC_orig_cols])
    sei_orig <- as.numeric(metafor_filtered[i, SE_orig_cols])
    
    # 🔥 usar exatamente os mesmos estudos do modelo
    yi_orig  <- yi_orig[keep]
    sei_orig <- sei_orig[keep]
    # evitar divisão por zero
    sei_orig[sei_orig == 0] <- NA
    valid <- !is.na(yi_orig) & !is.na(sei_orig)
    if (sum(valid) < 2) {
      logFC_meta <- NA
    } else {
      w <- 1 / (sei_orig[valid]^2)
      logFC_meta <- sum(yi_orig[valid] * w) / sum(w)
    }
    c(logFC_meta, fit$se, fit$pval, fit$I2, fit$tau2)
  }, error = function(e){
    c(NA, NA, NA, NA, NA)
  })
  return(out)
}

# rodar para todos os genes
results <- t(sapply(1:nrow(metafor_filtered), meta_gene))
falhas <- sum(is.na(results[,1]))
cat(paste0(falhas," (",round((falhas/nrow(metafor_filtered))*100,2),"%) genes falharam.\n"))

# transformar em dataframe
results <- as.data.frame(results)

colnames(results) <- c(
  "logFC_meta",
  "SE_meta",
  "p_meta",
  "I2",
  "tau2"
)
results$Gene <- rownames(metafor_filtered)

# anotar em gene symbol
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = results$Gene,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)
results$Symbol <- gene_symbols

# remove genes NA, sem anotação
results <- results[!is.na(results$Symbol), ]



# ============== TRATAMENTO DE DADOS E ESTIMAÇÃO DE DEGS ==============
# Nesta seção foi utilizado os resultados da meta-análise para calcular o valor
# Z, Aplicar o método de False Discovery Rate sobre o valor p e extimado os genes
# diferencialmente expressos.
# Critérios: |Log2FC| > 1 & p-value(FDR) <= 0,05 


# Aplicar FDR
results$FDR <- p.adjust(results$p_meta, method = "BH")
results <- results[!is.na(results$FDR), ]

# filtrar DEGs finais
DEGs <- subset(results,
               abs(logFC_meta) > 1 & FDR < 0.05)
cat(paste0("Foram obtidos ",nrow(DEGs)," genes diferencialmente expressos!\n"))

# filtando pelo resultado de I²
DEGs_filtered <- subset(DEGs,DEGs$I2==0)
cat(paste0("Foram obtidos ",nrow(DEGs_filtered)," genes diferencialmente expressos com I² igual a 0!\n"))

# filtrando por up-regulated
DEGs_filtered <- subset(DEGs_filtered,DEGs_filtered$logFC_meta>=1)
cat(paste0("Foram obtidos ",nrow(DEGs_filtered)," genes diferencialmente expressos up-regulados e com I² igual a 0!\n"))

# Salvando dados

#todos os genes
write.csv(results, "Results/results_meta.csv", row.names = FALSE)
#genes diferencialmente expressos
write.csv(DEGs, "Results/DEGs_meta.csv", row.names = FALSE)


# ============== GRÁFICOS DE EXPRESSÃO ==============
# Aqui são gerados os gráficos de expressão como MA e volcano plot, além de histograma
# dos valores de I² identificados dos genes totais e DEGs

# Análise de heterogeneidade
# heterogeneidade global
summary(results$I2)

## ==== HISTOGRAMA I² GLOBAL ====
#pico 0: genes estáveis que não se expressam diferencialmente. pico perto de 100: genes de expressão variável 

png("Images/histograma_Genes_totais_Bexiga.png", width = 3000, height = 2000, res = 300)

hist(results$I2,
     breaks = 50,
     col = "#69b3a2",
     border = "white",
     main = "Distribuição de I² em genes totais de câncer de bexiga",
     xlab = "Heterogeneidade entre estudos (I²%)",
     ylab = "Frequência",
     cex.main = 1.2,
     cex.lab = 1,
     cex.axis = 0.9)

dev.off()

# heterogeneidade de DEGs
summary(DEGs$I2)

## ==== HISTOGRAMA I² DEGS ====

png("Images/histograma_DEGs_Bexiga.png", width = 3000, height = 2000, res = 300)

hist(DEGs$I2,
     breaks = 50,
     col = "#69b3a2",
     border = "white",
     main = "Distribuição de I² em genes diferencialmente expressos de câncer de bexiga",
     xlab = "Heterogeneidade entre estudos (I²%)",
     ylab = "Frequência",
     cex.main = 1.2,
     cex.lab = 1,
     cex.axis = 0.9)

dev.off()

## ==== VOLCANO PLOT ====

results$significance <- "Not significant"
results$significance[
  results$FDR < 0.05 & results$logFC_meta > 1
] <- "Upregulated"
results$significance[
  results$FDR < 0.05 & results$logFC_meta < -1
] <- "Downregulated"
results$log10FDR <- -log10(results$FDR)

png("Images/meta_volcano_plot_bexiga.png", width = 8, height = 8, units = "in", res = 300)

ggplot(results, aes(x = logFC_meta, y = log10FDR, color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Downregulated" = "#2C7BB6",
    "Upregulated" = "#D41159",
    "Not significant" = "lightgrey"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    x = "log2 Fold Change",
    y = "-log10(FDR)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

dev.off()



# ============== ENRIQUECIMENTO FUNCIONAL ==============

# lista de genes diferencialmente expressos
genes_degs <- unique(DEGs$Gene)

# lista de todos os genes como universo
genes_background <- unique(results$Gene)


## ==== GO ENRICHMENT ====

ego_bp <- enrichGO(
  gene          = genes_degs,
  universe      = genes_background,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

ego_cc <- enrichGO(
  gene = genes_degs,
  universe = genes_background,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

ego_mf <- enrichGO(
  gene = genes_degs,
  universe = genes_background,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

## ==== KEGG ====

ekegg <- enrichKEGG(
  gene         = genes_degs,
  universe     = genes_background,
  organism     = "hsa",
  pvalueCutoff = 0.05
)

# converter para símbolos
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

## ==== REACTOME ====

ereact <- enrichPathway(
  gene          = genes_degs,
  universe      = genes_background,
  organism      = "human",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

## alinhamento de genes robustos com vias funcionais enriquecidas

# colentando os genes robustos (heterogeneidade e up-regulados)
genes_robustos <- unique(DEGs_filtered$Symbol)

# função para verificar presença dos genes robustos nas vias
check_genes <- function(enrich_result, genes_robustos) {
  
  df <- as.data.frame(enrich_result)
  
  df$genes_presentes <- sapply(df$geneID, function(x) {
    genes <- unlist(strsplit(x, "/"))
    intersect(genes, genes_robustos)
  })
  
  df$n_genes_robustos <- sapply(df$genes_presentes, length)
  
  return(df)
}

bp_DEGs    <- check_genes(ego_bp, genes_robustos)
cc_DEGs    <- check_genes(ego_cc, genes_robustos)
mf_DEGs    <- check_genes(ego_mf, genes_robustos)
kegg_DEGs  <- check_genes(ekegg, genes_robustos)
react_DEGs <- check_genes(ereact, genes_robustos)

## ==== PLOTS ====

png("Images/GO_BP_dotplot.png", width = 3000, height = 2000, res = 300)
dotplot(ego_bp, showCategory = 20)
dev.off()

png("Images/GO_CC_dotplot.png", width = 3000, height = 2000, res = 300)
dotplot(ego_cc, showCategory = 20)
dev.off()

png("Images/GO_MF_dotplot.png", width = 3000, height = 2000, res = 300)
dotplot(ego_mf, showCategory = 20)
dev.off()

png("Images/KEGG_dotplot.png", width = 3000, height = 2000, res = 300)
dotplot(ekegg, showCategory = 20)
dev.off()

png("Images/REACTOME_dotplot.png", width = 3000, height = 2000, res = 300)
dotplot(ereact, showCategory = 20)
dev.off()

