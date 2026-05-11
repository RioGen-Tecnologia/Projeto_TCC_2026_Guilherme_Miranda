
# Análise integrada de projetos de microarray obtidas do Gene Expression Omnibus
# Guilherme Moret Miranda - Riogen
# 25/04/2026
# Esse script obtém, normaliza e análisa cada um dos projetos buscados previamente,
# combina os dados de cada projeto por meta-análise e realiza posteriores análises
# para determinar biomarcadores gênicos de câncer de bexiga.

# ============== CARREGANDO PACOTES ==============
library(here)
source(here("scripts", "setup.r"))

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
                           "input_metafor_TCC_2026.csv"),
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
    
    # usar exatamente os mesmos estudos do modelo
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
# confeindo rapidamente se a pasta de salvamento está pronta
out_dir <- file.path(results_dir, "DEGs_tables")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
rm(out_dir)

#todos os genes
write.csv(
  results,
  file = file.path(results_dir, "DEGs_tables", "global_meta.csv"),
  row.names = FALSE
)
#genes diferencialmente expressos
write.csv(
  DEGs,
  file = file.path(results_dir, "DEGs_tables", "DEGs_meta.csv"),
  row.names = FALSE
)


# ============== GRÁFICOS DE EXPRESSÃO ==============
# Aqui são gerados os gráficos de expressão como MA e volcano plot, além de histograma
# dos valores de I² identificados dos genes totais e DEGs

# Análise de heterogeneidade
# heterogeneidade global
summary(results$I2)

## ==== HISTOGRAMA I² GLOBAL ====
#pico 0: genes estáveis que não se expressam diferencialmente. pico perto de 100: genes de expressão variável 

png(file.path(figures_dir, "histograma_Genes_totais_Bexiga.png"), width = 3000, height = 2000, res = 300)

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

png(file.path(figures_dir, "histograma_DEGs_Bexiga.png"), width = 3000, height = 2000, res = 300)

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


png(file.path(figures_dir, "meta_volcano_plot_bexiga.png"), width = 8, height = 8, units = "in", res = 300)

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
# realiza enriquecimento funcional (Over-representation analysis) com genes
# diferencialmente expressos nas bases de dados KEGG, REACTOME e GO (biological
# process, Cell Component & Molecular function) para identificar candidatos a
# biomarcadores é identificado os genes presentes nas vias mais relevantes

# lista de genes diferencialmente expressos
genes_degs <- unique(DEGs$Gene)

# lista de todos os genes como universo
genes_background <- unique(results$Gene)


## ==== GO ENRICHMENT ====

# biological process
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

ego_bp <- simplify(
  ego_bp,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

# celular component
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

ego_cc <- simplify(
  ego_cc,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

# molecular function
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

ego_mf <- simplify(
  ego_mf,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
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
genes_degs_symbol <- unique(DEGs$Symbol)

# função para verificar presença dos genes robustos nas vias
check_genes <- function(enrich_result, genes_degs_symbol) {
  df <- as.data.frame(enrich_result)
  df$genes_presentes <- sapply(df$geneID, function(x) {
    genes <- unlist(strsplit(x, "/"))
    presentes <- intersect(genes, genes_degs_symbol)
    if(length(presentes) == 0) {
      return(NA)
    }
    paste(presentes, collapse = ", ")
  })
  df$n_genes_robustos <- sapply(
    strsplit(ifelse(is.na(df$genes_presentes),
                    "",
                    df$genes_presentes),
             ", "),
    function(x) sum(x != "")
  )
  return(df)
}

bp_DEGs    <- check_genes(ego_bp, genes_degs_symbol)
cc_DEGs    <- check_genes(ego_cc, genes_degs_symbol)
mf_DEGs    <- check_genes(ego_mf, genes_degs_symbol)
kegg_DEGs  <- check_genes(ekegg, genes_degs_symbol)
react_DEGs <- check_genes(ereact, genes_degs_symbol)


# unindo as vias enriquecidas
all_pathways <- bind_rows(
  bp_DEGs %>% mutate(Database = "GO_BP"),
  cc_DEGs %>% mutate(Database = "GO_CC"),
  mf_DEGs %>% mutate(Database = "GO_MF"),
  kegg_DEGs %>% mutate(Database = "KEGG"),
  react_DEGs %>% mutate(Database = "REACTOME")
)

# filtrando apenas vias contendo genes robustos
all_pathways_filtered <- all_pathways %>%
  filter(!is.na(genes_presentes)) %>%
  filter(n_genes_robustos > 0)
# ordenando vias mais significativas
all_pathways_filtered <- all_pathways_filtered %>%
  arrange(p.adjust)
# top vias por base de dados
top_pathways <- all_pathways_filtered %>%
  group_by(Database) %>%
  slice_min(order_by = p.adjust, n = 10)

# explodindo genes das vias
all_pathways_long <- all_pathways_filtered %>%
  separate_rows(genes_presentes, sep = ", ")

# score funcional
enrichment_scores <- all_pathways_long %>%
  group_by(genes_presentes) %>%
  summarise(
    n_vias = n(),
    databases = n_distinct(Database),
    mean_enrichment = mean(FoldEnrichment, na.rm = TRUE),
    best_padj = min(p.adjust, na.rm = TRUE)
  ) %>%
  arrange(desc(n_vias), best_padj)

## ==== EXPORTANDO DADOS ====

# confeindo rapidamente se a pasta de salvamento está pronta
out_dir <- file.path(results_dir, "enrichment")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
rm(out_dir)

write.csv(
  all_pathways_filtered,
  file.path(results_dir, "enrichment", "all_enriched_pathways.csv"),
  row.names = FALSE
)

write.csv(
  enrichment_scores,
  file.path(results_dir, "enrichment",  "gene_enrichment_scores.csv"),
  row.names = FALSE
)

## ==== PLOTS ====

# dotplots
png(file.path(figures_dir, "GO_BP_dotplot.png"), width = 3000, height = 2000, res = 300)
dotplot(ego_bp, showCategory = 20)
dev.off()

png(file.path(figures_dir, "KEGG_dotplot.png"), width = 3000, height = 2000, res = 300)
dotplot(ekegg, showCategory = 20)
dev.off()

png(file.path(figures_dir, "REACTOME_dotplot.png"), width = 3000, height = 2000, res = 300)
dotplot(ereact, showCategory = 20)
dev.off()

# cnetplots
gene_fc <- DEGs$logFC_meta
names(gene_fc) <- DEGs$Symbol

png(file.path(figures_dir, "GO_BP_cnetplot.png"),width = 3500,height = 3000,res = 450)
enrichplot::cnetplot(
  ego_bp,
  showCategory = 5,
  foldChange = gene_fc,
  node_label = "all",
  layout = "kk"
) +
  scale_color_gradientn(
    colours = c("#2C7BB6", "#D41159"),
    limits = c(-3, 3),
    oob = scales::squish,
    name = "logFC"
  ) +
  labs(size = "Gene count") +
  
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )
dev.off()

png(file.path(figures_dir, "REACTOME_cnetplot.png"),width = 3500,height = 3000,res = 450)
enrichplot::cnetplot(
  ereact,
  showCategory = 5,
  foldChange = gene_fc,
  node_label = "all",
  layout = "kk"
) +
  scale_color_gradientn(
    colours = c("#2C7BB6", "#D41159"),
    limits = c(-3, 3),
    oob = scales::squish,
    name = "logFC"
  ) +
  labs(size = "Gene count") +
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )
dev.off()



# ============== PPI NETWORK ==============
# análise de rede de interação proteína-proteína dos genes diferencialmente expressos
# com o objetivo de encontrar "hub genes", genes centrais na rede tumoral ou que
# parecem coordenar múltiplos processos.


# Incializando o STRING
string_db <- STRINGdb$new(version="12.0", 
                          species=9606, # _Homo sapiens_
                          score_threshold=400, 
                          input_directory="")

# Mapeando genes diferencialmente expressos pelos Entrez IDs
degs_mapped <- string_db$map(DEGs, "Gene", removeUnmappedRows = TRUE)

# Filtrando a lista por I² < 50% e garantindo expressão diferencial por FDR e LogFC
degs_filtered <- degs_mapped %>%
  filter(FDR < 0.05,
         I2 < 50,
         abs(logFC_meta) > 1.0)

# Obtendo as interações atrarvés dos IDs
hits <- degs_filtered$STRING_id
interactions <- string_db$get_interactions(hits)

# Cálculo de Hubs (Grau de Conectividade)
all_nodes <- c(interactions$from, interactions$to)
node_degree <- as.data.frame(table(all_nodes))
colnames(node_degree) <- c("STRING_id", "degree")

# Mesclar com os dados originais
hubs_table <- merge(degs_filtered, node_degree, by="STRING_id") %>%
  arrange(desc(degree))

# limpando
rm(degs_mapped,all_nodes,node_degree)
gc()

## checando se o diretório existe
out_dir <- file.path(results_dir, "ppi")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
rm(out_dir)


## ==== PLOT ====
png(file.path(figures_dir, "ppi_DEGs_TCC_2026.png"),width = 5000,height = 4000,res = 600)

par(mar = c(1,1,3,1),cex = 1.4,lwd = 2)
string_db$plot_network(hubs_table$STRING_id)

dev.off()

## ==== Resultados ====
write.csv(interactions, file.path(results_dir, "ppi",  "string_interactions_edges.csv"), row.names = FALSE)
write.csv(degs_filtered, file.path(results_dir, "ppi",  "string_nodes_metadata.csv"), row.names = FALSE)