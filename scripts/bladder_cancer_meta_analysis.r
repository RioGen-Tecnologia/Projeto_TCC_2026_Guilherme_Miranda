
# Análise integrada de projetos de microarray obtidas do Gene Expression Omnibus
# Guilherme Moret Miranda - Riogen
# 25/04/2026
# Esse script obtém, normaliza e análisa cada um dos projetos buscados previamente,
# combina os dados de cada projeto por meta-análise e realiza posteriores análises
# para determinar biomarcadores gênicos de câncer de bexiga.

# ============== CARREGANDO PACOTES ==============

message("===========================================================================")
message("CARREGANDO PACOTES")
message("===========================================================================")

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
  pattern = "samples_metadata_bladder_cancer_TCC-2026.csv",
  full.names = TRUE
)

# limpando
rm(d,dirs)
gc()

# ============== ANÁLISE DE PROJETOS ==============
# Esta sessão executa os scripts secundários que extrem os dados do GEO, normalizam,
# anotam e análisam estatísticamente por limma cada um dos projetos.

message("===========================================================================")
message("ANALISANDO DATASETS")
message("===========================================================================")

# Lista scripts que começam com "GSE"
scripts_projetos <- list.files(
  path = scripts_dir,
  pattern = "^GSE.*\\.[Rr]$",
  full.names = TRUE
)

# Loop de execução dos scripts secundários
arquivos_metafor <- list()
arquivos_deg <- list()

for (script in scripts_projetos) {
  message("Processando: ", basename(script))
  env <- new.env()
  source(script, local = env)
  id_atual <- gsub(".*(GSE[0-9]+).*", "\\1", basename(script))
  nome_objeto <- paste0("metafor_", id_atual)
  nome_deg <- paste0("deg_", id_atual)
  # ===== metafor =====
  if (exists(nome_objeto, envir = env)) {
    arquivos_metafor[[id_atual]] <- get(nome_objeto, envir = env)
  } else {
    warning(nome_objeto, " não encontrado.")
  }
  # ===== DEG =====
  if (exists(nome_deg, envir = env)) {
    arquivos_deg[[id_atual]] <- get(nome_deg, envir = env)
  } else {
    warning(nome_deg, " não encontrado.")
  }
  rm(env)
  gc()
}

#limpando
rm(scripts_projetos,script,nome_objeto,id_atual,nome_deg)
gc()

# ============== PREPARAÇÃO PARA META-ANÁLISE ==============
# Nesta seção os dados dos diferentes projetos são integrados.

message("===========================================================================")
message("PREPARANDO DADOS PARA META-ANÁLISE")
message("===========================================================================")

# extraindo os nomes de todos os genes de todos os datasets num único vetor
todos_genes <- unlist(lapply(arquivos_metafor, rownames))

# contando em quantos datasets cada gene aparece
contagem_genes <- table(todos_genes)

# filtrando apenas os genes que aparecem em 3 ou mais datasets
genes_filtrados <- names(contagem_genes[contagem_genes >= 3])
message("Número de genes retidos para análise: ", length(genes_filtrados))


# criando matrizes vazias (com NA) para alinhar os dados
num_genes <- length(genes_filtrados)
num_datasets <- length(arquivos_metafor)
nomes_datasets <- names(arquivos_metafor)

matriz_SMD   <- matrix(NA, nrow = num_genes, ncol = num_datasets, dimnames = list(genes_filtrados, nomes_datasets))
matriz_SE    <- matrix(NA, nrow = num_genes, ncol = num_datasets, dimnames = list(genes_filtrados, nomes_datasets))
matriz_logFC <- matrix(NA, nrow = num_genes, ncol = num_datasets, dimnames = list(genes_filtrados, nomes_datasets))

# preenchendo as matrizes iterando sobre os datasets
for (i in seq_along(arquivos_metafor)) {
  df <- arquivos_metafor[[i]]
  # descobrir quais dos genes filtrados existem neste dataset específico
  genes_em_comum <- intersect(genes_filtrados, rownames(df))
  # injetar os valores nas posições corretas da matriz
  matriz_SMD[genes_em_comum, i]   <- df[genes_em_comum, "SMD"]
  matriz_SE[genes_em_comum, i]    <- df[genes_em_comum, "SE"]
  matriz_logFC[genes_em_comum, i] <- df[genes_em_comum, "logFC"]
}


# limpando
rm(num_genes,num_datasets,nomes_datasets,df)
gc()

# ============== META-ANÁLISE ==============
# Salvamos o input da meta-análise e executamos a meta-análise devidamente utilizando
# o pacote metafor.
# Foi utilizado como critério genes que aparecem em no mínimo 3 projetos simultaneamente.

message("===========================================================================")
message("REALIZANDO META-ANÁLISE...")
message("===========================================================================")

# rodar para todos os genes
resultados_meta <- pblapply(1:length(genes_filtrados), function(i) {
  
  yi  <- matriz_SMD[i, ]
  sei <- matriz_SE[i, ]
  
  estudos_validos <- sum(!is.na(yi) & !is.na(sei))
  
  if (estudos_validos < 3) {
    return(c(SMD_pool = NA, SE_pool = NA, pval = NA, I2 = NA, tau2 = NA, num_estudos = estudos_validos))
  }
  
  # Usando REML (padrão mais moderno que conversamos)
  fit <- tryCatch({
    rma.uni(yi = yi, sei = sei, method = "REML")
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(fit)) {
    return(c(SMD_pool = NA, SE_pool = NA, pval = NA, I2 = NA, tau2 = NA, num_estudos = estudos_validos))
  }
  
  c(
    SMD_pool    = as.numeric(fit$beta),
    SE_pool     = fit$se,
    pval        = fit$pval,
    I2          = fit$I2,
    tau2        = fit$tau2,
    num_estudos = estudos_validos
  )
})

message("===========================================================================")
message("META-ANÁLISE FINALIZADA")
message("===========================================================================")

## transformando em uma tabela final
df_meta_final <- as.data.frame(do.call(rbind, resultados_meta))
df_meta_final$Gene <- genes_filtrados

# calculando a média do logFC ignorando os NAs
df_meta_final$logFC_mean <- rowMeans(matriz_logFC, na.rm = TRUE)

# Remover eventuais genes que falharam no modelo (NAs) antes de corrigir o P-valor
df_meta_final <- df_meta_final[!is.na(df_meta_final$pval), ]

# Correção FDR (Benjamini-Hochberg)
df_meta_final$FDR <- p.adjust(df_meta_final$pval, method = "BH")

# anotando os genes
simbolos_mapeados <- mapIds(
  org.Hs.eg.db,
  keys = as.character(df_meta_final$Gene),
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first" # Se houver mais de um símbolo, pega o primeiro
)

# criando a coluna symbol
df_meta_final$Symbol <- simbolos_mapeados[as.character(df_meta_final$Gene)]

# Reordenar colunas
df_meta_final <- df_meta_final[, c("Gene", "Symbol", "num_estudos", "logFC_mean", "SMD_pool", "SE_pool", "pval", "FDR", "I2", "tau2")]
df_meta_final <- df_meta_final[order(df_meta_final$pval), ]

print(head(df_meta_final))

# limpando
rm(matriz_SMD, matriz_SE, matriz_logFC, resultados_meta, todos_genes, contagem_genes, 
  genes_filtrados, i,estudos_validos,simbolos_mapeados)
gc()

# ============== TRATAMENTO DE DADOS E ESTIMAÇÃO DE DEGS ==============
# Nesta seção foi utilizado os resultados da meta-análise para calcular o valor
# Z, Aplicar o método de False Discovery Rate sobre o valor p e extimado os genes
# diferencialmente expressos.
# Critérios: |Log2FC| > 1 & p-value(FDR) <= 0,05 


# filtrar DEGs finais
DEGs <- subset(df_meta_final,
               abs(logFC_mean) > 1 & FDR < 0.05)
cat(paste0("Foram obtidos ",nrow(DEGs)," genes diferencialmente expressos!\n"))

# filtrando por up-regulated e I²<50%
DEGs_filtered <- subset(DEGs, I2 <= 50)

# Salvando dados
# confeindo rapidamente se a pasta de salvamento está pronta
out_dir <- file.path(results_dir, "DEGs_tables")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
rm(out_dir)

#todos os genes
write.csv(df_meta_final,
          file = file.path(results_dir,"meta_analysis","resultados_meta_analise.csv"),
          row.names = FALSE)

#genes diferencialmente expressos
write.csv(DEGs,
          file = file.path(results_dir,"meta_analysis","DEGs_meta_analise.csv"),
          row.names = FALSE)


# ============== GRÁFICOS DE EXPRESSÃO ==============
# Aqui são gerados os gráficos de expressão como MA e volcano plot, além de histograma
# dos valores de I² identificados dos genes totais e DEGs

#renomeando objeto resultados para facilitar o pipeline
results <- df_meta_final
rm(df_meta_final)

# Análise de heterogeneidade
# heterogeneidade global
summary(results$I2)

## ==== HISTOGRAMA I² GLOBAL ====
#pico 0: genes estáveis que não se expressam diferencialmente. pico perto de 100: genes de expressão variável

message("===========================================================================")
message("GERANDO GRÁFICOS DE EXPRESSÃO")
message("===========================================================================")

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
  results$FDR < 0.05 & results$logFC_mean > 1
] <- "Upregulated"
results$significance[
  results$FDR < 0.05 & results$logFC_mean < -1
] <- "Downregulated"
results$log10FDR <- -log10(results$FDR)


png(file.path(figures_dir, "meta_volcano_plot_bexiga.png"), width = 8, height = 8, units = "in", res = 300)

ggplot(results, aes(x = logFC_mean, y = log10FDR, color = significance)) +
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

message("===========================================================================")
message("REALIZANDO ENRIQUECIMENTO FUNCIONAL")
message("GENE ONTOLOGY - BIOLOGICAL PROCESS")
message("===========================================================================")

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

ego_bp <- clusterProfiler::simplify(
  ego_bp,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

## ==== KEGG ====

message("===========================================================================")
message("REALIZANDO ENRIQUECIMENTO FUNCIONAL")
message("KEGG")
message("===========================================================================")

ekegg <- clusterProfiler::enrichKEGG(
  gene         = genes_degs,
  universe     = genes_background,
  organism     = "hsa",
  pvalueCutoff = 0.05
)

# converter para símbolos
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

## ==== REACTOME ====

message("===========================================================================")
message("REALIZANDO ENRIQUECIMENTO FUNCIONAL")
message("REACTOME")
message("===========================================================================")

ereact <- enrichPathway(
  gene          = genes_degs,
  universe      = genes_background,
  organism      = "human",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

## alinhamento de genes robustos com vias funcionais enriquecidas

message("===========================================================================")
message("REALIZANDO ENRIQUECIMENTO FUNCIONAL")
message("ANALISANDO GENES E VIAS")
message("===========================================================================")

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
kegg_DEGs  <- check_genes(ekegg, genes_degs_symbol)
react_DEGs <- check_genes(ereact, genes_degs_symbol)


# unindo as vias enriquecidas
all_pathways <- bind_rows(
  bp_DEGs %>% mutate(Database = "GO_BP"),
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


## ==== plots ====

message("===========================================================================")
message("REALIZANDO ENRIQUECIMENTO FUNCIONAL")
message("GERANDO GRÁFICOS DE ENRIQUECIMENTO")
message("===========================================================================")

## dotplots
png(file.path(figures_dir, "GO_BP_dotplot.png"), width = 3000, height = 2000, res = 300)
dotplot(ego_bp, showCategory = 20)
dev.off()

png(file.path(figures_dir, "KEGG_dotplot.png"), width = 3000, height = 2000, res = 300)
dotplot(ekegg, showCategory = 20)
dev.off()

png(file.path(figures_dir, "REACTOME_dotplot.png"), width = 3000, height = 2000, res = 300)
dotplot(ereact, showCategory = 20)
dev.off()

## cnetplots
gene_fc <- DEGs$logFC_mean
names(gene_fc) <- DEGs$Symbol

# GO Biologiacal Process
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

# reactome
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


## ==== Exportando dados ====

message("===========================================================================")
message("REALIZANDO ENRIQUECIMENTO FUNCIONAL")
message("EXPORTANDO DADOS")
message("===========================================================================")

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

#limpando
rm(genes_degs,genes_background,ego_bp,ekegg,ereact,check_genes,bp_DEGs,kegg_DEGs,react_DEGs,all_pathways,
   all_pathways_long,gene_fc)
gc()


# ============== PPI NETWORK ==============
# análise de rede de interação proteína-proteína dos genes diferencialmente expressos
# com o objetivo de encontrar "hub genes", genes centrais na rede tumoral ou que
# parecem coordenar múltiplos processos.

message("===========================================================================")
message("ANÁLISE PROTEIN-PROTEIN INTERACTION")
message("INICIANDO STRINGdb")
message("===========================================================================")

# Incializando o STRING
string_db <- STRINGdb$new(version="12.0", 
                          species=9606, # _Homo sapiens_
                          score_threshold=400, 
                          input_directory="")

message("===========================================================================")
message("ANÁLISE PROTEIN-PROTEIN INTERACTION")
message("MAPEANDO GENES E INBTRERAÇÕES")
message("===========================================================================")

# Mapeando genes diferencialmente expressos pelos Entrez IDs
degs_mapped <- string_db$map(DEGs, "Gene", removeUnmappedRows = TRUE)

# Filtrando a lista por I² < 50% e garantindo expressão diferencial por FDR e LogFC
degs_filtered <- degs_mapped %>%
  filter(FDR < 0.05,
         I2 < 50,
         abs(logFC_mean) > 1.0)

# Obtendo as interações atrarvés dos IDs
hits <- degs_filtered$STRING_id
interactions <- string_db$get_interactions(hits)

# Cálculo de Hubs (Grau de Conectividade)
all_nodes <- c(interactions$from, interactions$to)
node_degree <- as.data.frame(table(all_nodes))
colnames(node_degree) <- c("STRING_id", "degree")

message("===========================================================================")
message("ANÁLISE PROTEIN-PROTEIN INTERACTION")
message("CALCULANDO HUB GENES")
message("===========================================================================")

# Mesclar com os dados originais
hubs_table <- merge(degs_filtered, node_degree, by="STRING_id") %>%
  arrange(desc(degree))

## checando se o diretório existe
out_dir <- file.path(results_dir, "ppi")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
rm(out_dir)

message("===========================================================================")
message("ANÁLISE PROTEIN-PROTEIN INTERACTION")
message("GERANDO GRÁFICO E EXPORTANDO DADOS")
message("===========================================================================")

## ==== plot ====
png(file.path(figures_dir, "ppi_DEGs_TCC_2026.png"),width = 5000,height = 4000,res = 600)

par(mar = c(1,1,3,1),cex = 1.4,lwd = 2)
string_db$plot_network(hubs_table$STRING_id)

dev.off()

## ==== Exportando dados ====
write.csv(interactions, file.path(results_dir, "ppi",  "string_interactions_edges.csv"), row.names = FALSE)
write.csv(hubs_table, file.path(results_dir, "ppi",  "string_nodes_metadata.csv"), row.names = FALSE)

# limpando
rm(degs_mapped,all_nodes,node_degree,string_db,degs_filtered,hits)
gc()


# ============== SIGNIFICÂNCIA EM DATASETS INDIVIDUAIS ==============
# é extraído e análisado em quantos datasets os genes foram considerados significativos
# e em quantos estavam presentes para evitar viés.
# Cálculo: (datasets signficativos / datasets totais) * (datasets presentes / datasets totais)

message("===========================================================================")
message("ANALISANDO A SIGNIFICANCIA E A PRESENÇA DE GENES ENTRE DATASETS... ")
message("===========================================================================")

# criar lista temporária com apenas gene + significance
lista_significance <- lapply(names(arquivos_deg), function(projeto) {
  df <- arquivos_deg[[projeto]][, c("ENTREZID", "significance")]
  colnames(df) <- c(
    "gene",
    paste0("significance_", projeto)
  )
  df
})

# merge progressivo por gene
significance_results <- Reduce(
  function(x, y) merge(x, y, by = "gene", all = TRUE),
  lista_significance
)

# colunas de significance
cols_significance <- grep(
  "^significance_",
  colnames(significance_results),
  value = TRUE
)

# número total de datasets (8)
n_total <- 7

# número de datasets significativos
n_significant <- apply(
  significance_results[, cols_significance],
  1,
  function(x) sum(x == "significant", na.rm = TRUE)
)

# número de datasets presentes (não-NA)
n_present <- apply(
  significance_results[, cols_significance],
  1,
  function(x) sum(!is.na(x))
)

# score final
significance_results$score <- (
  n_significant / n_total
) * (
  n_present / n_total
)

# visualizar
head(significance_results)

## checando se o diretório existe
out_dir <- file.path(results_dir, "significance")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
rm(out_dir)

#exportando tabela
write.csv(significance_results, file.path(results_dir, "significance",  "significance_results.csv"), row.names = FALSE)

# limpando
rm(n_present,n_significant,n_total,cols_significance,lista_significance)
gc()


# ============== Validação no TCGA ==============
# aqui é feito uma análise diferencial de dados padronizados de RNA-seq da base de
# dados Recount3 que possui amostras BLCA tumorais do TCGA e amostras não-tumorais
# de bexiga do GTex padronizadas por monorail. Os dados pré-processados foram
# extraídos e análisados por limma (voom) para servirem como validação dos resultados.

message("===========================================================================")
message("INICIANDO ANÁLISE DE VALIDAÇÃO...")
message("BASE DE DADOS RECOUNT3 (TCGA / GTEX)")
message("===========================================================================")

# rodando o script em novo ambiente e extraíndo os resultados
recount <- new.env()
source(file.path(scripts_dir,"TCGA + GTex validation.r"),local = recount)
validation_results <- get("results", envir = recount)
v <- get("v", envir = recount)
group_roc <- get("group_roc", envir = recount)
group <- get("group", envir = recount)
validation_exprs <- v$E
# remover Normal_Adj
keep_samples <- group != "Normal_Adj"
validation_exprs <- validation_exprs[, keep_samples]
group_roc <- group_roc[keep_samples]
rm(v,group,keep_samples,recount)
gc()

# filtando os resultados de validação com os genes da análise principal
validation_results_filtered <- validation_results[
  rownames(validation_results) %in% as.character(DEGs_filtered$Gene),
]

#salvando os dados
# confere rapidamente se a pasta de salvamento está pronta
out_dir <- file.path(results_dir, "TCGA + GTex validation")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
rm(out_dir)

write.csv(validation_results,
          file.path(results_dir, 'TCGA + GTex validation', "validation_analysis_results.csv"),row.names = TRUE)
write.csv(validation_results_filtered,
          file.path(results_dir,"TCGA + GTex validation","validation_analysis_comparison.csv"),row.names = TRUE)

# ============== ROC/AUC ==============
# calculado as curvas ROC para cada gene que refletem a sensibilidade e especificidade
# dos genes para identificar condição de tumor e não tumor e o AUC que é um valor
# resume e escala os dados ROC para serem usados como critério de qualidade
# do gene como candidato a biomarcador.

message("===========================================================================")
message("ANALISANDO CURVAS ROC E AUC...")
message("===========================================================================")

# Criar labels binários (tumor/não-tumor)
labels <- ifelse(group_roc == "tumor", 1, 0)

# estimando curvas AUC dos genes
roc_results <- lapply(rownames(validation_exprs), function(gene_id) {
  values <- as.numeric(validation_exprs[gene_id, ])
  r <- pROC::roc(labels, values, quiet = TRUE)
  data.frame(
    gene = gene_id,
    auc = as.numeric(pROC::auc(r))
  )
})

#resultado final
roc_df <- do.call(rbind, roc_results)

# filtrando para DEGs identificados
roc_filtered <- roc_df %>%
  filter(as.character(gene) %in% as.character(DEGs_filtered$Gene))

## checando se o diretório existe
out_dir <- file.path(results_dir, "auc")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
rm(out_dir)

#salvando resultados
write.csv(roc_filtered,
          file.path(results_dir, 'auc',"auc.csv"))


## criando um plot de exemplo
gene <- roc_filtered[1,1]
values <- as.numeric(validation_exprs[gene,])
roc_curve <- pROC::roc(labels, values)
auc_value <- auc(roc_curve)
png(file.path("figures",paste("ROC_curve_gene",gene,".png")),height = 1800, width = 1800,res = 300)
ggroc(roc_curve) +
  geom_abline(
    intercept = 1,
    slope = 1,
    linetype = "dashed",
    color = "gray"
  ) +
  ggtitle(paste("ROC curve - gene", gene)) +
  annotate(
    "text",
    x = 0.65,
    y = 0.2,
    label = paste("AUC =", round(auc_value, 3))
  ) +
  theme_minimal()
dev.off()

# limpando
rm(roc_results,roc_df,gene,values,roc_curve)
gc()

# ============== COMPILAÇÃO DE RESULTADOS ==============
# Unindo todos os resultados para fazer a pontuação.

message("===========================================================================")
message("COMPILANDO RESULTADOS")
message("===========================================================================")

# criando dataframe dos resultados compilados
compiled_results <- tibble("Entrez"=DEGs_filtered$Gene,
                   "Gene_Symbol"=DEGs_filtered$Symbol,
                   "SMD"=DEGs_filtered$SMD_pool,
                   "adjusted_p.value"=DEGs_filtered$FDR,
                   "I2"=DEGs_filtered$I2)

# adicionando rede PPI
compiled_results <- compiled_results %>%
  left_join(
    hubs_table %>%
      select(
        Gene,
        degree
      ) %>%
      rename(
        Entrez = Gene,
        PPI_degree = degree
      ),
    by = "Entrez"
  )

# adicionando resultados de signficancia
compiled_results <- compiled_results %>%
  left_join(
    significance_results %>%
      select(gene,score) %>%
      rename(
        Entrez = gene,
        significance_score = score
      ),
    by = "Entrez"
  )

# dados de validação
compiled_results <- compiled_results %>%
  left_join(
    validation_results_filtered %>%
      rownames_to_column("Entrez") %>%
      select(Entrez, SMD_val, FDR_val),
    by = c("Entrez" = "Entrez")
  )

# dados AUC
compiled_results <- compiled_results %>%
  left_join(
    roc_filtered %>%
      select(gene, auc),
    by = c("Entrez" = "gene")
  )

#convertendo a dataframe formal
compiled_results <- as.data.frame(compiled_results)

#removendo os NA de PPI
compiled_results$PPI_degree[is.na(compiled_results$PPI_degree)] <- 0

cat(paste0(sum(is.na(compiled_results$SMD_val)), " genes não foram identificados na validação.\n"))

#removendo genes ausentes na validação
compiled_results <- compiled_results[!is.na(compiled_results$SMD_val), ]

# ============== PONTUAÇÃO ==============
# fase final da análise. Baseado nas análises feitas como critérios cada gene será
# pontuado num ranking de melhores candidatos a biomarcadores.

message("===========================================================================")
message("IDENTIFCANDO BIOMARCADORES...")
message("===========================================================================")

# Função para normalizar variáveis de cada critério para escala 0–1
safe_rescale <- function(x, to = c(0,1)) {
  # se todos forem NA
  if(all(is.na(x))) {
    return(rep(NA_real_, length(x)))
  }
  # intervalo dos dados
  r <- range(x, na.rm = TRUE)
  # evita divisão por zero
  # caso todos valores sejam iguais
  if(isTRUE(all.equal(r[1], r[2]))) {
    return(rep(mean(to), length(x)))
  }
  # normalização
  rescale(x, to = to, from = r)
}

# função para reduzir impacto de outliers extremos.
clip_quant <- function(x,
                       probs = c(0.05, 0.95)) {
  # calcula percentis
  q <- quantile(
    x,
    probs = probs,
    na.rm = TRUE,
    names = FALSE
  )
  # limita os valores
  pmin(pmax(x, q[1]), q[2])
}


## FUNÇÃO DE PONTUAÇÃO
score_biomarkers <- function(df) {
  df %>%
    mutate(
      
      # CRITÉRIO 1: mangnitude biológica (SMD)
      # priorização de up-regulados com pmax(LogFC, 0)
      # Genes negativos viram 0
      SMD_pos = pmax(SMD, 0),
      SMD_val_pos = pmax(SMD_val, 0),
      
      # CRITÉRIO 2: valor p ajustado por FDR
      # valores menores refletem melhor evidência estatística
      # usou-se -log10(FDR)
      # 0.01   -> 2
      # 0.001  -> 3
      # 1e-20  -> 20
      meta_fdr_score = -log10(pmax(adjusted_p.value, 1e-300)),
      val_fdr_score = -log10(pmax(FDR_val, 1e-300)),
      
      # NORMALIZAÇÃO DAS VARIÁVEIS
      # aplicação da função em cada critério
      # Reduz outliers e converte para escala 0-1

      # magnitude biológica GEO
      meta_effect_n = safe_rescale( clip_quant(SMD_pos)),
      
      # magnitude biológica TCGA
      val_effect_n = safe_rescale(clip_quant(SMD_val_pos)),
      
      # significância estatística GEO
      meta_fdr_n = safe_rescale(clip_quant(meta_fdr_score)),
      
      # significância estatística TCGA
      val_fdr_n = safe_rescale(clip_quant(val_fdr_score)),
      
      # AUC diagnóstica
      auc_n = safe_rescale(clip_quant(auc)),
      
      # consistência entre datasets
      significance_n = safe_rescale(clip_quant(significance_score)),
      
      # conectividade em PPI
      ppi_n = safe_rescale(clip_quant(PPI_degree)),
      
      # heterogeneidade
      i2_n = safe_rescale(clip_quant(I2)),
      
      # CRITÉRIO 3: robustez entre estudos na meta-análise
      # Definido por I² que reflete a heterogeneidade entre datasets
      # valores menores são mais homogêneros
      # I² alto -> score baixo
      # I² baixo -> score alto
      robustness_n = 1 - i2_n,
      
      
      # CRITÉRIO 4: Concordância direcional
      # Analisa se os genes são up ou down regulados na análise e validação ao
      # mesmo tempo
      # +1 = mesma direção
      #  0 = direção oposta
      direction_bonus = ifelse(sign(SMD) == sign(SMD_val),1,0),
      
      # Bônus de prioridade a genes up regulados
      # interesse clínico
      # 1 = up-regulado nas duas análises
      # 0.5 = up-regulado em apenas uma
      # 0 = não up-regulado
      up_bonus = case_when(SMD > 0 & SMD_val > 0 ~ 1, SMD > 0 | SMD_val > 0 ~ 0.5,TRUE ~ 0),
      
      # SCORE FINAL
      # Os pesos podem ser ajustados
      biomarker_score =
        100 * (
          # magnitude GEO
          0.14 * meta_effect_n +
            # magnitude TCGA
            0.14 * val_effect_n +
            # significância GEO
            0.12 * meta_fdr_n +
            # significância TCGA
            0.12 * val_fdr_n +
            # consistência entre datasets
            0.12 * significance_n +
            # AUC diagnóstica
            0.18 * auc_n +
            # baixa heterogeneidade
            0.10 * robustness_n +
            # contexto biológico PPI
            0.03 * ppi_n +
            # concordância GEO ↔ TCGA
            0.03 * direction_bonus +
            # bônus translacional
            0.02 * up_bonus
          
        )) %>%
    
  # ordenba os resultados pelo score
  arrange(desc(biomarker_score))}

# execução da pontuação
ranked_results <- score_biomarkers(compiled_results)

# VISUALIZAÇÃO DOS TOP GENES
head(
  ranked_results[, c(
    "Gene_Symbol",
    "biomarker_score",
    "SMD",
    "adjusted_p.value",
    "I2",
    "significance_score",
    "SMD_val",
    "FDR_val",
    "auc",
    "PPI_degree"
  )]
)

## checando se o diretório existe
out_dir <- file.path(results_dir, "biomarker_results")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
rm(out_dir)

ranked_results_clean <- ranked_results[,c((1:10),26)]

message("=======================================")
message("EXPORTANDO RESULTADOS DE BIOMARCADOERES")
message("=======================================")

write.csv(ranked_results_clean,file.path(results_dir,"biomarker_results","Bladder_cancer_biomarker_rank.csv"))

# ===== SESSION INFO =====
writeLines(
  capture.output(sessionInfo()),
  file.path(results_dir, "sessionInfo.txt")
)

gc()


# ============== Gráficos de resultados ==============
# Gráficos dos resultados obtidos na análise

message("===========================================================================")
message("GERANDO GRÁFICOS DOS RESULTADOS!")
message("===========================================================================")

## heatmap de resultados
# representação gráfica de pontuação dos genes

# compilando a pontuação
ranked_results_points <- ranked_results[,c(2,(15:21),(23:26))]

# organizando os dados
# renomeando colunas
ranked_results_points <- ranked_results_points %>%
  rename(
    SMD = meta_effect_n,
    SMD_val = val_effect_n,
    FDR = meta_fdr_n,
    FDR_val = val_fdr_n,
    AUC = auc_n,
    consistencia = significance_n,
    ppi = ppi_n,
    I2 = robustness_n,
    concordancia_direcional = direction_bonus,
    bonus_up_regulado = up_bonus
  )

# reordenando colunas
ranked_results_points <- ranked_results_points %>%
  select(Gene_Symbol,SMD,FDR,I2,consistencia,ppi,SMD_val,
         FDR_val,concordancia_direcional,AUC,bonus_up_regulado,biomarker_score)

# dividindo os pontos de biomarcador por 100 para ficarem na escala 0-1
ranked_results_points$biomarker_score <- ranked_results_points$biomarker_score / 100

#criando heatmap de visualização de dados
heatmap_results <- ranked_results_points %>%
  rename(
    "SMD validação" = SMD_val,
    "FDR validação" = FDR_val,
    "I²" = I2,
    "concordãncia" = concordancia_direcional,
    "up regulado" = bonus_up_regulado,
    "pontuação final" = biomarker_score
  ) %>%
  column_to_rownames("Gene_Symbol") %>%
  as.matrix()
heatmap_results <- t(heatmap_results)

png(
  file.path("figures","biomarker_score_heatmap.png"),
  width = 4200,
  height = 1400,
  res = 300
)
Heatmap(
  heatmap_results,
  name = "Pontuação",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 9),
  col = colorRamp2(
    c(0, 0.5, 1),
    c("#D73027", "#FFFFBF", "#1A9850")
  )
)
dev.off()

## gráfico de radar
# mostrando os top genes num gráfico circular

radar_top_5 <- ranked_results_points[(1:5),c(1,2,3,6,7,8,10)] %>%
  rename(
    "SMD validação" = SMD_val,
    "FDR validação" = FDR_val,
    "PPI" = ppi
  )
radar_top_10 <- ranked_results_points[(6:10),c(1,2,3,6,7,8,10)] %>%
  rename(
    "SMD validação" = SMD_val,
    "FDR validação" = FDR_val,
    "PPI" = ppi
  )

png(
  file.path("figures","top_5_radar_chart.png"),
  width = 2500,height = 2000,res=300)
ggradar(radar_top_5,
        legend.position = "bottom")
dev.off()

png(
  file.path("figures","top_10_radar_chart.png"),
  width = 2500,height = 2000,res=300)
ggradar(radar_top_10,
        legend.position = "bottom")
dev.off()

## Heatmap de expressão
# foi feito um gráfico de heatmap para os 89 DEGs e top 20 genes candidatos a
# biomarcador por pontuação.

## organizando amostras
# ordenar agrupando non-tumor e tumor
sample_order <- order(group_roc)

# aplicar ordenação
exprs_ordered <- validation_exprs[, sample_order]
group_ordered <- group_roc[sample_order]


## HEATMAP DOS GENES FILTRADOS HOMOGÊNEOS

# genes em ENTREZ
genes_filtered <- as.character(DEGs_filtered$Gene)

# filtrar matriz
heatmap_filtered <- exprs_ordered[rownames(exprs_ordered) %in% genes_filtered,]

# ordenar linhas igual DEGs_filtered
heatmap_filtered <- heatmap_filtered[match(genes_filtered, rownames(heatmap_filtered)),]

# substituir ENTREZ por símbolo gênico
rownames(heatmap_filtered) <- DEGs_filtered$Symbol

# Z-score por gene
heatmap_filtered_scaled <- t(scale(t(heatmap_filtered)))

# remover possíveis NAs
heatmap_filtered_scaled <- heatmap_filtered_scaled[complete.cases(heatmap_filtered_scaled),]

# anotação de grupos
ha_filtered <- HeatmapAnnotation(
  Group = group_ordered,
  col = list(Group = c("non_tumor" = "#3A7D44", "tumor" = "#7B2CBF")))

# salvar figura
png(file.path("figures","heatmap_filtered_genes.png"),
    width = 3200,height = 4200,res = 400)

Heatmap(
  heatmap_filtered_scaled,
  name = "Z-score",
  top_annotation = ha_filtered,
  cluster_rows = TRUE,
  # NÃO clusterizar amostras
  cluster_columns = FALSE,
  # separar grupos visualmente
  column_split = group_ordered,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8),
  column_title = "TCGA Tumor vs GTEx Healthy",
  heatmap_legend_param = list(
    title = "Expression"
  ),
  col = colorRamp2(
    c(-2, 0, 2),
    c("#2C7BB6", "#F0F0F0", "#D41159")
  )
)

dev.off()

## Heatmap dos top 20 biomarcadores

# top 20 genes
top20 <- ranked_results_clean$Entrez[1:20] %>%as.character()

# filtrar matriz
heatmap_top20 <- exprs_ordered[rownames(exprs_ordered) %in% top20,]

# ordenar linhas
heatmap_top20 <- heatmap_top20[match(top20, rownames(heatmap_top20)),]

# símbolos gênicos
rownames(heatmap_top20) <- ranked_results_clean$Gene_Symbol[1:20]

# Z-score
heatmap_top20_scaled <- t(
  scale(t(heatmap_top20)))

# remover NAs
heatmap_top20_scaled <- heatmap_top20_scaled[complete.cases(heatmap_top20_scaled),]

# anotação
ha_top20 <- HeatmapAnnotation(
  Group = group_ordered,
  col = list(Group = c("non_tumor" = "#3A7D44", "tumor" = "#7B2CBF")))

# salvar figura
png(file.path("figures","heatmap_top20_genes.png"),
    width = 2600,height = 2200,res = 400)

Heatmap(
  heatmap_top20_scaled,
  name = "Z-score",
  top_annotation = ha_top20,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_split = group_ordered,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 10),
  column_title = "Top 20 Biomarkers",
  heatmap_legend_param = list(
    title = "Expression"
  ),
  col = colorRamp2(
    c(-2, 0, 2),
    c("#2C7BB6", "#F0F0F0", "#D41159")))

dev.off()

## limpando
rm(
  sample_order,
  exprs_ordered,
  genes_filtered,
  heatmap_filtered,
  top20,
  heatmap_top20
)
gc()


message("===========================================================================")
message("ANÁLISE FINALIZADA!")
message("===========================================================================")
