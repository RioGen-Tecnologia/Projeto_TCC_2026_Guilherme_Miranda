# ============== CARREGANDO PACOTES ==============

library(here) #pacote "aqui" auxilia com diretórios
library(GEOquery) #busca e download de datasets do Gene Omnibus
library(affy) #pacote de normalização affymetrix
library(oligo) #pacote de normalização oligo
library(limma) #análise estatística e normalização
library(AnnotationDbi) # pacote de execução de anotação
library(hgu133plus2.db) #pacote de base de anotação
library(hta20transcriptcluster.db) #pacote de base de anotação
library(hgu133a.db) #pacote de base de anotação
library(hgu133acdf) #pacote de base de anotação
library(hgu133plus2cdf) #pacote de base de anotação
library(illuminaHumanv2.db) #pacote de base de anotação
library(illuminaHumanv3.db) #pacote de base de anotação
library(illuminaHumanv4.db) #pacote de base de anotação
library(org.Hs.eg.db) #pacote de base de anotação
library(metafor) #pacote de meta-análise
library(ggplot2) #pacote de gráficos de expressão
library(dplyr) #gerenciamento de dataframes
library(tidyr) #gerenciamento de dataframes
library(clusterProfiler) #enriquecimento funcional (GO e KEGG)
library(ReactomePA) #enriquecimento funcional (Reactome)
library(enrichplot) #pacote de gráficos de enriquecimento
library(STRINGdb) #pacote de rede PPi
library(igraph) #complemento do grafico de rede PPI
library(recount3) #base de dados de TCGA + GTex padronizados
library(edgeR) #complemento aos dados do TCGA + GTex

# complementos a outros pacotes
library(R.utils) #ferramentas do R
library(Biobase)
library(BiocGenerics)
library(generics)
library(stats4)
library(IRanges)
library(S4Vectors)
library(R.oo)
library(R.methodsS3)
library(oligoClasses)
library(Biostrings)
library(DBI)
library(RSQLite)

# ==================================================
# Opções
# ==================================================

options(stringsAsFactors = FALSE)

set.seed(123)