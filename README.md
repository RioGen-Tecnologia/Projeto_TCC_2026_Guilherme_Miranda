# Meta-análise transcriptômica integrativa para descoberta de biomarcadores de diagnóstico no câncer de bexiga

## Visão Geral

Este projeto executa uma meta-análise integrada a nível transcriptômico de diferentes projetos obtidos da base de dados Gene Expression Omnibus (GEO) do NCBI.
Este pipeline inclue:
- Meta-análise de expressão diferencial com integração de diferentes datasets.
- Análise de enriquecimento funcional de genes diferencilamente expressos.
- Análise de rede de interação proteína-proteína.
- Validação externa por análise diferencial através da base de dados Recount3 (TCGA + GTex).
- Avaliação de biomarcadores por ROC/AUC.
- Prirização de biomarcadores multi-criterial.

## pipeline de análise

1. Coleta de datasets GEO:
   - GSE7476
   - GSE76211
   - GSE3167
   - GSE65635
   - GSE37817
   - GSE13507
   - GSE52519
2. Processamento, normalização e análise estatística de expressão diferencial individualmente a cada datasets
3. Aplicação de meta-análise de efeitos aleatórios entre estudos
4. Análise de enriquecimento funcional
   - GO
   - KEGG
   - Reactome
5. Contrução de rede de interação proteína-proteína (STRINGdb)
6. Avaliação de estado de significância de genes em datasets individuais
7. Validação externa TCGA + GTex (Recount3)
8. Avaliação de diagnóstico ROC/AUC apartir da análise de validação
9. Integração e ranqueamento de candidatos a biomarcadores de câncer de bexiga

---

# Ferramentas e pacotes do R utilizados

## Software utilizado

- RStudio
- Bioconductor

## Pacotes R principais

### Aquisição de dados
- GEOquery
- recount3
  
### Processamento e análise diferencial
- affy
- oligo
- limma
- edgeR

### Anotação
- AnnotationDbi
- org.Hs.eg.db
- hgu133plus2.db
- hta20transcriptcluster.db
- hgu133a.db
- illuminaHumanv2.db
- illuminaHumanv4.db

### Meta-análise
- metafor

### enriquecimento funcional
- clusterProfiler
- ReactomePA
- enrichplot

### Análise de rede PPI
- STRINGdb
- igraph

### avaliação de desempenho diagnóstico
- pROC

### Manipulação e visualização de dados
- tidyverse (dplyr, tidyr, tibble, ggplot2)

## Reproducibility

Complete package versions and session information are available through:

```r
sessionInfo()
