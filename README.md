# Meta-análise transcriptômica integrativa para descoberta de biomarcadores diagnósticos em câncer de bexiga

## Visão Geral

Este projeto executa uma meta-análise integrada a nível transcriptômico de diferentes projetos obtidos da base de dados Gene Expression Omnibus (GEO) do NCBI.
Este pipeline inlcui:
- Meta-análise de expressão diferencial com integração de diferentes datasets.
- Análise de enriquecimento funcional de genes diferencialmente expressos.
- Análise de rede de interação proteína-proteína.
- Validação externa por análise diferencial através da base de dados Recount3 (TCGA + GTex).
- Avaliação de biomarcadores por ROC/AUC.
- Priorização de biomarcadores multi-criterial.

## Objetivo biológico

O objetivo deste projeto é identificar biomarcadores transcriptômicos gênicos robustos associados ao câncer de bexiga através da integração de múltiplos datasets independentes.

## Pipeline de análise

1. Coleta de datasets GEO:
   - GSE7476
   - GSE76211
   - GSE3167
   - GSE65635
   - GSE37817
   - GSE13507
   - GSE52519
2. Processamento, normalização e análise diferencial individual de cada dataset
3. Aplicação de meta-análise de efeitos aleatórios entre estudos
4. Análise de enriquecimento funcional
   - GO
   - KEGG
   - Reactome
5. Construção de rede de interação proteína-proteína (STRINGdb)
6. Avaliação de estado de significância de genes em datasets individuais
7. Validação externa TCGA + GTex (Recount3)
8. Avaliação de diagnóstico ROC/AUC a partir da análise de validação
9. Integração e ranqueamento de candidatos a biomarcadores de câncer de bexiga

---

## Ferramentas de Software utilizadas

- R (4.6.0)
- RStudio (2026.4.0.526)
- Bioconductor (3.23)

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

### Enriquecimento funcional
- clusterProfiler
- ReactomePA
- enrichplot

### Análise de rede PPI
- STRINGdb
- igraph

### Avaliação de desempenho diagnóstico
- pROC

### Manipulação e visualização de dados
- tidyverse (dplyr, tidyr, tibble, ggplot2)

## Reprodutibilidade

Informação de versões completas de pacotes e sessão está disponível através de: 

```r
sessionInfo()
```

R Environment também está disponível

---

## Execução e Reprodutibilidade

### Windows
#### Pré-requisitos
- Certifique-se de ter o [RTools](https://cran.r-project.org/bin/windows/Rtools/) instalado no sistema (necessário para compilação de pacotes de bioinformática no Windows).

Clonar repositório pelo Powershell em diretório de preferência:

```powershell
# suporte de caminhos longos
git config --global core.longpaths true
git clone https://github.com/RioGen-Tecnologia/Projeto_TCC_2026_Guilherme_Miranda.git
cd Projeto_TCC_2026_Guilherme_Miranda
```

Executar R:
```powershell
# substitua "R-4.6.0" pela sua versão atual do R
& "C:\Program Files\R\R-4.6.0\bin\x64\R.exe"
```

Intalar dependências:

```r
renv::restore()
```

Executar programa:

```r
source("scripts/bladder_cancer_meta_analysis.r")
```

---

### Linux

Clonar repositório em diretório de preferência:

```bash
git clone https://github.com/RioGen-Tecnologia/Projeto_TCC_2026_Guilherme_Miranda.git
cd Projeto_TCC_2026_Guilherme_Miranda
```

Iniciar R

```bash
R
```

Instalar dependências:

```r
renv::restore()
```

Executar programa:

```r
source("scripts/bladder_cancer_meta_analysis.r")
```

---

## Resultados principais

...

## Autor

Guilherme Moret Miranda - RioGen Tecnologia
