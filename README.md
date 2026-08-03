# Meta-análise transcriptômica integrativa para descoberta de biomarcadores diagnósticos em câncer de bexiga

---

## Visão Geral

Este projeto executa uma meta-análise transcriptômica integrativa de diferentes datasets obtidos da base Gene Expression Omnibus (GEO) do NCBI.
Este pipeline inclui:
- Meta-análise de expressão diferencial com integração de diferentes datasets.
- Análise de enriquecimento funcional de genes diferencialmente expressos (ORA & GSEA).
- Análise de rede de interação proteína-proteína (PPI).
- Validação externa por análise diferencial através da base de dados Recount3 (TCGA + GTex).
- Avaliação de biomarcadores por ROC/AUC.
- Priorização de biomarcadores multi-criterial.

## Objetivo biológico

O objetivo deste projeto é identificar biomarcadores transcriptômicos gênicos robustos associados ao câncer de bexiga como marcadores diagnósticos através da integração de múltiplos datasets independentes.

---

## Pipeline de análise

1. Coleta de datasets GEO.
2. Processamento, normalização e análise diferencial individual de cada dataset

| Dataset  | Plataforma | Normalização | Anotação |
|----------|------------|--------------|----------|
| GSE7476  | GPL570     | RMA (affy)   | hgu133plus2.db |
| GSE76211 | GPL17586   | RMA (oligo)  | hta20transcriptcluster.db |
| GSE3167  | GPL96      | RMA (affy)   | hgu133a.db |
| GSE65635 | GPL14951   | neqc (limma) | illuminaHumanv4.db |
| GSE37817 | GPL6102    | já normalizado | illuminaHumanv2.db |
| GSE13507 | GPL6102    | normalizeBetweenArrays (limma) | anotação interna |
| GSE52519 | GPL6884    | neqc (limma) | anotação interna |

3. Integração e aplicação de meta-análise de efeitos aleatórios entre estudos
   * A meta-análise em modelo de efeitos aleatórios foi escolhida como estratégia principal de integração devido à heterogeneidade entre plataformas transcriptômicas e protocolos experimentais, permitindo combinar evidências estatísticas preservando efeitos específicos de cada estudo.
   * Estimativa combinada de:
     * g de Hedges
     * erro padrão
     * significância estatística (false discovery rate)
     * heterogeneidade entre estudos (I² e tau²) (Restricted Maximum Likelihood)

4. Análise de enriquecimento funcional
   - GO
   - KEGG
   - Reactome

5. Construção de rede de interação proteína-proteína (STRINGdb)
   * Objetivos
     * identificar hub genes
     * avaliar conectividade funcional
     * detectar genes potencialmente centrais em processos tumorais

6. Avaliação de consistência de genes entre datasets individuais
   * Objetivos
     * identificar se o gene é significativo dentro de cada dataset
     * Adicionar o critério da presença dos genes entre os datasets
     * avaliar robuistez do gene

7. Validação externa TCGA + GTex (Recount3)
   * Objetivos
     * Validar os genes identificados como candidatos a biomarcadores
     * Aumentar a robustez da seleção de genes

8. Avaliação de diagnóstico ROC/AUC a partir da análise de validação
   * Objetivos
     * Avaliar a capacidade discriminatória dos genes entre amostras tumorais e não-tumorais
     * Priorizar genes com maior potencial diagnóstico

9. Integração e ranqueamento de candidatos a biomarcadores de câncer de bexiga
   - Os genes são submetidos a um sistema de pontuação para determinar os melhores candidatos a biomarcadores
   * Cut-off de genes
     * |logFC| > 1
     * p-value ajustado por FDR < 0,05
     * I² < 50
     
   * Critérios:
     * Magnitude de expressão do gene (g de Hedges)
     * Significância estatística do gene (p-value ajustado por FDR)
     * Magnitude de expressão do gene pela análise de validação (g de Hedges)
     * Significância estatística do gene pela análise de validação (p-value ajustado por FDR)
     * AUC do gene pela análise de validação
     * Consistência de significância e presença de genes entre datasets
     * Heterogeneidade do gene (I²)
     * Grau de interações em rede PPI
     * Concordância de direção de expressão do gene entre descobrimento e validação
     * Direção de expressão do gene (up-regulados priorizados)

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
- affy (leitura e normalização)
- oligo (leitura e normalização)
- limma (análise estatística)
- edgeR (preparação para análise estatística voom)
- metafor (summarized mean difference)

### Anotação
- AnnotationDbi
- org.Hs.eg.db
- hgu133plus2.db
- hta20transcriptcluster.db
- hgu133a.db
- illuminaHumanv2.db
- illuminaHumanv4.db

### Meta-análise
- metafor (meta-análise)

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

## Organização do repositório


---

## Resultados principais



## Agradecimentos

Muito obrigado a minha família, amigos, colegas e especialmente a RioGen por todo o apoio durante essa jornada!

<p align="center">
  <img width="624" height="392" alt="equipe_riogen" src="https://github.com/user-attachments/assets/976a2956-573b-4764-a5ee-2f702238fd41" />
  <br>
  <em>RioGen - 15 de abril de 2026</em>
</p>

## Autor

Guilherme Moret Miranda - RioGen Tecnologia
