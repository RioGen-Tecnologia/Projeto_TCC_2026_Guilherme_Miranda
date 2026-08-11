"""
PubMed Biomarker Literature Search
==================================

Ferramenta exploratória para auxiliar a interpretação dos genes candidatos
a biomarcadores de câncer de bexiga.

O script:
1. Lê o ranking de biomarcadores produzido pelo pipeline principal em R.
2. Seleciona os Top N genes.
3. Consulta o NCBI Gene usando o Entrez ID.
4. Recupera símbolo, aliases e designações alternativas.
5. Pesquisa esses termos no PubMed em associação com câncer de bexiga.
6. Realiza buscas específicas para categorias de interesse.
7. Remove artigos duplicados.
8. Salva os resultados em results/pubmed/.

Esta análise é exploratória e complementar à discussão do trabalho.
Ela não constitui uma etapa estatística ou metodológica do pipeline principal.
"""

from pathlib import Path
from typing import Any
import time
import re

import pandas as pd
# pyrefly: ignore [missing-import]
from Bio import Entrez


# =============================================================================
# CONFIGURAÇÕES
# =============================================================================

# E-mail exigido pelo NCBI para utilização da API Entrez.
ENTREZ_EMAIL = "guiguimoret@hotmail.com"

# Número de genes a serem analisados.
TOP_N_GENES = 20

# Número máximo de artigos retornados por cada consulta.
MAX_RESULTS_PER_QUERY = 100

# Pequena pausa entre requisições para respeitar as políticas do NCBI.
REQUEST_DELAY = 0.35

# Termo principal utilizado nas buscas.
CANCER_TERM = '"bladder cancer"'

# Diretório raiz do projeto.
PROJECT_ROOT = Path(__file__).resolve().parent.parent

# Arquivo com o ranking dos biomarcadores.
RANK_FILE = (
    PROJECT_ROOT
    / "results"
    / "biomarker_results"
    / "Bladder_cancer_biomarker_rank.csv"
)

# Diretório destinado exclusivamente aos resultados da busca bibliográfica.
OUTPUT_DIR = PROJECT_ROOT / "results" / "pubmed"


# =============================================================================
# CATEGORIAS DE BUSCA
# =============================================================================

# Estas categorias são utilizadas para encontrar artigos que possam ser
# particularmente interessantes para a discussão dos genes.
#
# Elas NÃO representam uma classificação científica definitiva do artigo.

SEARCH_CATEGORIES = {
    "biomarker": (
        '"biomarker" OR "biomarkers" OR '
        '"biological marker" OR "tumor marker"'
    ),
    "diagnostic": (
        '"diagnostic" OR "diagnosis" OR '
        '"diagnostic marker" OR "diagnostic biomarker"'
    ),
    "roc_auc": (
        '"ROC" OR "AUC" OR "area under the curve" '
        'OR "receiver operating characteristic"'
    ),
    "urine": (
        '"urine" OR "urinary" OR "urine sample"'
    ),
    "prognostic": (
        '"prognostic" OR "prognosis" OR '
        '"survival" OR "overall survival"'
    ),
    "validation": (
        '"validation" OR "validated" OR '
        '"validation cohort" OR "independent cohort"'
    ),
    "expression": (
        '"expression" OR "differential expression" OR '
        '"gene expression"'
    ),
}


# =============================================================================
# FUNÇÕES AUXILIARES
# =============================================================================

def clean_text(value: Any) -> str:
    """Converte um valor para texto limpo."""
    if value is None:
        return ""

    text = str(value).strip()

    if text.lower() in {"nan", "none", "null"}:
        return ""

    return text


def unique_preserve_order(values: list[str]) -> list[str]:
    """Remove duplicatas mantendo a ordem original."""
    seen = set()
    result = []

    for value in values:
        value = clean_text(value)

        if not value:
            continue

        key = value.lower()

        if key not in seen:
            seen.add(key)
            result.append(value)

    return result


def normalize_search_term(term: str) -> str:
    """
    Prepara um termo para busca no PubMed.

    Mantém o termo entre aspas para favorecer busca como frase,
    especialmente para designações compostas.
    """
    term = clean_text(term)

    if not term:
        return ""

    term = term.replace('"', "")

    return f'"{term}"'


def pause():
    """Pausa entre requisições ao NCBI."""
    time.sleep(REQUEST_DELAY)


# =============================================================================
# NCBI GENE
# =============================================================================

def get_gene_information(entrez_id: str) -> dict[str, Any]:
    """
    Consulta o NCBI Gene utilizando o Entrez ID.

    Retorna:
    - símbolo oficial;
    - aliases;
    - outras designações;
    - nome de nomenclatura.
    """

    try:
        handle = Entrez.esummary(
            db="gene",
            id=entrez_id,
            retmode="xml",
        )

        record = Entrez.read(handle)
        handle.close()

        pause()

        summaries = record["DocumentSummarySet"]["DocumentSummary"]

        if not summaries:
            return {
                "symbol": "",
                "aliases": [],
                "designations": [],
                "nomenclature": "",
            }

        summary = summaries[0]

        symbol = clean_text(summary.get("Name"))
        nomenclature = clean_text(summary.get("NomenclatureName"))

        aliases_raw = clean_text(summary.get("OtherAliases"))
        designations_raw = clean_text(summary.get("OtherDesignations"))

        aliases = []

        if aliases_raw:
            aliases = [
                alias.strip()
                for alias in aliases_raw.split(",")
                if alias.strip()
            ]

        designations = []

        if designations_raw:
            designations = [
                designation.strip()
                for designation in designations_raw.split("|")
                if designation.strip()
            ]

        return {
            "symbol": symbol,
            "aliases": unique_preserve_order(aliases),
            "designations": unique_preserve_order(designations),
            "nomenclature": nomenclature,
        }

    except Exception as error:
        print(f"  Erro ao consultar NCBI Gene: {error}")

        return {
            "symbol": "",
            "aliases": [],
            "designations": [],
            "nomenclature": "",
        }


# =============================================================================
# PUBMED
# =============================================================================

def pubmed_search(query: str) -> list[str]:
    """
    Executa uma busca no PubMed e retorna os PMIDs encontrados.
    """

    try:
        handle = Entrez.esearch(
            db="pubmed",
            term=query,
            retmax=MAX_RESULTS_PER_QUERY,
            sort="relevance",
            retmode="xml",
        )

        record = Entrez.read(handle)
        handle.close()

        pause()

        return [str(pmid) for pmid in record["IdList"]]

    except Exception as error:
        print(f"  Erro na busca PubMed: {error}")
        return []


def build_general_query(gene_term: str) -> str:
    """Cria a busca geral gene + câncer de bexiga."""

    return (
        f"{normalize_search_term(gene_term)}[Title/Abstract] "
        f"AND {CANCER_TERM}[Title/Abstract]"
    )


def build_category_query(gene_term: str, category_terms: str) -> str:
    """Cria uma busca gene + câncer de bexiga + categoria."""

    return (
        f"{normalize_search_term(gene_term)}[Title/Abstract] "
        f"AND {CANCER_TERM}[Title/Abstract] "
        f"AND ({category_terms})[Title/Abstract]"
    )


# =============================================================================
# BUSCA DE UM GENE
# =============================================================================

def search_gene(
    gene: str,
    entrez_id: str,
    rank: int,
    score: float,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """
    Pesquisa um gene no PubMed.

    Retorna:
    - lista de registros de artigos;
    - lista de registros das buscas realizadas.
    """

    print()
    print("=" * 70)
    print(f"Gene: {gene}")
    print(f"Rank: {rank}")
    print(f"Score: {score:.4f}")
    print("=" * 70)

    # -------------------------------------------------------------------------
    # Recupera informações adicionais pelo Entrez ID
    # -------------------------------------------------------------------------

    print(f"  Consultando NCBI Gene para Entrez ID {entrez_id}...")

    gene_info = get_gene_information(entrez_id)

    ncbi_symbol = gene_info["symbol"]
    aliases = gene_info["aliases"]
    designations = gene_info["designations"]
    nomenclature = gene_info["nomenclature"]

    if ncbi_symbol:
        print(f"  Símbolo NCBI: {ncbi_symbol}")

    if aliases:
        print(f"  Aliases: {', '.join(aliases)}")

    if designations:
        print(f"  Designações: {len(designations)}")

    # -------------------------------------------------------------------------
    # Monta conjunto de termos para busca
    # -------------------------------------------------------------------------

    search_terms = unique_preserve_order(
        [gene, ncbi_symbol] + aliases + designations
    )

    # -------------------------------------------------------------------------
    # Busca geral
    #
    # Utilizamos todos os termos alternativos, mas somente para a busca geral.
    # -------------------------------------------------------------------------

    article_records: dict[str, dict[str, Any]] = {}
    search_records = []

    for term in search_terms:

        query = build_general_query(term)

        print(f'  Geral: {query}')

        pmids = pubmed_search(query)

        print(f"    Artigos encontrados: {len(pmids)}")

        search_records.append({
            "Rank": rank,
            "Gene_Symbol": gene,
            "Entrez": entrez_id,
            "Biomarker_Score": score,
            "Search_Term": term,
            "Search_Type": "general",
            "Category": "",
            "Query": query,
            "Articles_Found": len(pmids),
        })

        for pmid in pmids:

            if pmid not in article_records:
                article_records[pmid] = {
                    "PMID": pmid,
                    "Rank": rank,
                    "Gene_Symbol": gene,
                    "Entrez": entrez_id,
                    "Biomarker_Score": score,
                    "Matched_Terms": set(),
                    "Search_Types": set(),
                    "Categories": set(),
                }

            article_records[pmid]["Matched_Terms"].add(term)
            article_records[pmid]["Search_Types"].add("general")

    # -------------------------------------------------------------------------
    # Buscas específicas
    #
    # Para evitar uma explosão de consultas, as buscas especializadas são
    # realizadas principalmente utilizando o símbolo oficial do gene.
    # -------------------------------------------------------------------------

    primary_term = ncbi_symbol or gene

    for category, category_terms in SEARCH_CATEGORIES.items():

        query = build_category_query(
            primary_term,
            category_terms,
        )

        print(f"  {category}: {query}")

        pmids = pubmed_search(query)

        print(f"    Artigos encontrados: {len(pmids)}")

        search_records.append({
            "Rank": rank,
            "Gene_Symbol": gene,
            "Entrez": entrez_id,
            "Biomarker_Score": score,
            "Search_Term": primary_term,
            "Search_Type": "category",
            "Category": category,
            "Query": query,
            "Articles_Found": len(pmids),
        })

        for pmid in pmids:

            if pmid not in article_records:
                article_records[pmid] = {
                    "PMID": pmid,
                    "Rank": rank,
                    "Gene_Symbol": gene,
                    "Entrez": entrez_id,
                    "Biomarker_Score": score,
                    "Matched_Terms": set(),
                    "Search_Types": set(),
                    "Categories": set(),
                }

            article_records[pmid]["Matched_Terms"].add(primary_term)
            article_records[pmid]["Search_Types"].add("category")
            article_records[pmid]["Categories"].add(category)

    # -------------------------------------------------------------------------
    # Converte conjuntos em texto para exportação
    # -------------------------------------------------------------------------

    records = []

    for record in article_records.values():

        records.append({
            **record,
            "Matched_Terms": "; ".join(
                sorted(record["Matched_Terms"])
            ),
            "Search_Types": "; ".join(
                sorted(record["Search_Types"])
            ),
            "Categories": "; ".join(
                sorted(record["Categories"])
            ),
            "NCBI_Symbol": ncbi_symbol,
            "NCBI_Aliases": "; ".join(aliases),
            "NCBI_Designations": "; ".join(designations),
            "NCBI_Nomenclature": nomenclature,
        })

    return records, search_records


# =============================================================================
# PUBMED ARTICLE DETAILS
# =============================================================================

def extract_abstract(article: Any) -> str:
    """Extrai o abstract completo de um registro PubMed."""

    abstract = article.get("MedlineCitation", {}).get("Article", {}).get(
        "Abstract"
    )

    if not abstract:
        return ""

    parts = []

    for item in abstract.get("AbstractText", []):

        text = str(item)

        label = ""

        try:
            label = item.attributes.get("Label", "")
        except AttributeError:
            pass

        if label:
            parts.append(f"{label}: {text}")
        else:
            parts.append(text)

    return " ".join(parts)


def extract_article_data(article: Any) -> dict[str, Any]:
    """Extrai informações bibliográficas de um artigo PubMed."""

    citation = article.get("MedlineCitation", {})
    article_data = citation.get("Article", {})

    pmid = clean_text(citation.get("PMID"))

    title = clean_text(
        article_data.get("ArticleTitle")
    )

    journal_data = article_data.get("Journal", {})

    journal = clean_text(
        journal_data.get("Title")
    )

    year = ""

    journal_issue = journal_data.get("JournalIssue", {})

    pub_date = journal_issue.get("PubDate", {})

    if pub_date:
        year = clean_text(
            pub_date.get("Year")
        )

        if not year:
            medline_date = clean_text(
                pub_date.get("MedlineDate")
            )

            match = re.search(r"\b(19|20)\d{2}\b", medline_date)

            if match:
                year = match.group(0)

    authors = []

    author_list = article_data.get("AuthorList", [])

    for author in author_list:

        last_name = clean_text(
            author.get("LastName")
        )

        initials = clean_text(
            author.get("Initials")
        )

        if last_name:
            if initials:
                authors.append(
                    f"{last_name} {initials}"
                )
            else:
                authors.append(last_name)

    doi = ""

    article_id_list = article.get("PubmedData", {}).get(
        "ArticleIdList",
        []
    )

    for article_id in article_id_list:

        try:
            if article_id.attributes.get("IdType") == "doi":
                doi = clean_text(article_id)
                break
        except AttributeError:
            continue

    abstract = extract_abstract(article)

    return {
        "PMID": pmid,
        "Title": title,
        "Year": year,
        "Journal": journal,
        "Authors": "; ".join(authors),
        "DOI": doi,
        "Abstract": abstract,
    }


def fetch_articles(pmids: list[str]) -> list[dict[str, Any]]:
    """
    Recupera dados bibliográficos dos artigos em lote.
    """

    if not pmids:
        return []

    print()
    print("=" * 70)
    print(f"Recuperando informações de {len(pmids)} artigos...")
    print("=" * 70)

    try:
        handle = Entrez.efetch(
            db="pubmed",
            id=",".join(pmids),
            retmode="xml",
        )

        records = Entrez.read(handle)
        handle.close()

        pause()

        articles = []

        for article in records["PubmedArticle"]:
            articles.append(
                extract_article_data(article)
            )

        return articles

    except Exception as error:
        print(f"Erro ao recuperar artigos: {error}")
        return []


# =============================================================================
# CLASSIFICAÇÃO SIMPLES POR TEXTO
# =============================================================================

def text_contains_any(text: str, terms: list[str]) -> bool:
    """Verifica se o texto contém pelo menos um dos termos."""

    text = text.lower()

    return any(
        term.lower() in text
        for term in terms
    )


def classify_article(
    title: str,
    abstract: str,
) -> dict[str, bool]:
    """
    Identifica termos/categorias presentes no título ou abstract.

    Isso é apenas uma ferramenta de triagem.
    Não representa uma classificação científica definitiva.
    """

    text = f"{title} {abstract}".lower()

    return {
        "Biomarker": text_contains_any(
            text,
            [
                "biomarker",
                "biomarkers",
                "biological marker",
                "tumor marker",
            ],
        ),

        "Diagnostic": text_contains_any(
            text,
            [
                "diagnostic",
                "diagnosis",
                "diagnostic marker",
                "diagnostic biomarker",
            ],
        ),

        "ROC_AUC": text_contains_any(
            text,
            [
                "roc",
                "auc",
                "area under the curve",
                "receiver operating characteristic",
            ],
        ),

        "Urine": text_contains_any(
            text,
            [
                "urine",
                "urinary",
                "urine sample",
            ],
        ),

        "Prognostic": text_contains_any(
            text,
            [
                "prognostic",
                "prognosis",
                "survival",
                "overall survival",
            ],
        ),

        "Validation": text_contains_any(
            text,
            [
                "validation",
                "validated",
                "validation cohort",
                "independent cohort",
            ],
        ),

        "Expression": text_contains_any(
            text,
            [
                "expression",
                "differential expression",
                "gene expression",
            ],
        ),
    }


# =============================================================================
# MAIN
# =============================================================================

def main():

    print()
    print("=" * 70)
    print("PUBMED BIOMARKER LITERATURE SEARCH")
    print("=" * 70)

    # -------------------------------------------------------------------------
    # Configuração do Entrez
    # -------------------------------------------------------------------------

    Entrez.email = ENTREZ_EMAIL

    # -------------------------------------------------------------------------
    # Verificação do arquivo de entrada
    # -------------------------------------------------------------------------

    if not RANK_FILE.exists():

        raise FileNotFoundError(
            f"\nArquivo de ranking não encontrado:\n{RANK_FILE}\n"
        )

    # -------------------------------------------------------------------------
    # Cria diretório de saída
    # -------------------------------------------------------------------------

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    print()
    print(f"Arquivo de ranking:")
    print(f"  {RANK_FILE}")

    print()
    print(f"Diretório de saída:")
    print(f"  {OUTPUT_DIR}")

    # -------------------------------------------------------------------------
    # Carrega ranking
    # -------------------------------------------------------------------------

    ranking = pd.read_csv(RANK_FILE)

    required_columns = {
        "Gene_Symbol",
        "Entrez",
        "biomarker_score",
    }

    missing_columns = required_columns - set(ranking.columns)

    if missing_columns:

        raise ValueError(
            "Colunas obrigatórias ausentes no arquivo de ranking: "
            + ", ".join(sorted(missing_columns))
        )

    ranking = ranking.head(TOP_N_GENES).copy()

    print()
    print(f"Genes selecionados: {len(ranking)}")

    print()
    print(
        ranking[
            [
                "Gene_Symbol",
                "biomarker_score",
            ]
        ].to_string(index=False)
    )

    # -------------------------------------------------------------------------
    # Busca dos genes
    # -------------------------------------------------------------------------

    all_article_records = []
    all_search_records = []

    for _, row in ranking.iterrows():

        gene = clean_text(row["Gene_Symbol"])
        entrez_id = clean_text(row["Entrez"])

        try:
            score = float(row["biomarker_score"])
        except (TypeError, ValueError):
            score = 0.0

        rank = int(row.name) + 1

        article_records, search_records = search_gene(
            gene=gene,
            entrez_id=entrez_id,
            rank=rank,
            score=score,
        )

        all_article_records.extend(article_records)
        all_search_records.extend(search_records)

    # -------------------------------------------------------------------------
    # Remove artigos duplicados
    # -------------------------------------------------------------------------

    article_df = pd.DataFrame(all_article_records)

    if article_df.empty:

        print()
        print("Nenhum artigo foi recuperado.")

        return

    article_df = (
        article_df
        .drop_duplicates(
            subset=[
                "Gene_Symbol",
                "PMID",
            ]
        )
        .reset_index(drop=True)
    )

    # -------------------------------------------------------------------------
    # Recupera dados bibliográficos
    # -------------------------------------------------------------------------

    pmids = article_df["PMID"].dropna().unique().tolist()

    article_details = fetch_articles(pmids)

    details_df = pd.DataFrame(article_details)

    if not details_df.empty:

        article_df = article_df.merge(
            details_df,
            on="PMID",
            how="left",
        )

    # -------------------------------------------------------------------------
    # Classificação textual auxiliar
    # -------------------------------------------------------------------------

    classification_records = []

    for _, row in article_df.iterrows():

        classification = classify_article(
            clean_text(row.get("Title", "")),
            clean_text(row.get("Abstract", "")),
        )

        classification_records.append(
            classification
        )

    classification_df = pd.DataFrame(
        classification_records
    )

    article_df = pd.concat(
        [
            article_df.reset_index(drop=True),
            classification_df.reset_index(drop=True),
        ],
        axis=1,
    )

    # -------------------------------------------------------------------------
    # Ordenação das colunas
    # -------------------------------------------------------------------------

    preferred_columns = [
        "Rank",
        "Gene_Symbol",
        "Entrez",
        "Biomarker_Score",
        "PMID",
        "Year",
        "Title",
        "Journal",
        "Authors",
        "DOI",
        "Matched_Terms",
        "Search_Types",
        "Categories",
        "Biomarker",
        "Diagnostic",
        "ROC_AUC",
        "Urine",
        "Prognostic",
        "Validation",
        "Expression",
        "NCBI_Symbol",
        "NCBI_Aliases",
        "NCBI_Designations",
        "NCBI_Nomenclature",
        "Abstract",
    ]

    existing_columns = [
        column
        for column in preferred_columns
        if column in article_df.columns
    ]

    remaining_columns = [
        column
        for column in article_df.columns
        if column not in existing_columns
    ]

    article_df = article_df[
        existing_columns + remaining_columns
    ]

    # -------------------------------------------------------------------------
    # Salva artigos
    # -------------------------------------------------------------------------

    articles_file = (
        OUTPUT_DIR
        / "pubmed_top20_articles.csv"
    )

    article_df.to_csv(
        articles_file,
        index=False,
    )

    # -------------------------------------------------------------------------
    # Resumo por gene
    # -------------------------------------------------------------------------

    summary_records = []

    for _, row in ranking.iterrows():

        gene = clean_text(row["Gene_Symbol"])

        gene_articles = article_df[
            article_df["Gene_Symbol"] == gene
        ]

        summary_records.append({
            "Rank": int(row.name) + 1,
            "Gene_Symbol": gene,
            "Entrez": clean_text(row["Entrez"]),
            "Biomarker_Score": float(row["biomarker_score"]),

            "Total_Articles": len(gene_articles),

            "Biomarker_Articles": int(
                gene_articles["Biomarker"].sum()
            ),

            "Diagnostic_Articles": int(
                gene_articles["Diagnostic"].sum()
            ),

            "ROC_AUC_Articles": int(
                gene_articles["ROC_AUC"].sum()
            ),

            "Urine_Articles": int(
                gene_articles["Urine"].sum()
            ),

            "Prognostic_Articles": int(
                gene_articles["Prognostic"].sum()
            ),

            "Validation_Articles": int(
                gene_articles["Validation"].sum()
            ),

            "Expression_Articles": int(
                gene_articles["Expression"].sum()
            ),
        })

    summary_df = pd.DataFrame(summary_records)

    summary_file = (
        OUTPUT_DIR
        / "pubmed_top20_gene_summary.csv"
    )

    summary_df.to_csv(
        summary_file,
        index=False,
    )

    # -------------------------------------------------------------------------
    # Log das buscas
    # -------------------------------------------------------------------------

    search_log_df = pd.DataFrame(
        all_search_records
    )

    search_log_file = (
        OUTPUT_DIR
        / "pubmed_top20_search_log.csv"
    )

    search_log_df.to_csv(
        search_log_file,
        index=False,
    )

    # -------------------------------------------------------------------------
    # Relatório final
    # -------------------------------------------------------------------------

    print()
    print("=" * 70)
    print("BUSCA CONCLUÍDA")
    print("=" * 70)

    print(
        f"Genes pesquisados: {len(ranking)}"
    )

    print(
        f"Artigos únicos recuperados: "
        f"{len(article_df)}"
    )

    print()
    print("Artigos:")
    print(f"  {articles_file}")

    print()
    print("Resumo:")
    print(f"  {summary_file}")

    print()
    print("Log das buscas:")
    print(f"  {search_log_file}")

    print()
    print("Resumo por gene:")
    print()

    print(
        summary_df.to_string(index=False)
    )

    print()
    print("=" * 70)


if __name__ == "__main__":
    main()