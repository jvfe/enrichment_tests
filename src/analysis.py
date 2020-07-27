# %%
import pandas as pd
from intermine import registry
from intermine.webservice import Service

# %%
######## Data acquistion #########
# Intermine setup
service = Service("https://targetmine.mizuguchilab.org/targetmine/service")
gene_end = service.new_query("Gene")
mirna_end = service.new_query("MiRNA")
registry.getData("targetmine")

# %%
fields_of_interest = [
    "Gene.organism.name",
    "Gene.symbol",
    "Gene.primaryIdentifier",
    "networkProperties.isBottleneck",
    "networkProperties.isHub",
    "pathways.identifier",
    "pathways.name",
    "miRNAInteractions.targetGene.symbol",
    "miRNAInteractions.supportType",
]
gene_end.add_constraint(
    "miRNAInteractions.supportType", "=", "Functional MTI", code="B"
)
gene_query = gene_end.select(fields_of_interest)

# %%
fields_of_interest2 = [
    "primaryIdentifier",
    "secondaryIdentifier",
    "organism.name",
    "miRNAInteractions.targetGene.primaryIdentifier",
    "miRNAInteractions.targetGene.symbol",
    "miRNAInteractions.targetGene.organism.shortName",
    "miRNAInteractions.supportType",
]
mirna_end.add_constraint(
    "miRNAInteractions.supportType", "=", "Functional MTI", code="B"
)
mirna_query = mirna_end.select(fields_of_interest2)

# %%
def save_query_table(queryobj, outfile):
    with open(outfile, "w") as f:
        for g in queryobj.results("tsv"):
            f.write(g + "\n")


# %%
save_query_table(gene_query, "../data/gene_query.tsv")
save_query_table(mirna_query, "../data/mirna_query.tsv")
# %%
######## Intersections #########
gene_features = pd.read_table("../data/gene_query.tsv", names=fields_of_interest).query(
    "`Gene.organism.name` == 'Homo sapiens'"
)

mirna_features = pd.read_table(
    "../data/mirna_query.tsv", names=fields_of_interest2
).query("`organism.name` == 'Homo sapiens'")

data_genes = pd.read_csv("../data/total_results_isa.csv")

# %%
gene_enriched = pd.merge(
    gene_features,
    data_genes,
    left_on=["Gene.symbol"],
    right_on=["gene_name"],
    how="inner",
).drop(
    [
        "Gene.organism.name",
        "Gene.symbol",
        "Rank",
        "miRNAInteractions.targetGene.symbol",
        "miRNAInteractions.supportType",
    ],
    axis=1,
)

mirna_enriched = pd.merge(
    mirna_features,
    data_genes,
    left_on=["miRNAInteractions.targetGene.symbol"],
    right_on=["gene_name"],
    how="inner",
).drop(
    [
        "miRNAInteractions.targetGene.organism.shortName",
        "miRNAInteractions.supportType",
        "organism.name",
        "Rank",
        "gene_name",
    ],
    axis=1,
)

# %%
gene_enriched.to_csv("../data/isa_pathway_enriched.csv", index=False)
mirna_enriched.to_csv("../data/isa_mirna_enriched.csv", index=False)
