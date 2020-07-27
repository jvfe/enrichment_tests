# %%
import pandas as pd
from intermine import registry
from intermine.webservice import Service

# %%
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
gene_features = pd.read_table("../data/gene_query.tsv", names=fields_of_interest)
mirna_features = pd.read_table("../data/mirna_query.tsv", names=fields_of_interest2)
# %%
