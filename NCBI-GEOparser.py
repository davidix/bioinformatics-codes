import GEOparse
import pandas as pd
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import importr

# Download and parse data from GEO database
gse = GEOparse.get_GEO(geo="GSE123456")

# Extract gene expression data
df = gse.pivot_samples('VALUE')

# Convert gene expression data to R DataFrame
pandas2ri.activate()
r_df = pandas2ri.py2rpy(df)

# Perform differential gene expression analysis using DESeq2
deseq2 = importr('DESeq2')
dds = deseq2.DESeqDataSetFromMatrix(countData=r_df, colData=gse.phenotype_data, design=~gse.phenotype_data['condition'])
dds = deseq2.DESeq(dds)
res = deseq2.results(dds)

# Extract differentially expressed genes
de_genes = res[res.rx2('padj') < 0.05].index.tolist()

# Perform pathway analysis using KEGG
kegg = importr('KEGGREST')
kegg_genes = kegg.bconv(de_genes, 'id', 'ncbi-geneid', organism='hsa')
pathways = kegg.pathwayFromGene(kegg_genes, 'hsa')
