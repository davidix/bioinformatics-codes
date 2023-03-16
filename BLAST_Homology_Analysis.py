from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import Phylo

# Set input parameters
query_sequence = "query.fasta"
database = "nr"
evalue_cutoff = 0.001
alignment_output = "alignment.fasta"
tree_output = "tree.nw"

# Perform BLAST search
result_handle = NCBIWWW.qblast("blastn", database, query_sequence, expect=evalue_cutoff, hitlist_size=1000)

# Parse BLAST results
blast_records = NCBIXML.parse(result_handle)
blast_record = next(blast_records)

# Extract homologous sequences
homologous_sequences = []
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < evalue_cutoff:
            homologous_sequences.append(alignment.title.split()[1])

# Download homologous sequences
homologous_records = []
for homologous_sequence in homologous_sequences:
    homologous_record = SeqIO.read(NCBIWWW.qblast("fasta", database, homologous_sequence), "fasta")
    homologous_records.append(homologous_record)

# Align sequences using Clustal Omega
ClustalOmegaCommandline("clustalo", infile=query_sequence, outfile=alignment_output)()

for homologous_record in homologous_records:
    ClustalOmegaCommandline("clustalo", infile=homologous_record.seq, outfile=alignment_output, append=True)()

# Generate phylogenetic tree using neighbor-joining method
alignment = AlignIO.read(alignment_output, "fasta")
calculator = DistanceCalculator('identity')
dm = calculator.get_distance(alignment)
constructor = DistanceTreeConstructor(calculator, 'nj')
tree = constructor.build_tree(alignment)
Phylo.write(tree, tree_output, "newick")
