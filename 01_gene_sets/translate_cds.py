#!/usr/bin/env python3
"""
translate_cds.py — Translate CDS FASTA to protein FASTA.
Keeps longest isoform per gene (splits on '.' to get gene ID).

Usage:
  python translate_cds.py input_cds.fa output_prot.fa output_cds_filtered.fa
"""
import sys
from Bio import SeqIO
from Bio.Seq import Seq

cds_path, prot_out, cds_out = sys.argv[1], sys.argv[2], sys.argv[3]

# Load all, keep longest per gene
genes = {}
for rec in SeqIO.parse(cds_path, "fasta"):
    gene_id = rec.id.split('.')[0]  # g694.t1 -> g694
    seq = str(rec.seq).upper()
    if gene_id not in genes or len(seq) > len(str(genes[gene_id].seq)):
        rec.id = gene_id
        rec.name = gene_id
        rec.description = gene_id
        genes[gene_id] = rec

print(f"Loaded {len(genes)} unique genes from CDS")

# Write filtered CDS and translated protein
prot_records = []
cds_records = []
skipped = 0
for gene_id, rec in genes.items():
    seq = str(rec.seq).upper()
    # Remove trailing stop codon if present
    if len(seq) % 3 != 0:
        # Trim to multiple of 3
        seq = seq[:len(seq) - (len(seq) % 3)]
    prot_seq = str(Seq(seq).translate()).rstrip('*')
    if len(prot_seq) < 30:
        skipped += 1
        continue
    from Bio.SeqRecord import SeqRecord
    prot_rec = SeqRecord(Seq(prot_seq), id=gene_id, name=gene_id, description=gene_id)
    prot_records.append(prot_rec)
    cds_records.append(rec)

SeqIO.write(prot_records, prot_out, "fasta")
SeqIO.write(cds_records, cds_out, "fasta")
print(f"Wrote {len(prot_records)} proteins ({skipped} skipped, <30 aa)")
