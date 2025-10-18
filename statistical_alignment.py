from Bio import SeqIO
from Bio.Seq import Seq
from Bio import motifs 
import pandas as pd
import random
import numpy as np
import argparse
import os
import ast
with open("upstream_seqs.txt", "r") as f:
    upstream_seqs_final = ast.literal_eval(f.read().strip())


parser = argparse.ArgumentParser(
    prog='ProgramName',
    description='Script to perform statistical alignment using genome and GFF files',
    epilog='Example: python statistical_alignment.py --genome_file genome.fna --gff_file annotations.gff'
)

parser.add_argument('--genome_file', required=True, help='path to genome FASTA file')
parser.add_argument('--gff_file', required=True, help='path to genome FASTA file')
parser.add_argument('--promoter_seqs', required=True, help='path to output file for normalized scores')



args = parser.parse_args()
genome_file = args.genome_file
genome_record = SeqIO.read(genome_file, "fasta")
genome_seq = str(genome_record.seq)

gff_file = args.gff_file
gff_cols = ["seqid","source","type","start","end","score","strand","phase","attributes"]
gff_df = pd.read_csv(gff_file, sep="\t", comment="#", names=gff_cols)

def extracting_bases(genome_seq, gff_df,upstream_seqs_final):
    genes_df = gff_df[(gff_df["source"] == "Protein Homology")]
    upstream_seqs = []

    for _, row in genes_df.iterrows():
        start = int(row["start"]) - 1 
        if row["strand"] == "+":  
            upstream_start = max(0, start - 15)
            upstream_end = max(0, start - 5)
            seq = genome_seq[upstream_start:upstream_end]
        else:  # antisense strand
            # upstream for antisense is downstream on genome
            upstream_start = start + 5  
            upstream_end = start + 15   
            seq = genome_seq[upstream_start:upstream_end]
            seq = str(Seq(seq).reverse_complement())  

        upstream_seqs.append(seq)
        
    
    promoter_seqs = []

    for seq in upstream_seqs_final[:100]:  # use only 1st 100 sequences
        for i in range(len(seq) - 5):
            window = seq[i:i+6]
            if all(base in "AT" for base in window)  and window[1]=='A' and window[-1]=='T':  # check if all 6 are W (A/T)
                promoter_seqs.append(window)


    return promoter_seqs

def create_ppm(promoter_seqs):
    promoter_motifs = motifs.create([Seq(i) for i in promoter_seqs])

    for i in range(6):
        for base in "ATCG":
            if promoter_motifs.counts[base][i] == 0:
                promoter_motifs.counts[base][i] = 0.001
    bases = 'ACGT'

    promoter_motifs_counts = [[float(promoter_motifs.counts[b][i]) for b in bases] for i in range(6)]

    for i in range(6):
        for j in range(4):
            promoter_motifs_counts[i][j] /= sum(promoter_motifs_counts[i])
            round(promoter_motifs_counts[i][j],3)

    ppm = promoter_motifs.counts.normalize()

    consensus = promoter_motifs.consensus
    consensus_score = ppm["A",0]*ppm["A",1]*ppm["A",2]*ppm["A",3]*ppm["A",4]*ppm["T",5]

    return ppm, consensus, consensus_score

def statistical_alignment(upstream_seqs_final, ppm, consensus_score):
    norm_score = {}

    for idx, seq in enumerate(upstream_seqs_final[100:]):
        seq_scores = {}
        for i in range(len(seq) - 5):
            window = seq[i:i+6]
            window_score = (ppm[window[0], 0] * ppm[window[1], 1] * ppm[window[2], 2] * ppm[window[3], 3] * ppm[window[4], 4] * ppm[window[5], 5])
            seq_scores[window] = np.log(window_score/consensus_score)
        norm_score[idx] = seq_scores

    for i in range(1000):
        for j in norm_score[i]:
            if norm_score[i][j] > -1.757020252640843:
                print(f"Sequence index: {i}, Window: {j}, Normalized Score: {norm_score[i][j]}")