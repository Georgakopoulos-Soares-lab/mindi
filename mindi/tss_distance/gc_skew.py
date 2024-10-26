from Bio.SeqIO.FastaIO import SimpleFastaParser
from pybedtools import BedTool
import numpy as np
import pandas as pd

def parse_accession(accession):
    for record in SimpleFastaParser(accession):
        yield record[0].split(' ')[0], record[1].upper().strip()

if __name__ == "__main__":
    import sys
    oric  = sys.argv[1]
    genome = sys.argv[2]
    WINDOW_SIZE = 5_000
    GC_SKEW_WINDOW = 30
    oric = pd.read_table(oric,
                         header=None, 
                         comment="#",
                         # usecols=range(3),
                         usecols=[0, 3, 4],
                         names=["seqID", "start", "end"]
    )
    oric["start"] = (oric["start"]+oric["end"])//2
    oric["end"] = oric["start"] + WINDOW_SIZE + 1
    oric["start"] = np.maximum(oric["start"] - WINDOW_SIZE - GC_SKEW_WINDOW + 1, 0)
    oric_bed = BedTool.from_dataframe(oric)
    # extended_oric_bed = oric_bed.slop(
    #                                  l=WINDOW_SIZE + GC_SKEW_WINDOW, 
    #                                  r=WINDOW_SIZE, 
    #                                  g="genome.txt")
    oric_bed.sequence(genome, fo='ORIC_sequences.txt')
    # calculate GC skew in windows 
    gc_skew_all = []
    with open("ORIC_sequences.txt", "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                continue
            sequence = line 
            gc_skew_records = []
            for i in range(len(sequence)-GC_SKEW_WINDOW+1):
                chunk = sequence[i:i+GC_SKEW_WINDOW]
                g_content = chunk.count("G")
                c_content = chunk.count("C")
                gc_content = g_content + c_content
                if gc_content > 0:
                    gc_skew = (g_content - c_content)/(g_content + c_content) 
                else:
                    gc_skew = 0
                gc_skew_records.append(gc_skew)
            gc_skew_all.append(gc_skew_records)
    gc_skew = pd.DataFrame(gc_skew_all).sum()
    gc_skew_enriched = gc_skew / np.mean(gc_skew)
    gc_skew_enriched.index = gc_skew_enriched.index.map(lambda x: x - WINDOW_SIZE)
    gc_skew_enriched.to_csv("gc_skew_ORIC_chm13v2.0.txt", mode="w", index=True, header=None)
