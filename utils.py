from Bio import SeqIO
import gzip
import os
from pathlib import Path

def parse_fasta(accession: os.PathLike[str]) -> tuple[str]:
    accession = Path(accession).resolve()
    if accession.name.endswith(".gz"):
        file = gzip.open(accession, 'rt')
    else:
        file = open(accession, mode='r', encoding='UTF-8')
        
    for record in SeqIO.parse(file, 'fasta'):
        yield str(record.id), str(record.seq).lower()

    file.close()
