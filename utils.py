from Bio import SeqIO
import time
import gzip
import os
from pathlib import Path
import threading
import logging

class ProgressTracker:

    def __init__(self, total_accessions: int, log_filename: str = "biologs/tracker.log") -> None:
        self.log_filename = Path(log_filename).resolve()
        self.log_filename.parent.mkdir(exist_ok=True)

        self.sleeping_time = 600 # seconds sleeping time

        self.counter = 0
        self.total_accessions = total_accessions

        logging.basicConfig(
                            filename=self.log_filename,
                            level=logging.INFO, 
                            format="%(levelname)s:%(message)s"
                        )
    
    def _get_progress(self) -> float:
        return self.counter * 1e2 / self.total_accessions

    def track_progress(self) -> None:
        while True:
            logging.info(f"Progress level: {self._get_progress()}%.")
            time.sleep(self.sleeping_time)



def parse_fasta(accession: os.PathLike[str]) -> tuple[str]:
    accession = Path(accession).resolve()
    if accession.name.endswith(".gz"):
        file = gzip.open(accession, 'rt')
    else:
        file = open(accession, mode='r', encoding='UTF-8')
        
    for record in SeqIO.parse(file, 'fasta'):
        yield str(record.id), str(record.seq).lower()

    file.close()
