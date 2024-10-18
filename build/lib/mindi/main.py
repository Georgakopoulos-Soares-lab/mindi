from mindi.mindi import MindiTool
import pandas as pd
from pathlib import Path
import os
import concurrent.futures 
from functools import reduce
import numpy as np

if __name__ == "__main__":

    def extract_mirrors(accessions: list[str], tempdir: str = "extractions_MR") -> pd.DataFrame:
        def _extract_mirrors(accession: str):
            hunter = MindiTool()
            hunter.set_tempdir(tempdir)
            mindi_table = hunter.extract_mirror(accession)
            mindi_table.validate(accession)
            df = mindi_table.to_dataframe()
            return mindi_table
        
        return [_extract_mirrors(accession) for accession in accessions]
#                 reduce(
#                      lambda acc, val: pd.concat([acc, _extract_mirrors(val)]), 
#                      accessions, 
#                      pd.DataFrame()
#                      )
#    
#
    indir = Path("accessions")
    extractions_dir = Path("extractions")
    extractions_dir.mkdir(exist_ok=True)

    accessions = [acc for acc in indir.glob("**/GC*.fna")]
    print(f"Total accessions detected {accessions}.")

    max_workers = 5
    jobs = [job.tolist() for job in np.array_split(accessions, max_workers)]
    print(f"Redirecting jobs to {max_workers} buckets!")

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = executor.map(extract_mirrors, jobs)

    for result in results:
        for df in result:
            df.moveto(extractions_dir)


    breakpoint()






