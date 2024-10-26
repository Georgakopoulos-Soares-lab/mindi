from mindi.mindi import MindiTool
import pandas as pd
from pathlib import Path
import concurrent.futures 
import numpy as np

if __name__ == "__main__":

    extract = lambda mode, mindi: {"MR": mindi.extract_MR,
                                   "IR": mindi.extract_IR,
                                   "G4": mindi.extract_regex,
                                   "DR": mindi.extract_DR,
                                   }[mode]

    def extract_repeats(accessions: list[str], mode: str, tempdir: str = "tmpdir") -> pd.DataFrame:
        def _extract_repeats(accession: str):
            hunter = MindiTool()
            hunter.set_tempdir(tempdir)
            extractor = extract(mode, hunter)
            mindi_table = extractor(accession)
            # mindi_table.validate(accession)
            df = mindi_table.to_dataframe()
            return mindi_table

        return [_extract_repeats(accession) for accession in accessions]
#                 reduce(
#                      lambda acc, val: pd.concat([acc, _extract_MRs(val)]), 
#                      accessions, 
#                      pd.DataFrame()
#                      )
#
    indir = Path("mindi/accessions")
    extractions_dir = Path("mindi/extractions")
    extractions_dir.mkdir(exist_ok=True)

    accessions = [acc for acc in indir.glob("**/GC*.fna")]
    print(f"Total accessions detected {accessions}.")

    max_workers = 5
    jobs = [job.tolist() for job in np.array_split(accessions, max_workers)]
    print(f"Redirecting jobs to {max_workers} buckets!")

    modes = ["MR", "IR", "G4"] #, "DR"]
    # modes = ["DR"]
    for mode in modes:
        extractions_dir = Path(f"mindi/extractions_{mode}")
        extractions_dir.mkdir(exist_ok=True)
        def repeats(accession):
            return extract_repeats(accession, mode=mode)

        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            results = executor.map(repeats, jobs)

        for result in results:
            for df in result:
                df.moveto(extractions_dir)
