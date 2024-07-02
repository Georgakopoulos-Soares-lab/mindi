from pathlib import Path
import argparse
from termcolor import colored

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    
    parser.add_argument("--files_path", type=str, default="filtered_assemblies.txt")
    parser.add_argument("--extractions_path", type=str, default="repeat_out/IR_extracted_accessions")

    args = parser.parse_args()
    all_files_path = Path(args.files_path).resolve()
    extractions_path = Path(args.extractions_path).resolve()

    accessions = dict()
    extract_id = lambda accession: '_'.join(Path(accession).name.split("_")[:2])

    with open(all_files_path, mode="r") as f:

        for line in f:
            line = line.strip()

            if line.count("\t") > 0:
                line = line.split("\t")[0]
            
            accession_id = extract_id(line)
            accessions.update({accession_id: line})

    
    extractions_path = Path(extractions_path).resolve()
    extracted_accessions = {extract_id(file) for file in extractions_path.glob("*.csv")}
    print(f"Total accessions detected: {len(accessions)}.")
    print(colored(f"Total extracted accessions detected: {len(extracted_accessions)}.", "green"))
    
    total_remaining = 0
    with open("remaining_files.txt", mode="w") as f:
        for accession_id, acc in accessions.items():
            if accession_id not in extracted_accessions:
                f.write(str(acc) + "\n")
                total_remaining += 1

    print(colored(f"Total remaining accessions detected: {total_remaining}.", "red"))
