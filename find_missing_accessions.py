from pathlib import Path
from termcolor import colored

if __name__ == "__main__":

    all_files_path = "filtered_assemblies.txt"
    extractions_path = "repeat_out/IR_extracted_accessions"

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
    extracted_accessions = {extract_id(file) for file in extractions_path.glob("*.IR.csv")}
    print(f"Total accessions detected: {len(accessions)}.")
    print(colored(f"Total extracted accessions detected: {len(extracted_accessions)}.", "green"))
    
    total_remaining = 0
    with open("remaining_files.txt", mode="w") as f:
        for accession_id, acc in accessions.items():
            if accession_id not in extracted_accessions:
                f.write(str(acc) + "\n")
                total_remaining += 1

    print(colored(f"Total remaining accessions detected: {total_remaining}.", "red"))








    
