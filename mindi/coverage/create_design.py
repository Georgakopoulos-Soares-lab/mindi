import argparse
from collections import defaultdict
from pathlib import Path
import pandas as pd
from termcolor import colored

def create_workflow_design(accession_list: str, extraction_path: str, suffix: str, tree: str) -> None:
    print(colored(f"Redirecting output to `design.csv`...", "green"))
    design = defaultdict(list)
    extract_id = lambda accession: '_'.join(Path(accession).name.split('_')[:2])
    extractions = {extract_id(file): file for file in extraction_path.glob(f"*{suffix}")}
    with open(accession_list, mode="r", encoding="UTF-8") as f:
        for line in f:
            line = line.strip()
            accession_path = Path(line).resolve()
            # gff_path = Path(line.replace("fna", "gff"))
            gff_path = accession_path.parent.joinpath(extract_id(accession_path) + ".gff").resolve()
            accession_id = extract_id(accession_path)
            if not gff_path.is_file():
                print(f"Skipping `{accession_path}` since it doesn't have a GFF file.")
                continue
            if accession_id not in extractions:
                print(f"Skipping `{accession_path}` since it doesn't have an extraction file.")
                continue
            extraction = extractions[accession_id]
            gff_path = str(gff_path)
            design["accession_id"].append(accession_id)
            design["accession"].append(accession_path)
            design["extraction"].append(extraction)
            design["gff"].append(gff_path)
    design = pd.DataFrame(design)
    design = merge_with_tree_of_life(design, tree)
    design.to_csv("design.csv", sep=",", mode="w", index=False, header=True)
    print(colored(f"Succesfully created `design.csv` file!", "green"))

def merge_with_tree_of_life(design, tree):
    tree_df = pd.read_table(tree)[["#assembly_accession", 
                                   "organism_name", 
                                   "group", 
                                   "genome_size", 
                                   "phylum", 
                                   "kingdom", 
                                   "superkingdom"]]\
                    .set_index("#assembly_accession")
    merged_design = design.merge(tree_df, 
                                 how="left", 
                                 left_on="accession_id",
                                 right_on="#assembly_accession",
                                 right_index=True
                                 )
    return merged_design
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Creates workflow design file.""")
    parser.add_argument("--accession", type=str, default="")
    parser.add_argument("--extraction", type=str, default="")
    parser.add_argument("--suffix", type=str, default=".txt")
    parser.add_argument("--tree", type=str, default="tree_of_life_new.txt.gz")
    args = parser.parse_args()
    accession_list = args.accession
    extraction_path = Path(args.extraction).resolve()
    tree = args.tree
    suffix = args.suffix
    create_workflow_design(accession_list, extraction_path, suffix=suffix, tree=tree)
