import subprocess
from tqdm import tqdm
from utils import load_bucket
import tempfile
import shutil
import os
from pathlib import Path

def agatify(gff_files: list[os.PathLike[str]], destination: os.PathLike[str]) -> None:
    extract_name = lambda accession: accession.split('.gff')[0]
    total_gff_files = len(gff_files)
    destination = Path(destination).resolve()
    destination.mkdir(exist_ok=True)

    for file in tqdm(gff_files, total=total_gff_files, leave=True, position=0):
        dest_name = extract_name(file) + ".agat.gff"

        tempdir = tempfile.TemporaryDirectory(prefix=extract_name(file) + ".")

        cur_dir = os.getcwd()
        os.chdir(tempdir.name)
        command = f"agat_convert_sp_gxf2gxf.pl -g {file} -o {dest_name}"
        subprocess.run(
                command,
                shell=True,
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
                )

        shutil.move(dest_name, destination)
        os.chdir(cur_dir)
        tempdir.cleanup()


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--bucket_id", type=int, default=0)
    parser.add_argument("--destination", type=str, default="gff_agatified")
    parser.add_argument("--schedule_path", type=str, default="schedule_100.json")
    
    args = parser.parse_args()
    destination = args.destination
    bucket_id = args.bucket_id
    schedule_path = args.schedule_path
    

    gff_files = load_bucket(bucket_id, schedule_path)
    agatify(gff_files, destination=destination)

