import json
import os
import random
from abc import abstractmethod

class Scheduler:
    
    # INTERFACE

    @abstractmethod
    def schedule(self, files: list[os.PathLike[str]], total_buckets: int) -> list[list[str]]:
        pass
    
    def saveas(self, scheduled_files: dict[str, list[str]], dest: os.PathLike[str]) -> None:
        with open(dest, mode="w", encoding="UTF-8") as f:
            json.dump(scheduled_files, f, indent=4)

class MiniBucketScheduler(Scheduler):

    def schedule(self, files: list[os.PathLike[str]], total_buckets: int) -> dict[str, list[str]]:
        scheduled_files = [[] for _ in range(total_buckets)]
        bucket_burden = [0 for _ in range(total_buckets)]
        
        print(f"Initializing scheduling for {len(files)} total files.")
        print(f"Assigning files to ⟶  {total_buckets} buckets.")

        for file in files:

            # pick minimum bucket
            file_size = os.path.getsize(file)
            
            # fetch bucket with minimum computational burden
            minimum_bucket_burden = min(bucket_burden)

            # fetch bucket array ids corresponding to minimum computational burden
            minimum_bucket_ids = [i for i in range(total_buckets) if bucket_burden[i] == minimum_bucket_burden]
            
            # fetch random id from the collection of minimum burden buckets
            random_mini_bucket_pos = random.choice(minimum_bucket_ids)

            # pick the bucket
            mini_bucket = scheduled_files[random_mini_bucket_pos]
        
            # append to it the file
            mini_bucket.append(file)

            # increase the corresponding bucket burden by the current files total byte size
            bucket_burden[random_mini_bucket_pos] += file_size


        print("Scheduling has been completed succesfully.")

        for i in range(total_buckets):
            files = scheduled_files[i]
            burden = sum(map(os.path.getsize, files))
            print(f"Bucket {i+1}: {len(files)} files; {bucket_burden[i] / 1000: .2f} predicted Kbytes; {burden / 1000: .2f} total Kbytes.")

        return {bucket_id: job for bucket_id, job in enumerate(scheduled_files)}

class MeanBucketScheduler(Scheduler):

    def schedule(self, files: list[os.PathLike[str]], total_buckets: int) -> list[list[str]]:
        raise NotImplementedYet()

if __name__ == "__main__":
        
    from pathlib import Path
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--accessions_path", type=str, default='filtered_assemblies.txt')
    parser.add_argument("--total_buckets", type=int, default=2)
    
    
    args = parser.parse_args()
    accessions_path = Path(args.accessions_path).resolve()
    total_buckets = args.total_buckets

    accessions = []
    with open(accessions_path, mode="r") as f:
        for line in f:
            line = line.strip()

            if line.count("\t") > 1:
                line = line.split("\t")[0]

            accessions.append(line)

    scheduler = MiniBucketScheduler()
    scheduled_files = scheduler.schedule(accessions, total_buckets=total_buckets)
    
    destination = f"new_schedule_{total_buckets}.json"
    scheduler.saveas(scheduled_files, destination)











        


        
        
