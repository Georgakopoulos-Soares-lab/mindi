import os
import csv
from abc import abstractmethod
import pandas as pd

class PWMExtractor:

    def __init__(self, mode: str) -> None:
        self.mode = mode
        self.nucleotides = "agctn"

    @staticmethod
    def invert(nucleotide: str) -> str:
        match nucleotide:
            case 'a':
                return 't'
            case 't':
                return 'a'
            case 'g':
                return 'c'
            case 'c':
                return 'g'
            case 'n':
                # this is used as a spacer variable;
                return 'n'
            case _:
                raise ValueError(f'Unknown nucleotide {nucleotide}.')

    def extract_RHS(self, sequence_of_arm: str) -> str:
        match self.mode:
            case 'DR':
                return sequence_of_arm
            case 'MR':
                return sequence_of_arm[::-1]
            case 'IR':
                return ''.join(PWMExtractor.invert(n) for n in sequence_of_arm)[::-1]
            case _:
                raise ValueError(f'Unknown mode {mode}.')
    
    def extract_PWM(self, intersect_df: pd.DataFrame, window_size: int, total_counts: bool = False) -> dict[str, list[int]]:
        total_counts = {n: [0 for _ in range(2*window_size+1)] for n in self.nucleotides}
        # total_counts = [0 for _ in range(2*window_size+1)]

        total_overlap = 0
        for _, row in intersect_df.iterrows():
                strand = row['strand']

                if strand == "?":
                    continue

                start = int(row['start'])
                end = int(row['end'])
                motif_start = int(row['motif_start'])
                motif_end = int(row['motif_end'])

                if self.mode == "STR":
                    sequence = row['sequence']
                    category = row['sru']
                elif self.mode == "DR" or self.mode == "MR" or self.mode == "IR":
                    sequence_of_arm = row['sequenceOfArm']
                    category = int(row['spacerLength'])
                    right_hand_side = self.extract_RHS(sequence_of_arm)
                    sequence = sequence_of_arm + 'n' * category + right_hand_side
                else:
                    raise ValueError(f"Invalid PWM mode {self.mode}.")


                if strand == "-":
                    sequence = "".join(PWMExtractor.invert(n) for n in sequence)[::-1]

                overlap = int(row['overlap']) 
                total_overlap += overlap
                origin = end - window_size - 1
                L = max(0, window_size - (origin - motif_start))
                U = min(2 * window_size + 1, window_size - (origin - motif_end))
                
                assert L <= U
                overlap_start = max(motif_start, start) 
                overlap_end = min(motif_end, end)
                overlap_length = overlap_end - overlap_start
                assert overlap == overlap_length

                overlapping_sequence = sequence[max(0, start-motif_start): min(end-motif_start, len(sequence))]
                assert len(overlapping_sequence) == overlap == U-L
                
                for idx, pos in enumerate(range(L, U)):
                    if strand == "-":
                        index = 2 * window_size - pos
                    else:
                        index = pos
                    
                    nucl = overlapping_sequence[idx]
                    total_counts[nucl][index] += 1
        
        assert total_overlap == sum(sum(v) for v in total_counts.values())
        return total_counts
