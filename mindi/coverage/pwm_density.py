import numpy as np
import pandas as pd
from typing import Callable
from abc import abstractmethod 

class StrandEvaluator:
    
    @abstractmethod
    def determine_strand(self, motif: str) -> str:
        raise NotImplementError()

class HDNAStrandEvaluator(StrandEvaluator):
    
    def determine_strand(self, motif: str) -> str:
        ag = motif.count("a") + motif.count("g")
        ct = len(motif) - ag
        if ag >= ct:
            return "+"
        return "-"
        
class GQuadruplexStrandEvaluator(StrandEvaluator):
    
    def determine_strand(self, motif: str) -> str:
        return "+" if motif.count("g") >= motif.count("c") else "-"

class GCRichStrandEvaluator(StrandEvaluator):

    def determine_strand(self, motif: str) -> str:
        gc = motif.count("g") + motif.count("c")
        ta = len(motif) - gc 
        if gc >= ta:
            return "+"
        return "-"

strand_evaluators = {
                     "HDNA": HDNAStrandEvaluator().determine_strand,
                     "G4": GQuadruplexStrandEvaluator().determine_strand,
                     "GCRich": GCRichStrandEvaluator().determine_strand
                    }

class PWMExtractor:

    def __init__(self) -> None:
        self.nucleotides = "agct"

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
            case _:
                raise ValueError(f'Unknown nucleotide {nucleotide}.')

    def determine_template_strand(self, compartment_strand: str, motif_strand: str) -> str:
        if compartment_strand == motif_strand:
            return "non_template"
        return "template"

    def extract_template_density(self, intersect_df: pd.DataFrame, 
                                        window_size: int, 
                                        determine_strand: Callable[[str], str],
                                        ) -> dict[str, list[int]]:
        total_counts = {
                        "template": np.zeros(2*window_size+1),
                        "non_template": np.zeros(2*window_size+1)
                    }

        total_overlap = 0
        for _, row in intersect_df.iterrows():
                compartment_strand = row['strand']
                if compartment_strand == "?":
                    continue
                start = int(row['start'])
                end = int(row['end'])
                motif_start = int(row['motif_start'])
                motif_end = int(row['motif_end'])
                motif_strand = determine_strand(row['sequence'])

                template = self.determine_template_strand(compartment_strand, motif_strand)
                temp_counts = np.zeros(2*window_size+1)

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

                temp_counts[L:U] += 1
                if compartment_strand == "-":
                    temp_counts = temp_counts[::-1]
                total_counts[template] += temp_counts
        assert total_overlap == np.sum(total_counts["template"]) + np.sum(total_counts["non_template"]), f"Overlap: {total_overlap} vs. Calculated overlap {np.sum(total_counts)}."
        return total_counts

    @staticmethod
    def bayes_estimator(counts: int, total_counts: int, total_bins: int) -> float:
        return (counts+1)/(total_counts+total_bins)

    @staticmethod
    def expected_value(counts: int, total_counts: int, total_bins: int) -> float:
        return PWMExtractor.bayes_estimator(counts, total_counts, total_bins)

    @staticmethod
    def enrichment(counts: int, total_counts: int, total_bins: int) -> float:
        return counts/(total_counts * PWMExtractor.bayes_estimator(counts, total_counts, total_bins))

    def extract_density(self, intersect_df: pd.DataFrame, window_size: int) -> list[int]:
        total_counts = np.zeros(2*window_size+1)
        total_overlap = 0
        for _, row in intersect_df.iterrows():
                compartment_strand = row['strand']
                if compartment_strand == "?":
                    continue
                start = int(row['start'])
                end = int(row['end'])
                motif_start = int(row['motif_start'])
                motif_end = int(row['motif_end'])

                temp_counts = np.zeros(2*window_size+1)
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
                assert overlap == U-L
                temp_counts[L:U] += 1

                if compartment_strand == "-":
                    temp_counts = temp_counts[::-1]
                total_counts += temp_counts
        total_sum = int(np.sum(total_counts))
        assert total_overlap == total_sum, f"Overlap: {total_overlap} vs. Calculated overlap {total_sum}."
        total_counts = list(total_counts)
        return total_counts
        # enrichments = [PWMExtractor.enrichment(counts, total_sum, 2 * window_size+1) for counts in total_counts]
        return enrichments

    def extract_PWM(self, 
                    intersect_df: pd.DataFrame, 
                    window_size: int, 
                    ) -> dict[str, list[int]]:
        total_counts = {n: [0 for _ in range(2*window_size+1)] for n in self.nucleotides}
        total_overlap = 0
        for _, row in intersect_df.iterrows():
                compartment_strand = row['strand']
                if compartment_strand == "?":
                    continue
                start = int(row['start'])
                end = int(row['end'])
                motif_start = int(row['motif_start'])
                motif_end = int(row['motif_end'])
                sequence = row['sequence']
                if compartment_strand == "-":
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
                    if compartment_strand == "-":
                        index = 2 * window_size - pos
                    else:
                        index = pos

                    nucl = overlapping_sequence[idx]
                    total_counts[nucl][index] += 1

        assert total_overlap == sum(sum(v) for v in total_counts.values())
        return total_counts
