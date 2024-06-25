from typing import Optional

def hunt_tail(sequence: str, minrepeat: int, chromosome: Optional[str] = None) -> list[dict]:
    start = 1
    consensus = sequence[0]
    rows = []
    total_length = 1
    nucleotides = {"a", "g", "c", "t"}

    for i, nucleotide in enumerate(sequence[1:] + "*", 1):
        if consensus == nucleotide and nucleotide in nucleotides:
            total_length += 1
        else:
            if total_length >= minrepeat:
                end = i
                tandem = sequence[start-1: end]
                rows.append([start,
                             end,
                             end - start + 1,
                             consensus,
                             tandem,
                             len(consensus),
                             tandem.count(consensus)
                             ]
                    )

            consensus = nucleotide
            total_length = 1
            start = i+1
    return rows


if __name__ == "__main__":

    import argparse

    default_sequence ="aaatttttttttaaaaaaattatatatatatatttatatataggcgcgggcggcggggtt"
    parser = argparse.ArgumentParser()
    parser.add_argument("--minrepeat", type=int, default=4)
    parser.add_argument("--sequence", type=str, default=default_sequence) 

    args = parser.parse_args()
    minrepeat = args.minrepeat
    sequence = args.sequence 


    tandem_tail = hunt_tail(sequence=sequence, minrepeat=minrepeat)
    breakpoint()

