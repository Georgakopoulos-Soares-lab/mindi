from typing import Optional

def hunt_tail(sequence: str, minrepeat: int, seqID: Optional[str] = None) -> list[dict]:
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
                tandem_sequence = sequence[start-1: end]
                send_dict = dict()

                if seqID:
                    send_dict.update({"seqID": seqID})

                send_dict.update({
                      "start": start,
                      "end": end,
                      "sequence": tandem_sequence,
                      "length": end - start + 1,
                      "consensus": consensus,
                      "sru": len(consensus),
                      "consensus_repeats": tandem.count(consensus),
                    })

                yield send_dict

            consensus = nucleotide
            total_length = 1
            start = i+1


if __name__ == "__main__":

    import argparse

    default_sequence ="aaatttttttttaaaaaaattatatatatatatttatatataggcgcgggcggcggggtt"
    parser = argparse.ArgumentParser()
    parser.add_argument("--minrepeat", type=int, default=4)
    parser.add_argument("--sequence", type=str, default=default_sequence) 

    args = parser.parse_args()
    minrepeat = args.minrepeat
    sequence = args.sequence 


    for tandem_tail in hunt_tail(sequence=sequence, minrepeat=minrepeat):
        print(tandem_tail)
    breakpoint()

