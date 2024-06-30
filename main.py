from mindi import MindiHunter
import os

if __name__ == "__main__":

    hunter = MindiHunter()
    hunter.set_tempdir("extractions_IR")
    
    accession = "accessions/GCF_000002515.2_ASM251v1_genomic.fna.gz"
    accession = "accessions/GCA_019677325.2_ASM1967732v2_genomic.fna.gz"
    accession = "accessions/GCF_023122875.1_ASM2312287v1_genomic.fna.gz"

    mindi_table = hunter.extract_inverted(accession)






