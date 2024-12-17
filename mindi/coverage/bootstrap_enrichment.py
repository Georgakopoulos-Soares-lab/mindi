import pandas as pd 
from collections import defaultdict
from dataclasses import dataclass, field
from tqdm import tqdm
from scipy.stats import ks_2samp

@dataclass(frozen=True)
class params:
    N: int = field(default=1000)
    window_size: int = field(default=500)
    alpha: float = field(default=0.05)

def bootstrap_enrichment(design: str, 
                         enrichment_file: str,
                         domain: str,
                         params,
                         output: tuple[str, str]) -> None:
        def bootstrap(df: pd.DataFrame, N: int = 1_000, alpha: float = 0.05) -> tuple:
            bootstrapped_samples = []
            # for _ in range(N):
            for _ in tqdm(range(N), leave=True, position=0):
                sample_df = df.sample(frac=1.0, replace=True)
                sample_df = sample_df.sum(axis=0)
                mean = sample_df.mean()
                resampled_sum = sample_df / mean
                bootstrapped_samples.append(resampled_sum)
            bootstrapped_samples = pd.DataFrame(bootstrapped_samples)
            average = bootstrapped_samples.mean(axis=0)
            lower_ci = bootstrapped_samples.quantile((1-alpha)/2)
            upper_ci = bootstrapped_samples.quantile((1+alpha)/2)
            return average, lower_ci, upper_ci

        design_df = pd.read_csv(design, usecols=["accession_id", 
                                                 "group", 
                                                 "phylum", 
                                                 "kingdom", 
                                                 "superkingdom"])
        enrichment_df = pd.read_parquet(enrichment_file, engine="fastparquet")\
                        .merge(
                                design_df,
                                right_on="accession_id",
                                left_on="#assembly_accession",
                                how="inner"
                              )\
                        .query(f"superkingdom == '{domain}'")
        if enrichment_df.shape[0] == 0:
            raise ValueError(f"Empty dataframe for domain `{domain}`.")
        
        combinations = [
                        ("Occurrences_template", "protein_coding"),
                         ("Occurrences_template", "non_coding"),
                         ("Occurrences_non_template", "protein_coding"),
                         ("Occurrences_non_template", "non_coding")
                        ]

        ## Domain level bootstrap
        confidence_intervals = defaultdict(list)
        for comb in combinations: 
            typ, biotype = comb
            temp_df = enrichment_df[(enrichment_df["template|non_template"] == typ) \
                                & (enrichment_df["biotype"] == biotype)]\
                        [[str(i) for i in range(-params.window_size, params.window_size+1)]]
            average, lower_ci, upper_ci = bootstrap(temp_df, N=params.N, alpha=params.alpha)
            confidence_intervals[f"average_{typ}_{biotype}_{domain}"] = average
            confidence_intervals[f"lower_ci_{typ}_{biotype}_{domain}"] = lower_ci
            confidence_intervals[f"upper_ci_{typ}_{biotype}_{domain}"] = upper_ci

        #        for biotype in ["protein_coding", "non_coding"]:
        #            avg_template = confidence_intervals[f"average_Occurrences_template_{biotype}_{domain}"]
        #            avg_non_template = confidence_intervals[f"average_Occurrences_non_template_{biotype}_{domain}"]
        #            significance = ks_2samp(avg_template, avg_non_template, alternative='two-sided')
        #

        confidence_intervals = pd.DataFrame(confidence_intervals).T
        confidence_intervals.to_csv(output[0], sep=",", mode="w", header=True, index=True)

        ## Phylum calculation
        ## keep assembly accessions with a present phylum
        enrichment_df_phylum = enrichment_df.dropna(subset=['phylum'])
        phylum_to_kingdom = enrichment_df_phylum[["phylum", "kingdom"]].set_index("phylum")["kingdom"].to_dict()
        enrichment_df_phylum = enrichment_df_phylum.groupby(["phylum", "template|non_template", "biotype"])\
            .agg({str(i): "sum" \
                      for i in range(-params.window_size, params.window_size+1)
                })
        enrichment_df_phylum = enrichment_df_phylum.apply(lambda row: row / row.mean(), axis=1)
        enrichment_df_phylum.reset_index(inplace=True)
        enrichment_df_phylum["kingdom"] = enrichment_df_phylum["phylum"].map(phylum_to_kingdom)
        enrichment_df_phylum.loc[:, "domain"] = domain
        enrichment_df_phylum.to_csv(output[1], sep=",", mode="w", header=True, index=True)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="""Utility""")
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--N", type=int, default=1_000)
    parser.add_argument("--design", type=str, default="design.csv")
    parser.add_argument("--enrichment", type=str, default="")
    parser.add_argument("--output", type=str, default="test_A.csv,test_B.csv")
    parser.add_argument("--window_size", type=int, default=500)
    parser.add_argument("--domain", type=str, default="Bacteria", choices=["Eukaryota", "Archaea", "Viruses", "Bacteria"])
    args = parser.parse_args()
    N = args.N 
    alpha = args.alpha
    design = args.design 
    domain = args.domain
    enrichment_file = args.enrichment
    output = tuple(args.output.split(','))
    param = params(alpha=alpha, 
                   N=N, 
                   window_size=params.window_size)

    bootstrap_enrichment(design=design, 
                         enrichment_file=enrichment_file,
                         params=params,
                         domain=domain,
                         output=output)
