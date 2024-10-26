import pandas as pd

def find_top_10(df: pd.DataFrame) -> pd.Series:
    return df['sequenceOfArm'].value_counts().sort_values(ascending=False).head(10)

def find_spacer_distribution(df: pd.DataFrame) -> pd.Series:
    return 1e2 * df['spacerLength'].value_counts(normalize=True).sort_values(ascending=False)

def find_gc_content(df: pd.DataFrame) -> pd.Series:
    df['gc_content'] = 1e2 * (df['arm_g'] + df['arm_c']) / df['sequenceLength']

def subcompartment_coverage(df: pd.DataFrame) -> pd.DataFrame:
    pass


