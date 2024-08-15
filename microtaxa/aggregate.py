import pandas as pd
from typing import List
from .template import Processor


class Aggregate(Processor):

    blast_tabular_tsvs: List[str]
    min_percent_identity: float

    def main(
            self,
            blast_tabular_tsvs: List[str],
            min_percent_identity: float):

        self.blast_tabular_tsvs = blast_tabular_tsvs
        self.min_percent_identity = min_percent_identity

        for tsv in self.blast_tabular_tsvs:
            df = ProcessBlastTabularTsv(self.settings).main(tsv=tsv, min_percent_identity=self.min_percent_identity)
            print(df)


class ProcessBlastTabularTsv(Processor):

    tsv: str
    min_percent_identity: float

    df: pd.DataFrame
    count_df: pd.DataFrame
    percent_id_df: pd.DataFrame
    percent_id_std_df: pd.DataFrame

    def main(
            self,
            tsv: str,
            min_percent_identity: float) -> pd.DataFrame:

        self.tsv = tsv
        self.min_percent_identity = min_percent_identity

        self.df = read_blast(self.tsv)

        self.df = self.df[self.df['Percent Identity'] >= self.min_percent_identity]

        self.df = self.df.sample(  # shuffle
            frac=1.
        ).reset_index(
            drop=True
        ).sort_values(
            by='Percent Identity',  # high to low
            ascending=False
        ).drop_duplicates(
            subset='Query ID',
            keep='first'
        )

        return self.df


def read_blast(tsv: str) -> pd.DataFrame:
    return pd.read_csv(tsv, sep='\t', header=None, names=[
        'Query ID',
        'Subject ID',
        'Percent Identity',
        'Alignment Length',
        'Number of Mismatches',
        'Number of Gap Openings',
        'Query Start',
        'Query End',
        'Subject Start',
        'Subject End',
        'E-value',
        'Bit Score'
    ])
