import pandas as pd
from os.path import basename
from typing import List, Dict, Tuple
from .utils import FastaParser
from .template import Processor


class Aggregate(Processor):

    blast_tabular_tsvs: List[str]
    min_percent_identity: float
    ref_fa: str
    query_fastas: List[str]

    count_df: pd.DataFrame
    percent_id_mean_df: pd.DataFrame
    percent_id_std_df: pd.DataFrame
    subject_id_to_taxon: Dict[str, str]
    sample_id_to_total_count: Dict[str, int]

    def main(
            self,
            blast_tabular_tsvs: List[str],
            min_percent_identity: float,
            ref_fa: str,
            query_fastas: List[str]) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:

        self.blast_tabular_tsvs = blast_tabular_tsvs
        self.min_percent_identity = min_percent_identity
        self.ref_fa = ref_fa
        self.query_fastas = query_fastas

        self.count_df = pd.DataFrame()
        self.percent_id_mean_df = pd.DataFrame()
        self.percent_id_std_df = pd.DataFrame()
        for tsv in self.blast_tabular_tsvs:
            self.process_one(tsv=tsv)
        self.count_df.fillna(0, inplace=True)

        self.set_subject_id_to_taxon()
        self.set_sample_id_to_total_count()
        self.calculate_unmapped_read_counts()
        self.label_rows_with_taxon()

        return self.count_df, self.percent_id_mean_df, self.percent_id_std_df

    def process_one(self, tsv: str):
        query_df = ReadBlastTsv(self.settings).main(tsv=tsv, min_percent_identity=self.min_percent_identity)

        sample_id = basename(tsv)[:-len('.tsv')]

        count_series = query_df.groupby('Subject ID').size()
        count_series.name = sample_id  # column name after join
        count_series.index.name = None
        self.count_df = self.count_df.join(count_series, how='outer')

        percent_id_mean_series = query_df.groupby('Subject ID')['Percent Identity'].mean()
        percent_id_mean_series.name = sample_id  # column name after join
        percent_id_mean_series.index.name = None
        self.percent_id_mean_df = self.percent_id_mean_df.join(percent_id_mean_series, how='outer')

        percent_id_std_series = query_df.groupby('Subject ID')['Percent Identity'].std()
        percent_id_std_series.name = sample_id  # column name after join
        percent_id_std_series.index.name = None
        self.percent_id_std_df = self.percent_id_std_df.join(percent_id_std_series, how='outer')

    def set_subject_id_to_taxon(self):
        self.subject_id_to_taxon = {}
        with FastaParser(self.ref_fa) as parser:
            for head, seq in parser:
                subject_id = head.split(' ')[0]
                self.subject_id_to_taxon[subject_id] = head

    def set_sample_id_to_total_count(self):
        self.sample_id_to_total_count = {}
        for fa in self.query_fastas:
            sample_id = basename(fa)[:-len('.fasta')]
            count = 0
            with FastaParser(fa) as parser:
                for _ in parser:
                    count += 1
            self.sample_id_to_total_count[sample_id] = count

    def calculate_unmapped_read_counts(self):
        for sample_id in self.count_df.columns:
            total = self.sample_id_to_total_count[sample_id]
            unmapped = total - self.count_df[sample_id].sum()
            self.count_df.at['Others', sample_id] = unmapped

    def label_rows_with_taxon(self):
        self.count_df.rename(index=self.subject_id_to_taxon, inplace=True)
        self.percent_id_mean_df.rename(index=self.subject_id_to_taxon, inplace=True)
        self.percent_id_std_df.rename(index=self.subject_id_to_taxon, inplace=True)


class ReadBlastTsv(Processor):

    tsv: str
    min_percent_identity: float

    df: pd.DataFrame
    query_df: pd.DataFrame

    def main(
            self,
            tsv: str,
            min_percent_identity: float) -> pd.DataFrame:

        self.tsv = tsv
        self.min_percent_identity = min_percent_identity

        self.df = pd.read_csv(tsv, sep='\t', header=None, names=[
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

        self.df = self.df[[
            'Query ID',
            'Subject ID',
            'Percent Identity',
        ]]

        self.df = self.df[self.df['Percent Identity'] >= self.min_percent_identity]

        self.df = self.df.sample(  # shuffle
            frac=1.
        ).sort_values(
            by='Percent Identity',  # high to low
            ascending=False
        )

        self.query_df = self.df.drop_duplicates(
            subset='Query ID',  # unique Query ID
            keep='first'
        ).reset_index(
            drop=True
        )

        return self.query_df
