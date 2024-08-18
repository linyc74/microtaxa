import pandas as pd

from microtaxa.aggregate import Aggregate
from .setup import TestCase


class TestMicroTaxa(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        count_df, percent_id_mean_df, percent_id_std_df = Aggregate(self.settings).main(
            blast_tabular_tsvs=[
                f'{self.indir}/glsearch/EPI-001.tsv',
                f'{self.indir}/glsearch/EPI-002.tsv',
                f'{self.indir}/glsearch/EPI-003.tsv',
                f'{self.indir}/glsearch/TP-184.tsv',
                f'{self.indir}/glsearch/TP-190.tsv',
                f'{self.indir}/glsearch/TP-202.tsv',
            ],
            min_percent_identity=90.0,
            ref_fa=f'{self.indir}/reference.fasta',
            query_fastas=[
                f'{self.indir}/fasta/EPI-001.fasta',
                f'{self.indir}/fasta/EPI-002.fasta',
                f'{self.indir}/fasta/EPI-003.fasta',
                f'{self.indir}/fasta/TP-184.fasta',
                f'{self.indir}/fasta/TP-190.fasta',
                f'{self.indir}/fasta/TP-202.fasta',
            ]
        )

        self.assertDataFrameEqual(
            count_df,
            pd.read_csv(f'{self.indir}/count-table.csv', index_col=0)
        )
        self.assertDataFrameEqual(
            percent_id_mean_df,
            pd.read_csv(f'{self.indir}/percent-identity-mean.csv', index_col=0)
        )
        self.assertDataFrameEqual(
            percent_id_std_df,
            pd.read_csv(f'{self.indir}/percent-identity-std.csv', index_col=0)
        )
