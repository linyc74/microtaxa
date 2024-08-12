from microtaxa.microtaxa import MicroTaxa
from .setup import TestCase


class TestMicroTaxa(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    # def tearDown(self):
    #     self.tear_down()

    def test_main(self):
        MicroTaxa(self.settings).main(
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            fq_dir=f'{self.indir}/fq-dir',
            fq1_suffix='_R1.fastq.gz',
            fq2_suffix=None,
            ref_fa=f'{self.indir}/16S.fa',
            min_percent_id=97.0,
            clip_r1_5_prime=0,
            clip_r2_5_prime=0
        )
