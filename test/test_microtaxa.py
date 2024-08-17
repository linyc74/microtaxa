from microtaxa.microtaxa import MicroTaxa
from .setup import TestCase


class TestMicroTaxa(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        MicroTaxa(self.settings).main(
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            fq_dir=f'{self.indir}/fq-dir',
            fq1_suffix='_R1.fastq.gz',
            fq2_suffix='_R2.fastq.gz',
            ref_fa=f'{self.indir}/reference.fasta',
            min_percent_identity=97.0,
            clip_r1_5_prime=0,
            clip_r2_5_prime=0,
            colormap='viridis',
            invert_colors=False
        )
