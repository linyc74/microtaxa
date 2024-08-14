from microtaxa.trimming import TrimGaloreSingleEnd, TrimGalorePairedEnd
from .setup import TestCase


class TestTrimGaloreSingleEnd(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        actual = TrimGaloreSingleEnd(self.settings).main(
            fq=f'{self.indir}/H-1_S24_R1.fastq.gz',
            clip_5_prime=1
        )
        expected = f'{self.workdir}/trimmed-fastq/H-1_S24_R1_trimmed.fq.gz'
        self.assertFileExists(expected, actual)


class TestTrimGalorePairedEnd(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        fq1, fq2 = TrimGalorePairedEnd(self.settings).main(
            fq1=f'{self.indir}/H-1_S24_R1.fastq.gz',
            fq2=f'{self.indir}/H-1_S24_R2.fastq.gz',
            clip_r1_5_prime=1,
            clip_r2_5_prime=1
        )
        self.assertFileExists(f'{self.workdir}/trimmed-fastq/H-1_S24_R1_val_1.fq.gz', fq1)
        self.assertFileExists(f'{self.workdir}/trimmed-fastq/H-1_S24_R2_val_2.fq.gz', fq2)
