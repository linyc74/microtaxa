from microtaxa.merge import MergePairedEndReads
from .setup import TestCase


class TestMergePairedEndReads(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_single_end(self):
        copied_fq1 = MergePairedEndReads(self.settings).main(
            sample_id='S01',
            fastq_pair=(f'{self.indir}/R1.fastq.gz', None)
        )
        self.assertFileExists(f'{self.workdir}/merged-fastq/S01.fastq.gz', copied_fq1)

    def test_paired_end(self):
        merged_fq = MergePairedEndReads(self.settings).main(
            sample_id='S01',
            fastq_pair=(f'{self.indir}/R1.fastq.gz', f'{self.indir}/R2.fastq.gz')
        )
        self.assertFileExists(f'{self.workdir}/merged-fastq/S01.fastq.gz', merged_fq)
