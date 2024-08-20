from typing import Optional, Tuple
from .template import Processor


class MergePairedEndReads(Processor):

    DSTDIR_NAME = 'merged-fastq'
    MIN_OVERLAP = 10

    sample_id: str
    fastq_pair: Tuple[str, Optional[str]]

    dstdir: str
    output_fastq: str

    def main(
            self,
            sample_id: str,
            fastq_pair: Tuple[str, Optional[str]]) -> str:

        self.sample_id = sample_id
        self.fastq_pair = fastq_pair

        self.make_dstdir()
        self.set_output_fastq()

        fq2 = self.fastq_pair[1]
        if fq2 is None:
            self.copy_fq1()
        else:
            self.merge_fq1_fq2()

        return self.output_fastq

    def make_dstdir(self):
        self.dstdir = f'{self.workdir}/{self.DSTDIR_NAME}'
        self.call(f'mkdir -p {self.dstdir}')

    def set_output_fastq(self):
        self.output_fastq = f'{self.dstdir}/{self.sample_id}.fastq.gz'

    def copy_fq1(self):
        fq1 = self.fastq_pair[0]
        if fq1.endswith('.gz'):
            self.call(f'cp {fq1} {self.output_fastq}')
        else:
            self.call(f'gzip -c {fq1} > {self.output_fastq}')

    def merge_fq1_fq2(self):
        temp_dir = f'{self.workdir}/pear-temp'
        self.call(f'mkdir -p {temp_dir}')

        output_prefix = f'{temp_dir}/{self.sample_id}'
        log = f'{self.outdir}/pear.log'
        args = [
            'pear',
            f'--forward-fastq {self.fastq_pair[0]}',
            f'--reverse-fastq {self.fastq_pair[1]}',
            f'--output {output_prefix}',
            f'--min-overlap {self.MIN_OVERLAP}',
            f'--threads {self.threads}',
            f'1>> "{log}"',
            f'2>> "{log}"'
        ]
        self.call(self.CMD_LINEBREAK.join(args))

        self.call(f'gzip -c {output_prefix}.assembled.fastq > {self.output_fastq}')
        self.call(f'rm -r {temp_dir}')
