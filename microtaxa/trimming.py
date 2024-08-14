import os
from os.path import basename
from typing import Tuple
from .template import Processor


class TrimGalore(Processor):

    DSTDIR_NAME = 'trimmed-fastq'
    QUALITY = 20
    LENGTH = 20
    MAX_N = 0
    CUTADAPT_TOTAL_CORES = 2
    # According to the help message of trim_galore, 2 cores for cutadapt -> actually up to 9 cores

    dstdir: str

    def make_dstdir(self):
        self.dstdir = f'{self.workdir}/{self.DSTDIR_NAME}'
        self.call(f'mkdir -p "{self.dstdir}"')

    def move_fastqc_report(self):
        dstdir = f'{self.outdir}/fastqc'
        os.makedirs(dstdir, exist_ok=True)
        for suffix in [
            'fastqc.html',
            'fastqc.zip',
            'trimming_report.txt'
        ]:
            self.call(f'mv {self.dstdir}/*{suffix} "{dstdir}/"')


class TrimGalorePairedEnd(TrimGalore):

    fq1: str
    fq2: str
    clip_r1_5_prime: int
    clip_r2_5_prime: int

    out_fq1: str
    out_fq2: str

    def main(
            self,
            fq1: str,
            fq2: str,
            clip_r1_5_prime: int,
            clip_r2_5_prime: int) -> Tuple[str, str]:

        self.fq1 = fq1
        self.fq2 = fq2
        self.clip_r1_5_prime = clip_r1_5_prime
        self.clip_r2_5_prime = clip_r2_5_prime

        self.make_dstdir()
        self.execute()
        self.move_fastqc_report()
        self.set_out_fq1()
        self.set_out_fq2()

        return self.out_fq1, self.out_fq2

    def execute(self):
        args = [
            'trim_galore',
            '--paired',
            f'--quality {self.QUALITY}',
            '--phred33',
            f'--cores {self.CUTADAPT_TOTAL_CORES}',
            f'--fastqc_args "--threads {self.threads}"',
            '--illumina',
            f'--length {self.LENGTH}',
            f'--max_n {self.MAX_N}',
            '--trim-n',
            '--gzip',
            f'--output_dir {self.dstdir}'
        ]

        if self.clip_r1_5_prime > 0:
            args.append(f'--clip_R1 {self.clip_r1_5_prime}')

        if self.clip_r2_5_prime > 0:
            args.append(f'--clip_R2 {self.clip_r2_5_prime}')

        log = f'{self.outdir}/trim_galore.log'
        args += [
            self.fq1,
            self.fq2,
            f'1>> "{log}"',
            f'2>> "{log}"'
        ]

        self.call(self.CMD_LINEBREAK.join(args))

    def set_out_fq1(self):
        f = basename(self.fq1)
        f = self.__strip_file_extension(f)
        self.out_fq1 = f'{self.dstdir}/{f}_val_1.fq.gz'

    def set_out_fq2(self):
        f = basename(self.fq2)
        f = self.__strip_file_extension(f)
        self.out_fq2 = f'{self.dstdir}/{f}_val_2.fq.gz'

    def __strip_file_extension(self, f):
        for suffix in [
            '.fq',
            '.fq.gz',
            '.fastq',
            '.fastq.gz',
        ]:
            if f.endswith(suffix):
                f = f[:-len(suffix)]  # strip suffix
        return f


class TrimGaloreSingleEnd(TrimGalore):

    fq: str
    clip_5_prime: int

    out_fq: str

    def main(self, fq: str, clip_5_prime: int) -> str:
        self.fq = fq
        self.clip_5_prime = clip_5_prime

        self.make_dstdir()
        self.execute()
        self.move_fastqc_report()
        self.set_out_fq()

        return self.out_fq

    def execute(self):
        args = [
            'trim_galore',
            f'--quality {self.QUALITY}',
            '--phred33',
            f'--cores {self.CUTADAPT_TOTAL_CORES}',
            f'--fastqc_args "--threads {self.threads}"',
            '--illumina',
            f'--length {self.LENGTH}',
            f'--max_n {self.MAX_N}',
            '--trim-n',
            '--gzip',
            f'--output_dir {self.dstdir}',
        ]

        if self.clip_5_prime > 0:
            args.append(f'--clip_R1 {self.clip_5_prime}')

        log = f'{self.outdir}/trim_galore.log'
        args += [
            self.fq,
            f'1>> "{log}"',
            f'2>> "{log}"'
        ]

        self.call(self.CMD_LINEBREAK.join(args))

    def set_out_fq(self):
        f = basename(self.fq)
        f = self.__strip_file_extension(f)
        self.out_fq = f'{self.dstdir}/{f}_trimmed.fq.gz'

    def __strip_file_extension(self, f):
        for suffix in [
            '.fq',
            '.fq.gz',
            '.fastq',
            '.fastq.gz',
        ]:
            if f.endswith(suffix):
                f = f[:-len(suffix)]  # strip suffix
        return f
