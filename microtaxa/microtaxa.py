import pandas as pd
from os.path import basename
from typing import Optional, List, Tuple
from .template import Processor


class MicroTaxa(Processor):

    sample_sheet: str
    fq_dir: str
    fq1_suffix: str
    fq2_suffix: Optional[str]
    ref_fa: str
    min_percent_id: float
    clip_r1_5_prime: int
    clip_r2_5_prime: int

    fastq_pairs: List[Tuple[str, Optional[str]]]
    merged_fastqs: List[str]
    fastas: List[str]
    glsearch_tsvs: List[str]

    def main(
            self,
            sample_sheet: str,
            fq_dir: str,
            fq1_suffix: str,
            fq2_suffix: Optional[str],
            ref_fa: str,
            min_percent_id: float,
            clip_r1_5_prime: int,
            clip_r2_5_prime: int):

        self.sample_sheet = sample_sheet
        self.fq_dir = fq_dir
        self.fq1_suffix = fq1_suffix
        self.fq2_suffix = fq2_suffix
        self.ref_fa = ref_fa
        self.min_percent_id = min_percent_id
        self.clip_r1_5_prime = clip_r1_5_prime
        self.clip_r2_5_prime = clip_r2_5_prime

        self.set_fastq_pairs()
        self.merge_paired_end_reads()
        self.convert_fastqs_to_fastas()
        self.run_glsearches()

    def set_fastq_pairs(self):
        self.fastq_pairs = []
        sample_ids = pd.read_csv(self.sample_sheet, index_col=0).index.tolist()
        for s in sample_ids:
            fq1 = f'{self.fq_dir}/{s}{self.fq1_suffix}'
            if self.fq2_suffix is None:
                fq2 = None
            else:
                fq2 = f'{self.fq_dir}/{s}{self.fq2_suffix}'
            self.fastq_pairs.append((fq1, fq2))

    def merge_paired_end_reads(self):
        if self.fq2_suffix is not None:
            self.merged_fastqs = MergePairedEndReads(self.settings).main(self.fastq_pairs)
        else:
            self.merged_fastqs = [fq for fq, _ in self.fastq_pairs]

    def convert_fastqs_to_fastas(self):
        self.fastas = []
        for fq in self.merged_fastqs:
            fa = FastqToFasta(self.settings).main(fq)
            self.fastas.append(fa)

    def run_glsearches(self):
        self.glsearch_tsvs = []
        for fa in self.fastas:
            tsv = Glsearch(self.settings).main(
                query_fa=fa,
                library_fa=self.ref_fa,
                min_percent_id=self.min_percent_id)
            self.glsearch_tsvs.append(tsv)


class MergePairedEndReads(Processor):

    fastq_pairs: List[Tuple[str, str]]

    merged_fastqs: List[str]

    def main(self, fastq_pairs: List[Tuple[str, str]]) -> List[str]:

        self.fastq_pairs = fastq_pairs

        self.merged_fastqs = []

        return self.merged_fastqs


class FastqToFasta(Processor):

    DSTDIR_NAME = 'fasta'

    fastq: str

    fasta: str

    def main(self, fastq: str) -> str:
        self.fastq = fastq

        self.make_dstdir()
        self.set_fasta_path()
        self.run_seqtk()

        return self.fasta

    def make_dstdir(self):
        self.call(f'mkdir -p {self.workdir}/{self.DSTDIR_NAME}')

    def set_fasta_path(self):
        fname = basename(self.fastq)
        if fname.endswith('.gz'):
            fname = fname[:-3]
        if fname.endswith('.fq'):
            fname = fname[:-3]
        if fname.endswith('.fastq'):
            fname = fname[:-6]
        self.fasta = f'{self.workdir}/{self.DSTDIR_NAME}/{fname}.fasta'

    def run_seqtk(self):
        lines = [
            'seqtk seq',
            f'-a {self.fastq}',
            f'> {self.fasta}'
        ]
        self.call(self.CMD_LINEBREAK.join(lines))


class Glsearch(Processor):

    DSTDIR_NAME = 'glsearch'

    E_VALUE = 10e-6

    query_fa: str
    library_fa: str
    min_percent_id: float

    output_tsv: str

    def main(
            self,
            query_fa: str,
            library_fa: str,
            min_percent_id: float) -> str:

        self.query_fa = query_fa
        self.library_fa = library_fa
        self.min_percent_id = min_percent_id

        self.make_dstdir()

        fname = basename(self.query_fa)[:-len('.fasta')]
        self.output_tsv = f'{self.workdir}/{self.DSTDIR_NAME}/{fname}.tsv'
        args = [
            'glsearch36',
            '-3',  # forward strand only
            f'-E {self.E_VALUE}',
            '-m 8',  # BLAST tabular output format
            '-n',  # DNA/RNA query
            f'-O {self.output_tsv}',
            self.query_fa,
            self.library_fa,
            # f'1>> {self.outdir}/glsearch.log',
            # f'2>> {self.outdir}/glsearch.log'
        ]
        self.call(self.CMD_LINEBREAK.join(args))

        return self.output_tsv

    def make_dstdir(self):
        self.call(f'mkdir -p {self.workdir}/{self.DSTDIR_NAME}')


