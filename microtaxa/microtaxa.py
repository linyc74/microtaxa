import pandas as pd
from os.path import basename
from typing import Optional, List, Tuple
from .template import Processor
from .merge import MergePairedEndReads
from .trimming import TrimGalorePairedEnd, TrimGaloreSingleEnd


class MicroTaxa(Processor):

    sample_sheet: str
    fq_dir: str
    fq1_suffix: str
    fq2_suffix: Optional[str]
    ref_fa: str
    min_percent_identity: float
    clip_r1_5_prime: int
    clip_r2_5_prime: int

    sample_ids: List[str]
    fastq_pairs: List[Tuple[str, Optional[str]]]
    trimmed_fastq_pairs: List[Tuple[str, Optional[str]]]
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
            min_percent_identity: float,
            clip_r1_5_prime: int,
            clip_r2_5_prime: int):

        self.sample_sheet = sample_sheet
        self.fq_dir = fq_dir
        self.fq1_suffix = fq1_suffix
        self.fq2_suffix = fq2_suffix
        self.ref_fa = ref_fa
        self.min_percent_identity = min_percent_identity
        self.clip_r1_5_prime = clip_r1_5_prime
        self.clip_r2_5_prime = clip_r2_5_prime

        self.read_sample_sheet()
        self.trim_galore()
        self.merge_paired_end_reads()
        self.convert_fastqs_to_fastas()
        self.run_glsearches()

    def read_sample_sheet(self):
        self.sample_ids = []
        self.fastq_pairs = []
        self.sample_ids = pd.read_csv(self.sample_sheet, index_col=0).index.tolist()
        for s in self.sample_ids:
            fq1 = f'{self.fq_dir}/{s}{self.fq1_suffix}'
            if self.fq2_suffix is None:
                fq2 = None
            else:
                fq2 = f'{self.fq_dir}/{s}{self.fq2_suffix}'
            self.fastq_pairs.append((fq1, fq2))

    def trim_galore(self):
        self.trimmed_fastq_pairs = []
        for fq1, fq2 in self.fastq_pairs:
            if fq2 is None:
                trimmed_fq = TrimGaloreSingleEnd(self.settings).main(
                    fq=fq1,
                    clip_5_prime=self.clip_r1_5_prime)
                self.trimmed_fastq_pairs.append((trimmed_fq, None))
            else:
                trimmed_fq1, trimmed_fq2 = TrimGalorePairedEnd(self.settings).main(
                    fq1=fq1,
                    fq2=fq2,
                    clip_r1_5_prime=self.clip_r1_5_prime,
                    clip_r2_5_prime=self.clip_r2_5_prime)
                self.trimmed_fastq_pairs.append((trimmed_fq1, trimmed_fq2))

    def merge_paired_end_reads(self):
        self.merged_fastqs = []
        for sample_id, fastq_pair in zip(self.sample_ids, self.trimmed_fastq_pairs):
            fq = MergePairedEndReads(self.settings).main(
                sample_id=sample_id,
                fastq_pair=fastq_pair
            )
            self.merged_fastqs.append(fq)

    def convert_fastqs_to_fastas(self):
        self.fastas = []
        for fq in self.merged_fastqs:
            fa = FastqToFasta(self.settings).main(fastq=fq)
            self.fastas.append(fa)

    def run_glsearches(self):
        self.glsearch_tsvs = []
        for fa in self.fastas:
            tsv = Glsearch(self.settings).main(
                query_fa=fa,
                library_fa=self.ref_fa,
                min_percent_identity=self.min_percent_identity)
            self.glsearch_tsvs.append(tsv)


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

    E_VALUE = 10e-30

    query_fa: str
    library_fa: str
    min_percent_identity: float

    output_tsv: str

    def main(
            self,
            query_fa: str,
            library_fa: str,
            min_percent_identity: float) -> str:

        self.query_fa = query_fa
        self.library_fa = library_fa
        self.min_percent_identity = min_percent_identity

        self.make_dstdir()

        fname = basename(self.query_fa)[:-len('.fasta')]
        self.output_tsv = f'{self.workdir}/{self.DSTDIR_NAME}/{fname}.tsv'
        args = [
            'glsearch36',
            '-3',  # forward strand only
            '-m 8',  # BLAST tabular output format
            '-n',  # DNA/RNA query
            f'-E {self.E_VALUE}',
            f'-T {self.threads}',
            self.query_fa,
            self.library_fa,
            f'1> {self.output_tsv}',
            f'2>> {self.outdir}/glsearch.log'
        ]
        self.call(self.CMD_LINEBREAK.join(args))

        return self.output_tsv

    def make_dstdir(self):
        self.call(f'mkdir -p {self.workdir}/{self.DSTDIR_NAME}')
