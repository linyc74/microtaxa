import os
from shutil import rmtree
from .template import Settings
from .microtaxa import MicroTaxa
from .utils import get_temp_path


def entrypoint(
        sample_sheet: str,
        fq_dir: str,
        fq1_suffix: str,
        fq2_suffix: str,
        ref_fa: str,
        clip_r1_5_prime: int,
        clip_r2_5_prime: int,
        min_percent_identity: float,
        e_value: float,
        colormap: str,
        invert_colors: bool,
        publication_figure: bool,
        outdir: str,
        threads: int,
        debug: bool):

    settings = Settings(
        workdir=get_temp_path(prefix='./microtaxa_workdir_'),
        outdir=outdir,
        threads=threads,
        debug=debug,
        mock=False,
        for_publication=publication_figure)

    for d in [settings.workdir, settings.outdir]:
        os.makedirs(d, exist_ok=True)

    MicroTaxa(settings).main(
        ref_fa=ref_fa,
        sample_sheet=sample_sheet,
        fq_dir=fq_dir,
        fq1_suffix=fq1_suffix,
        fq2_suffix=fq2_suffix,
        min_percent_identity=min_percent_identity,
        e_value=e_value,
        clip_r1_5_prime=clip_r1_5_prime,
        clip_r2_5_prime=clip_r2_5_prime,
        colormap=colormap,
        invert_colors=invert_colors)

    if not debug:
        rmtree(settings.workdir)
