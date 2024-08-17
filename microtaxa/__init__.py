import os
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
        colormap: str,
        invert_colors: bool,
        outdir: str,
        threads: int,
        debug: bool):

    settings = Settings(
        workdir=get_temp_path(prefix='./microtaxa_workdir_'),
        outdir=outdir,
        threads=threads,
        debug=debug,
        mock=False)

    for d in [settings.workdir, settings.outdir]:
        os.makedirs(d, exist_ok=True)

    MicroTaxa(settings).main(
        ref_fa=ref_fa,
        sample_sheet=sample_sheet,
        fq_dir=fq_dir,
        fq1_suffix=fq1_suffix,
        fq2_suffix=fq2_suffix,
        min_percent_identity=min_percent_identity,
        clip_r1_5_prime=clip_r1_5_prime,
        clip_r2_5_prime=clip_r2_5_prime,
        colormap=colormap,
        invert_colors=invert_colors)
