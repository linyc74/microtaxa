import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from typing import Tuple
from .template import Processor
from .normalization import CountNormalization
from .grouping import TagGroupNamesOnSampleColumns


DSTDIR_NAME = 'heatmap'


class PlotHeatmaps(Processor):

    def main(
            self,
            count_df: pd.DataFrame,
            percent_id_mean_df: pd.DataFrame,
            percent_id_std_df: pd.DataFrame,
            sample_sheet: str):

        PlotOneHeatmap(self.settings).main(
            df=count_df,
            sample_sheet=sample_sheet,
            log_pseudocount=True,
            normalize_by_sample_reads=False,
            output_fname='log-pseudocount'
        )

        PlotOneHeatmap(self.settings).main(
            df=percent_id_mean_df.fillna(0),
            sample_sheet=sample_sheet,
            log_pseudocount=False,
            normalize_by_sample_reads=False,
            output_fname='percent-identity-mean'
        )

        PlotOneHeatmap(self.settings).main(
            df=percent_id_std_df.fillna(0),
            sample_sheet=sample_sheet,
            log_pseudocount=False,
            normalize_by_sample_reads=False,
            output_fname='percent-identity-std'
        )


class PlotOneHeatmap(Processor):

    df: pd.DataFrame
    sample_sheet: str
    log_pseudocount: bool
    normalize_by_sample_reads: bool
    output_fname: str

    dstdir: str

    def main(
            self,
            df: pd.DataFrame,
            sample_sheet: str,
            log_pseudocount: bool,
            normalize_by_sample_reads: bool,
            output_fname: str):

        self.df = df
        self.sample_sheet = sample_sheet
        self.log_pseudocount = log_pseudocount
        self.normalize_by_sample_reads = normalize_by_sample_reads
        self.output_fname = output_fname

        self.dstdir = f'{self.outdir}/{DSTDIR_NAME}'
        os.makedirs(self.dstdir, exist_ok=True)

        self.df = CountNormalization(self.settings).main(
            df=self.df,
            log_pseudocount=self.log_pseudocount,
            by_sample_reads=self.normalize_by_sample_reads)

        Clustermap(self.settings).main(
            data=self.df,
            sample_sheet=self.sample_sheet,
            output_prefix=f'{self.dstdir}/{self.output_fname}')


class Clustermap(Processor):

    CLUSTER_COLUMNS = True
    COLORMAP = 'PuBu'
    Y_LABEL_CHAR_WIDTH = 0.14 / 2.54
    X_LABEL_CHAR_WIDTH = 0.14 / 2.54
    CELL_WIDTH = 0.4 / 2.54
    CELL_HEIGHT = 0.4 / 2.54
    DENDROGRAM_SIZE = 1.0 / 2.54
    COLORBAR_WIDTH = 0.01
    COLORBAR_HORIZONTAL_POSITION = 1.
    FONTSIZE = 7
    LINE_WIDTH = 0.5
    DPI = 600

    data: pd.DataFrame
    sample_sheet: str
    output_prefix: str

    x_label_padding: float
    y_label_padding: float
    figsize: Tuple[float, float]
    grid: sns.matrix.ClusterGrid

    def main(
            self,
            data: pd.DataFrame,
            sample_sheet: str,
            output_prefix: str):

        self.data = data.copy()
        self.sample_sheet = sample_sheet
        self.output_prefix = output_prefix

        self.tag_group_names_on_sample_columns()
        self.shorten_taxon_names_for_publication()
        self.set_figsize()
        self.clustermap()
        self.config_clustermap()
        self.save_fig()
        self.save_csv()

    def tag_group_names_on_sample_columns(self):
        self.data = TagGroupNamesOnSampleColumns(self.settings).main(
            df=self.data,
            sample_sheet=self.sample_sheet)

    def shorten_taxon_names_for_publication(self):
        if not self.settings.for_publication:
            return

        def shorten_silva(s: str) -> str:
            """
            SILVA taxonomy format:

            AY188352.1.1546 Bacteria;Bacillota;Bacilli;Lactobacillales;Streptococcaceae;Streptococcus;Streptococcus salivarius

            Shortened format:

            AY188352.1.1546 Streptococcus salivarius
            """
            prefix = s.split(' ')[0]
            suffix = s.split(';')[-1]
            return f'{prefix} {suffix}'

        rename = {}
        for idx in self.data.index:
            if ';' in idx:
                rename[idx] = shorten_silva(idx)

        self.data.rename(index=rename, inplace=True)

    def set_figsize(self):
        self.__set_x_y_label_padding()
        w = (len(self.data.columns) * self.CELL_WIDTH) + self.y_label_padding
        h = (len(self.data.index) * self.CELL_HEIGHT) + self.x_label_padding
        self.figsize = (w, h)

    def __set_x_y_label_padding(self):
        max_x_label_length = pd.Series(self.data.columns).apply(len).max()
        self.x_label_padding = max_x_label_length * self.X_LABEL_CHAR_WIDTH + self.DENDROGRAM_SIZE

        max_y_label_length = pd.Series(self.data.index).apply(len).max()
        self.y_label_padding = max_y_label_length * self.Y_LABEL_CHAR_WIDTH + self.DENDROGRAM_SIZE

    def clustermap(self):
        plt.rcParams['font.size'] = self.FONTSIZE
        plt.rcParams['axes.linewidth'] = self.LINE_WIDTH
        w, h = self.figsize
        dendrogram_ratio = (self.DENDROGRAM_SIZE / w, self.DENDROGRAM_SIZE / h)
        self.grid = sns.clustermap(
            data=self.data,
            cmap=self.COLORMAP,
            figsize=self.figsize,
            xticklabels=True,  # include every x label
            yticklabels=True,  # include every y label
            col_cluster=self.CLUSTER_COLUMNS,
            dendrogram_ratio=dendrogram_ratio,
            linewidth=0.25)
        self.__set_plotted_data()

    def __set_plotted_data(self):
        self.data = self.grid.__dict__['data2d']

    def config_clustermap(self):
        self.__set_x_y_axes()
        self.__set_colorbar()

    def __set_x_y_axes(self):
        heatmap = self.grid.ax_heatmap

        n_rows = len(self.data)
        heatmap.set_ylim(n_rows, 0)  # first and last row not be chopped half

        heatmap.tick_params(
            axis='x',
            bottom=True,
            top=False,
            labelbottom=True,
            labeltop=False,
            labelrotation=90
        )
        heatmap.tick_params(
            axis='y',
            left=False,
            right=True,
            labelleft=False,
            labelright=True,
            labelrotation=0,
        )

    def __set_colorbar(self):
        colorbar = self.grid.cax
        p = colorbar.get_position()
        colorbar.set_position([
            self.COLORBAR_HORIZONTAL_POSITION,  # left
            p.y0,  # bottom
            self.COLORBAR_WIDTH,  # width
            p.height  # height
        ])

    def save_fig(self):
        # must use grid.savefig(), but not plt.savefig()
        # plt.savefig() crops out the colorbar

        plt.tight_layout()
        dpi = self.__downsize_dpi_if_too_large()
        for ext in ['pdf', 'png']:
            self.grid.savefig(f'{self.output_prefix}.{ext}', dpi=dpi)
        plt.close()

    def __downsize_dpi_if_too_large(self) -> int:
        longer_side = max(self.figsize)
        dpi = self.DPI
        while longer_side * dpi >= 2**15:  # 2^16 is the limit of matplotlib, but 2^15 is safer after some tests
            dpi = int(dpi/2)  # downsize
        return dpi

    def save_csv(self):
        self.data.to_csv(f'{self.output_prefix}.tsv', sep='\t', index=True)
