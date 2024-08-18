import pandas as pd
from microtaxa.heatmap import PlotHeatmaps, PlotOneHeatmap
from .setup import TestCase


class TestPlotHeatmaps(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        self.settings.for_publication = True
        PlotHeatmaps(self.settings).main(
            count_df=pd.read_csv(f'{self.indir}/count-table.csv', index_col=0),
            percent_id_mean_df=pd.read_csv(f'{self.indir}/percent-identity-mean.csv', index_col=0),
            percent_id_std_df=pd.read_csv(f'{self.indir}/percent-identity-std.csv', index_col=0),
            sample_sheet=f'{self.indir}/sample-sheet.csv',
        )


class TestPlotOneHeatmap(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        self.settings.for_publication = True
        PlotOneHeatmap(self.settings).main(
            df=pd.read_csv(f'{self.indir}/count-table.csv', index_col=0),
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            log_pseudocount=True,
            normalize_by_sample_reads=False,
            output_fname='log-pseudocount',
        )
