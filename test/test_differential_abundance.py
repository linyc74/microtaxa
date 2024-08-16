from microtaxa.differential_abundance import DifferentialAbundance
from .setup import TestCase


class TestDifferentialAbundance(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    # def tearDown(self):
    #     self.tear_down()

    def test_normal_format(self):
        DifferentialAbundance(self.settings).main(
            taxon_csv=f'{self.indir}/count-table-normal.csv',
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0)],
        )

    def __test_silva_format(self):
        DifferentialAbundance(self.settings).main(
            taxon_csv=f'{self.indir}/count-table-silva.csv',
            sample_sheet=f'{self.indir}/sample-sheet.csv',
            colors=[(0.2, 0.5, 0.7, 1.0), (0.9, 0.1, 0.1, 1.0)],
        )
