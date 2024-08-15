from microtaxa.aggregate import Aggregate
from .setup import TestCase


class TestMicroTaxa(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    # def tearDown(self):
    #     self.tear_down()

    def test_main(self):
        Aggregate(self.settings).main(
            blast_tabular_tsvs=[
                f'{self.indir}/glsearch/EPI-001.tsv',
                f'{self.indir}/glsearch/EPI-002.tsv',
                f'{self.indir}/glsearch/EPI-003.tsv',
                f'{self.indir}/glsearch/TP-184.tsv',
                f'{self.indir}/glsearch/TP-190.tsv',
                f'{self.indir}/glsearch/TP-202.tsv',
            ],
            min_percent_identity=90.0,
        )
