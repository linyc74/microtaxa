from microtaxa.microtaxa import MicroTaxa
from .setup import TestCase


class TestMicroTaxa(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        MicroTaxa(self.settings).main()
