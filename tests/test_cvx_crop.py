"""Test functionalities of gap filling."""
import cobra.io
import pathlib

path = pathlib.Path(__file__).parent


def load_model():
    return cobra.io.read_sbml_model(
        path.joinpath('example_model.xml').__str__()
    )


class TestClass:
    # def __init__(self):
    model = load_model()

    def test_no_grow_e(self):
        # E  nutrient condition is no growth

        media = {'EX_E_e': 100}
        with self.model as m:
            m.medium = media
            obj_func = m.slim_optimize()
            assert obj_func == 0.0

    def test_no_grow_b(self):
        """B  nutrient condition is no growth"""

        media = {'EX_B_e': 100}
        with self.model as m:
            m.medium = media
            obj_func = m.slim_optimize()
            assert obj_func == 0.0

    def test_no_grow_b_and_e(self):
        # B + E  nutrient condition is no growth

        media = {'EX_E_e': 100, 'EX_B_e': 100}
        with self.model as m:
            m.medium = media
            obj_func = m.slim_optimize()
            assert obj_func == 200.0

    def test_grow_a_and_e(self):
        """A + E  nutrient condition is growth"""
        media = {'EX_E_e': 100, 'EX_A_e': 100}
        with self.model as m:
            m.medium = media
            obj_func = m.slim_optimize()
            assert obj_func > 0.0

    def test_grow_a_and_b(self):
        """A + B  nutrient condition is growth"""

        media = {'EX_B_e': 100, 'EX_A_e': 100}
        with self.model as m:
            m.medium = media
            obj_func = m.slim_optimize()
            assert obj_func > 0.0

    def test_grow_a(self):
        """ A alone is growth """

        media = {'EX_A_e': 100}
        with self.model as m:
            m.medium = media
            obj_func = m.slim_optimize()
            assert obj_func == 100.0

    def test_no_grow_4(self):
        """A + v3 knockout is no-growth"""
        media = {'EX_A_e': 100}
        with self.model as m:
            m.remove_reactions([m.reactions.get_by_id('R_A_to_C')])
            m.medium = media
            obj_func = m.slim_optimize()

            assert obj_func == 0.0

    def test_no_grow_5(self):
        """A + E + v3 knockout is growth"""

        media = {
            'EX_A_e': 100,
            'EX_E_e': 100
        }
        with self.model as m:
            m.remove_reactions([m.reactions.get_by_id('R_A_to_C')])
            m.medium = media
            obj_func = m.slim_optimize()
            assert obj_func > 0.0
