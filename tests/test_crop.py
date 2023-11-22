"""Test functionalities of gap filling."""
import cobra.io
import pathlib

path = pathlib.Path(__file__).parent


def load_model():
    return cobra.io.read_sbml_model(
        path.joinpath('example_model.xml').__str__()
    )


class TestExampleModel:
    """ Tests that the example model behaves according to the conditions that we want.

    Basically this is step 1 in the process. We have a model that fails in a designed way. That way
    as we develop crop, we can have a "perfect" model that we know behaves the way we would like. We found one
    example already where the model grew but we expected it not to, thus we added a test to make sure.

    """
    # def __init__(self):
    model = load_model()

    def test_no_grow_e(self):
        """ E  nutrient condition is no growth """

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
        """ B + E  nutrient condition is no growth

        This is the function that we expected to fail, but turns out the model can convert to double growth.

        Returns:

        """

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
