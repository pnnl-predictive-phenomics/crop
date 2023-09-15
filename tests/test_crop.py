"""Test functionalities of gap filling."""
import cobra.io
import pathlib

path = pathlib.Path(__file__).parent


def test_no_grow_e():
    # B + E  nutrient condition is no growth
    model = cobra.io.read_sbml_model(
        path.joinpath('example_model.xml').__str__()
    )
    media = {'EX_E_e': 100}
    model.medium = media
    obj_func = model.slim_optimize()
    assert obj_func == 0.0


def test_no_grow_a():
    # A alone is growth
    model = cobra.io.read_sbml_model(
        path.joinpath('example_model.xml').__str__()
    )
    media = {'EX_A_e': 100}
    model.medium = media
    obj_func = model.slim_optimize()
    assert obj_func == 0.0


def test_no_grow_4():
    # A + v3 knockout is no-growth
    model = cobra.io.read_sbml_model(
        path.joinpath('example_model.xml').__str__()
    )
    media = {'EX_A_e': 100}
    model.remove_reactions(model.reactions.get_by_id('R_A_to_C'))
    model.medium = media
    obj_func = model.slim_optimize()
    assert obj_func == 0.0


def test_no_grow_5():
    #A + E + v3 knockout is growth
    model = cobra.io.read_sbml_model(
        path.joinpath('example_model.xml').__str__()
    )
    media = {
        'EX_A_e': 100,
        'EX_E_e': 100
    }
    model.remove_reactions(model.reactions.get_by_id('R_A_to_C'))
    model.medium = media
    obj_func = model.slim_optimize()
    assert obj_func > 0.0
