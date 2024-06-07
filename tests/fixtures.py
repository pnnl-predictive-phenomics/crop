import cobra as cb
import optlang as op
import pandas as pd
import pytest
from optlang import Variable

import crop
from crop import (
    PhenotypeObservation,
    generate_decision_variables,
    get_steady_state_dual_constraint,
    model_from_stoich_matrix,
)


# pytest fixtures
@pytest.fixture
def expected_model():
    return cb.io.load_json_model("./tests/ABC_toy_model_3.json")


@pytest.fixture
def actual_model():
    S_index = ["A_int", "B_int", "C_int"]
    S_dict = {
        "A_SRC->A_int": [1, 0, 0],
        "A_int->B_int": [-1, 1, 0],
        "B_int->C_int": [0, -1, 1],
        "C_int->C_SNK": [0, 0, -1],
        "A_int->C_int": [-1, 0, 1],
        "B_SRC->B_int": [0, 1, 0],
    }
    S_table = pd.DataFrame(S_dict, index=S_index)
    mets, rxns = S_table.index, S_table.columns
    n_mets = len(mets)
    n_rxns = len(rxns)

    growth_objective = {rxn_id: 1 if rxn_id == "C_int->C_SNK" else 0 for rxn_id in rxns}
    upper_flux_bounds = {rxn_id: 0 for rxn_id in rxns}
    lower_flux_bounds = {rxn_id: 0 for rxn_id in rxns}

    model = model_from_stoich_matrix(
        S_table,
        "ABC_toy_model_3",
        obj=growth_objective,
        lower_flux_bounds=lower_flux_bounds,
        upper_flux_bounds=upper_flux_bounds,
    )
    return model


@pytest.fixture
def phenotype_observations():

    # flux bounds
    min_growth_flux = 10
    max_nogrowth_flux = 5

    # media conditions
    A_not_B_medium = {"A_SRC->A_int": min_growth_flux, "B_SRC->B_int": 0.0}

    B_not_A_medium = {"A_SRC->A_int": 0.0, "B_SRC->B_int": min_growth_flux}

    # knockout reactions
    reaction_knockouts = {
        "A_int->C_int": True,
    }

    phenotype_observations = {
        "A_not_B": PhenotypeObservation(
            medium=A_not_B_medium,
            reaction_knockouts={},
            gene_knockouts=dict(),
            growth_phenotype=True,
        ),
        "B_not_A": PhenotypeObservation(
            medium=B_not_A_medium,
            reaction_knockouts={},
            gene_knockouts=dict(),
            growth_phenotype=True,
        ),
        "A_not_B_ko_AC": PhenotypeObservation(
            medium=A_not_B_medium,
            reaction_knockouts=reaction_knockouts,
            gene_knockouts=dict(),
            growth_phenotype=False,
        ),
        "B_not_A_ko_AC": PhenotypeObservation(
            medium=B_not_A_medium,
            reaction_knockouts=reaction_knockouts,
            gene_knockouts=dict(),
            growth_phenotype=True,
        ),
    }
    return phenotype_observations


@pytest.fixture
def growth_objective():
    S_index = ["A_int", "B_int", "C_int"]
    S_dict = {
        "A_SRC->A_int": [1, 0, 0],
        "A_int->B_int": [-1, 1, 0],
        "B_int->C_int": [0, -1, 1],
        "C_int->C_SNK": [0, 0, -1],
        "A_int->C_int": [-1, 0, 1],
        "B_SRC->B_int": [0, 1, 0],
    }
    S_table = pd.DataFrame(S_dict, index=S_index)
    mets, rxns = S_table.index, S_table.columns
    growth_objective = {rxn_id: 1 if rxn_id == "C_int->C_SNK" else 0 for rxn_id in rxns}
    return growth_objective


@pytest.fixture
def metabolite_dual():
    S_index = ["A_int", "B_int", "C_int"]
    S_dict = {
        "A_SRC->A_int": [1, 0, 0],
        "A_int->B_int": [-1, 1, 0],
        "B_int->C_int": [0, -1, 1],
        "C_int->C_SNK": [0, 0, -1],
        "A_int->C_int": [-1, 0, 1],
        "B_SRC->B_int": [0, 1, 0],
    }
    S_table = pd.DataFrame(S_dict, index=S_index)
    mets, rxns = S_table.index, S_table.columns
    metabolite_dual = {met_id: op.Variable(f"m_nogrowth_{met_id}") for met_id in mets}
    return metabolite_dual


@pytest.fixture
def flux_dual():
    S_index = ["A_int", "B_int", "C_int"]
    S_dict = {
        "A_SRC->A_int": [1, 0, 0],
        "A_int->B_int": [-1, 1, 0],
        "B_int->C_int": [0, -1, 1],
        "C_int->C_SNK": [0, 0, -1],
        "A_int->C_int": [-1, 0, 1],
        "B_SRC->B_int": [0, 1, 0],
    }
    S_table = pd.DataFrame(S_dict, index=S_index)
    mets, rxns = S_table.index, S_table.columns
    flux_dual = {rxn_id: op.Variable(f"r_nogrowth_{rxn_id}") for rxn_id in rxns}
    return flux_dual


@pytest.fixture
def flux():
    S_index = ["A_int", "B_int", "C_int"]
    S_dict = {
        "A_SRC->A_int": [1, 0, 0],
        "A_int->B_int": [-1, 1, 0],
        "B_int->C_int": [0, -1, 1],
        "C_int->C_SNK": [0, 0, -1],
        "A_int->C_int": [-1, 0, 1],
        "B_SRC->B_int": [0, 1, 0],
    }
    S_table = pd.DataFrame(S_dict, index=S_index)
    mets, rxns = S_table.index, S_table.columns
    flux = {rxn_id: op.Variable(f"v_growth_{rxn_id}") for rxn_id in rxns}
    return flux


def metabolite_equal(
    expected: cb.core.Metabolite,
    actual: cb.core.Metabolite,
    attributes: list[str] = ["id", "name", "compartment", "formula"],
) -> bool:
    return all(getattr(expected, attr) == getattr(actual, attr) for attr in attributes)


def reaction_equal(
    expected: cb.core.Reaction,
    actual: cb.core.Reaction,
    attributes: list[str] = ["id", "name", "lower_bound", "upper_bound", "objective_coefficient"],
) -> bool:
    return all(getattr(expected, attr) == getattr(actual, attr) for attr in attributes)


def constraint_equal(expected: op.Constraint, actual: op.Constraint) -> bool:
    """Checks if two Optlang constraint expressions are equal."""
    actual_expression = actual.to_json()["expression"]
    expected_expression = expected.to_json()["expression"]
    return actual_expression == expected_expression


### fixtures for testing

### copied from notebook - double check!!

# # stoichiometry
# S_index = ['A_int', 'B_int', 'C_int']
# S_dict = {'A_SRC->A_int': [1,0,0],
#           'A_int->B_int': [-1,1,0],
#           'B_int->C_int': [0,-1,1],
#           'C_int->C_SNK': [0,0,-1],
#           'A_int->C_int': [-1,0,1],
#           'B_SRC->B_int': [0,1,0],
#           }

# S_table = pd.DataFrame( S_dict ,index=S_index)


# mets, rxns = S_table.index, S_table.columns
# n_mets = len(mets)
# n_rxns = len(rxns)


# # constant Omega
# upper_flux_bound = 1e3

# # flux bounds
# min_growth_flux = 10
# max_nogrowth_flux = 5

# # media conditions
# A_not_B_medium = {
#  'A_SRC->A_int': min_growth_flux,
#  'B_SRC->B_int': 0.0
#  }

# B_not_A_medium = {
#  'A_SRC->A_int': 0.0,
#  'B_SRC->B_int': min_growth_flux
#  }

# # knockout reactions
# reaction_knockouts = {
#  'A_int->C_int': True,
#  }

# ### TODO - IN PROGRESS
# # list of U_nogrowth and U_growth (reaction size)
# U_growth = {rxn_id:upper_flux_bound for rxn_id in rxns}
# U_nogrowth = {rxn_id:upper_flux_bound for rxn_id in rxns}
# U_growth['A_SRC->A_int'] = min_growth_flux
# U_growth['B_SRC->B_int'] = 0
# growth_objective = {rxn_id:1 if rxn_id=='C_int->C_SNK' else 0 for rxn_id in rxns}
# upper_flux_bounds = {rxn_id:0 for rxn_id in rxns}
# lower_flux_bounds = {rxn_id:0 for rxn_id in rxns}

# phenotype_observations = {
#     'A_not_B': PhenotypeObservation(medium=A_not_B_medium, reaction_knockouts={}, gene_knockouts=dict(), growth_phenotype=True),
#     'B_not_A': PhenotypeObservation(medium=B_not_A_medium, reaction_knockouts={}, gene_knockouts=dict(), growth_phenotype=True),
#     'A_not_B_ko_AC':PhenotypeObservation(medium=A_not_B_medium, reaction_knockouts=reaction_knockouts, gene_knockouts=dict(), growth_phenotype=False),
#     'B_not_A_ko_AC':PhenotypeObservation(medium=B_not_A_medium, reaction_knockouts=reaction_knockouts, gene_knockouts=dict(), growth_phenotype=True)
# }

# cobra_model = model_from_stoich_matrix(S_table, "ABC_toy_model_3", obj=growth_objective, )

# # list of weights (reaction size)
# likelihoods = [0.5 for _ in rxns]

# # list of z (reaction size)
# reaction_indicator = {rxn_id:op.Variable(f"z_{rxn_id}", type="binary") for rxn_id in rxns}

# # flux dual, metabolite dual, and flux variables (replaced with function)
# flux_dual = {}
# flux = {}
# metabolite_dual = {}
# for observation_id, observation in phenotype_observations.items():
#     if observation.growth_phenotype: # primal problem when growth is observed
#         flux[observation_id] = {rxn_id:op.Variable(f"v_{observation_id}_{rxn_id}") for rxn_id in rxns}
#     else: # dual problem when no growth is observed
#         flux_dual[observation_id] = {rxn_id:op.Variable(f"r_{observation_id}_{rxn_id}") for rxn_id in rxns}
#         metabolite_dual[observation_id] = {met_id:op.Variable(f"m_{observation_id}_{met_id}") for met_id in mets}

# flux, flux_dual, metabolite_dual = generate_decision_variables(cobra_model, phenotype_observations)

# gibbs_like_constraint = get_steady_state_dual_constraint(cobra_model, phenotype_observations, growth_objective)

# # list of z (reaction size)
# reaction_indicator = {rxn_id:op.Variable(f"z_{rxn_id}", type="binary") for rxn_id in rxns}

# # list of r (no growth)
# flux_dual = {rxn_id:op.Variable(f"r_nogrowth_{rxn_id}") for rxn_id in rxns}

# # list of v (growth)
# flux = {rxn_id:op.Variable(f"v_growth_{rxn_id}") for rxn_id in rxns}

# # list of m (metabolite size)
# metabolite_dual = {met_id:op.Variable(f"m_nogrowth_{met_id}") for met_id in mets}

# # list of e_c (reaction size)

##### Fixtures #####

##### Tests #####
