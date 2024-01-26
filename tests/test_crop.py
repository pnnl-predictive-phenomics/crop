import cobra as cb
import optlang as op
from optlang import Variable
import pandas as pd
import pytest

import crop
from crop import (
    PhenotypeObservation,
    generate_decision_variables,
    get_steady_state_dual_constraint,
    model_from_stoich_matrix,
)

from .fixtures import (
    expected_model,
    actual_model,
    phenotype_observations,
    growth_objective,
    metabolite_dual,
    flux_dual,
    metabolite_equal,
    reaction_equal,
    constraint_equal
)


# steps
# 1. port ABC toy model into COBRA model and test it (metabolites, S, etc)
# 2. load and validate 'observations' using dataclass
# 3. generate and validate constraints for CROP problem
# 4. generate and validate predictions from CROP

# note: might be easier to start w/ one condition and FBA constraints first
#       then move to multiple conditions and full MILP (2 test cases)



### COBRA model tests
def test_cobra_model( expected_model:cb.core.model.Model, actual_model:cb.core.model.Model):
    """Test that actual cobra model matches expected cobra model - for ABC toy model 3"""

    # test metabolites match
    for expected_metabolite in expected_model.metabolites:
        assert metabolite_equal(expected_metabolite,actual_model.metabolites.get_by_id(expected_metabolite.id) )
    
    # test reactions match
    for expected_reaction in expected_model.reactions:
        assert reaction_equal(expected_reaction,actual_model.reactions.get_by_id(expected_reaction.id) )
    
    # check stoichiometry
    pd.testing.assert_frame_equal(
        cb.util.array.create_stoichiometric_matrix(actual_model, array_type='DataFrame'),
        cb.util.array.create_stoichiometric_matrix(expected_model, array_type='DataFrame')
    )


def test_get_steady_state_dual_constraints(expected_model,
                                          actual_model, 
                                          phenotype_observations, 
                                          growth_objective, 
                                          metabolite_dual, 
                                          flux_dual
                                          ):
    """Test getting the dual steady state flux variable constraint (i.e., :math:`S^Tm + e_{C\rightarrow} = r_{nogrowth}`)"""
    raise NotImplementedError


def test_get_dual_flux_variable_constraints(expected_model,
                                          actual_model, 
                                          phenotype_observations, 
                                          growth_objective, 
                                          flux_dual
                                          ):
    """Test getting the dual flux constraints (i.e., LB <= r_nogrowth <= UB)"""
    raise NotImplementedError


def test_get_steady_state_flux_constraints(expected_model,
                                          actual_model, 
                                          phenotype_observations, 
                                          growth_objective, 
                                          flux
                                          ):
    """Test getting the steady state flux constraints (i.e., Sv=0) """
    raise NotImplementedError

def test_get_flux_constraints(expected_model,
                                actual_model, 
                                phenotype_observations, 
                                growth_objective, 
                                flux
                                ):
    """Test getting the flux constraints (i.e., LB <= v_growth <= UB) """
    raise NotImplementedError


# z = 1 for A, B
def test_get_fixed_reaction_constraints(expected_model,
                                          actual_model, 
                                          phenotype_observations, 
                                          growth_objective, 
                                          flux
                                        ):
    "Testing getting the fixed reaction constraints (i.e., z=1 for A_SRC->A_int)"
    raise NotImplementedError


def get_CROP_predictions(expected_model,
                        actual_model, 
                        phenotype_observations, 
                        growth_objective
                        ):
    "Test that CROP removes the correct reactions for the ABC toy model 3 system"
    raise NotImplementedError




# # steady state flux constraint (Sv_growth=0)
# def test_steady_state_flux_constraint():
#     actual = ...  # get optlang steady state flux constraint for ABC toy model 3
#     expected = ...  # steady state flux constraints for ABC toy model 3
#     raise NotImplementedError


# # flux variable constraint (0<=v_growth_i<=U_growth_i*z_i ...)
# def test_flux_constraint():
#     actual = ...  # get optlang flux constraints for ABC toy model 3
#     expected = ...  # flux constraints for ABC toy model 3
#     raise NotImplementedError


# # reaction constraints (z = 1 for import reactions)
# def test_reaction_constraints():
#     actual = ...  # get optlang reaction selector constraints for ABC toy model 3
#     expected = ...  # reaction selector constraints for ABC toy model 3
#     raise NotImplementedError


# # mixed integer test
# def test_mixed_integer_constraint():
#     actual = ...  # get Z type
#     expected = ...  # Z should be mixed integer (binary)
#     raise NotImplementedError


# ### CROP predictions for ABC toy model 3
# def test_CROP_prediction_ABCmodel3():
#     """
#     Test case 1:
#     Current model has A-->B reaction
#     True model doesn't have A-->B reaction

#     Conditions
#     1. observe growth with A in media but not B
#     2. observe growth with B in media but not A
#     3. observe no growth with A in media but not B - with a knockout for A-->C
#     4. observe growth with B in the media but not A - with a knockout for A-->C

#     We expect CROP to remove A-->B reaction (z=0 for A-->B reaction)
#     """
#     raise NotImplementedError


# ### CROP predictions for ABC toy model 3
# def test_CROP_prediction_ABCmodel3_consistency():
#     """
#     Test case 1:
#     Current model has A-->B reaction
#     True model doesn't have A-->B reaction

#     Conditions
#     1. observe growth with A in media but not B
#     2. observe growth with B in media but not A
#     3. observe no growth with A in media but not B - with a knockout for A-->C
#     4. observe growth with B in the media but not A - with a knockout for A-->C

#     We expect CROP to remove A-->B reaction (z=0 for A-->B reaction)
#     Furthermore, this solution should satisfy the constraints given
#     """
#     raise NotImplementedError


# ### Other tests


# # test decision variable generation
# def test_generate_decision_variables(
#     model: cb.core.model.Model, phenotype_observations: PhenotypeObservation
# ):
#     # check that flux, flux dual, and metabolite dual variables are as generated expected
#     raise NotImplementedError


# # test flux bounds function(s)
# def test_generate_flux_bounds_from_growth_observation(phenotype_observation: PhenotypeObservation):
#     raise NotImplementedError


# def test_generate_flux_dual_bounds_from_nogrowth_observation(
#     phenotype_observation: PhenotypeObservation,
# ):
#     raise NotImplementedError


# # test PhenotypeObservation data class
# def test_phenotype_observation():
#     raise NotImplementedError


# # test ConsistentReproductionOfPhenotype class
# def test_ConsistentReproductionOfPhenotype_class():
#     raise NotImplementedError
