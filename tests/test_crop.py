import pytest
import pandas as pd
import optlang as op
import cobra as cb

import crop
from crop import PhenotypeObservation
from crop import get_steady_state_dual_constraint
from crop import model_from_stoich_matrix
from crop import generate_decision_variables


### fixtures for testing 

### copied from notebook - double check!!

# stoichiometry 
S_index = ['A_int', 'B_int', 'C_int']
S_dict = {'A_SRC->A_int': [1,0,0],
          'A_int->B_int': [-1,1,0],
          'B_int->C_int': [0,-1,1],
          'C_int->C_SNK': [0,0,-1],
          'A_int->C_int': [-1,0,1],
          'B_SRC->B_int': [0,1,0],
          }

S_table = pd.DataFrame( S_dict ,index=S_index)
mets, rxns = S_table.index, S_table.columns
n_mets = len(mets)
n_rxns = len(rxns)


# constant Omega
upper_flux_bound = 1e3

# flux bounds
min_growth_flux = 10
max_nogrowth_flux = 5

# media conditions
A_not_B_medium = {
 'A_SRC->A_int': min_growth_flux,
 'B_SRC->B_int': 0.0
 }

B_not_A_medium = {
 'A_SRC->A_int': 0.0,
 'B_SRC->B_int': min_growth_flux
 }

# knockout reactions
reaction_knockouts = {
 'A_int->C_int': True,
 }

### TODO - IN PROGRESS
# list of U_nogrowth and U_growth (reaction size)
U_growth = {rxn_id:upper_flux_bound for rxn_id in rxns}
U_nogrowth = {rxn_id:upper_flux_bound for rxn_id in rxns}
U_growth['A_SRC->A_int'] = min_growth_flux
U_growth['B_SRC->B_int'] = 0
growth_objective = {rxn_id:1 if rxn_id=='C_int->C_SNK' else 0 for rxn_id in rxns}
upper_flux_bounds = {rxn_id:0 for rxn_id in rxns}
lower_flux_bounds = {rxn_id:0 for rxn_id in rxns}

phenotype_observations = {
    'A_not_B': PhenotypeObservation(medium=A_not_B_medium, reaction_knockouts={}, gene_knockouts=dict(), growth_phenotype=True),
    'B_not_A': PhenotypeObservation(medium=B_not_A_medium, reaction_knockouts={}, gene_knockouts=dict(), growth_phenotype=True),
    'A_not_B_ko_AC':PhenotypeObservation(medium=A_not_B_medium, reaction_knockouts=reaction_knockouts, gene_knockouts=dict(), growth_phenotype=False),
    'B_not_A_ko_AC':PhenotypeObservation(medium=B_not_A_medium, reaction_knockouts=reaction_knockouts, gene_knockouts=dict(), growth_phenotype=True)
}

cobra_model = model_from_stoich_matrix(S_table, "ABC_toy_model_3", obj=growth_objective, )

# list of weights (reaction size)
likelihoods = [0.5 for _ in rxns]

# list of z (reaction size)
reaction_indicator = {rxn_id:op.Variable(f"z_{rxn_id}", type="binary") for rxn_id in rxns}

# flux dual, metabolite dual, and flux variables (replaced with function)
flux_dual = {}
flux = {}
metabolite_dual = {}
for observation_id, observation in phenotype_observations.items():
    if observation.growth_phenotype: # primal problem when growth is observed
        flux[observation_id] = {rxn_id:op.Variable(f"v_{observation_id}_{rxn_id}") for rxn_id in rxns}
    else: # dual problem when no growth is observed
        flux_dual[observation_id] = {rxn_id:op.Variable(f"r_{observation_id}_{rxn_id}") for rxn_id in rxns}
        metabolite_dual[observation_id] = {met_id:op.Variable(f"m_{observation_id}_{met_id}") for met_id in mets}

flux, flux_dual, metabolite_dual = generate_decision_variables(cobra_model, phenotype_observations)

gibbs_like_constraint = get_steady_state_dual_constraint(cobra_model, phenotype_observations, growth_objective)

# list of z (reaction size)
reaction_indicator = {rxn_id:op.Variable(f"z_{rxn_id}", type="binary") for rxn_id in rxns}

# list of r (no growth)
flux_dual = {rxn_id:op.Variable(f"r_nogrowth_{rxn_id}") for rxn_id in rxns}

# list of v (growth)
flux = {rxn_id:op.Variable(f"v_growth_{rxn_id}") for rxn_id in rxns}

# list of m (metabolite size)
metabolite_dual = {met_id:op.Variable(f"m_nogrowth_{met_id}") for met_id in mets}

# # list of e_c (reaction size)

##### Tests #####

# growth objective (e.g. C_int->C_SNK)
def test_growth_objective():
    actual = ... # get growth objective for ABC toy model 3
    expected = ...  # growth objective constraints for ABC toy model 3
    raise NotImplementedError


### COBRA model tests
def test_cobra_model_metabolites(model:cb.core.model.Model):
    # expect 'A_int' 'B_int' 'C_int'
    actual = ... # get list cobra model metabolite ids for ABC toy model 3
    expected = ...  # list of metabolite ids for ABC toy model 3
    raise NotImplementedError

def test_cobra_model_reactions(model:cb.core.model.Model):
    # expect 'A_SRC->int' 'A_int->B_int' 'B_int->C_int' 'C_int->C_SNK' 
    #        'A_int->C_int' 'B_SRC->B_int' 
    actual = ... # get list cobra model reactions for ABC toy model 3
    expected = ...  # list of reactions for ABC toy model 3
    raise NotImplementedError

def test_cobra_model_stoichiometric_matrix(model:cb.core.model.Model):
    # expect ...
    actual = ... # get cobra model S matrix for ABC toy model 3
    expected = ...  # S matrix for ABC toy model 3
    raise NotImplementedError


### MILP constraints
# test steady state dual constraint variables (Gibbs-like S^T*m + e_c = r_nogrowth)
def test_get_steady_state_dual_constraint(model, phenotype_observations, growth_objective):
    actual_constraint = get_steady_state_dual_constraint(model, 
                                                         phenotype_observations, 
                                                         growth_objective)
    expected_constraint = {}  # print out constraints and check by hand
    assert actual_constraint == expected_constraint 

# dual flux variable constraint (r_nogrowth_i < Omega*(1-z_i) ...)
def test_get_dual_flux_variable_constraint():
    actual = ... # get optlang dual flux constraint for ABC toy model 3
    expected = ...  # dual flux constraints for ABC toy model 3
    raise NotImplementedError

# steady state flux constraint (Sv_growth=0)
def test_steady_state_flux_constraint():
    actual = ... # get optlang steady state flux constraint for ABC toy model 3
    expected = ...  # steady state flux constraints for ABC toy model 3
    raise NotImplementedError 

# flux variable constraint (0<=v_growth_i<=U_growth_i*z_i ...)
def test_flux_constraint():
    actual = ... # get optlang flux constraints for ABC toy model 3
    expected = ...  # flux constraints for ABC toy model 3
    raise NotImplementedError
    
# reaction constraints (z = 1 for import reactions)
def test_reaction_constraints():
    actual = ... # get optlang reaction selector constraints for ABC toy model 3
    expected = ...  # reaction selector constraints for ABC toy model 3
    raise NotImplementedError

# mixed integer test
def test_mixed_integer_constraint():
    actual = ... # get Z type 
    expected = ...  # Z should be mixed integer (binary)
    raise NotImplementedError


### CROP predictions for ABC toy model 3
def test_CROP_prediction_ABCmodel3():
    """
    Test case 1:
    Current model has A-->B reaction
    True model doesn't have A-->B reaction

    Conditions
    1. observe growth with A in media but not B
    2. observe growth with B in media but not A
    3. observe no growth with A in media but not B - with a knockout for A-->C
    4. observe growth with B in the media but not A - with a knockout for A-->C

    We expect CROP to remove A-->B reaction (z=0 for A-->B reaction)
    """
    raise NotImplementedError

### CROP predictions for ABC toy model 3
def test_CROP_prediction_ABCmodel3_consistency():
    """
    Test case 1:
    Current model has A-->B reaction
    True model doesn't have A-->B reaction

    Conditions
    1. observe growth with A in media but not B
    2. observe growth with B in media but not A
    3. observe no growth with A in media but not B - with a knockout for A-->C
    4. observe growth with B in the media but not A - with a knockout for A-->C

    We expect CROP to remove A-->B reaction (z=0 for A-->B reaction)
    Furthermore, this solution should satisfy the constraints given
    """
    raise NotImplementedError

### Other tests

# test decision variable generation
def test_generate_decision_variables(model:cb.core.model.Model, phenotype_observations:PhenotypeObservation):
    # check that flux, flux dual, and metabolite dual variables are as generated expected
    raise NotImplementedError
    
# test flux bounds function(s)
def test_generate_flux_bounds_from_growth_observation(phenotype_observation:PhenotypeObservation):
    raise NotImplementedError
    
def test_generate_flux_dual_bounds_from_nogrowth_observation(phenotype_observation:PhenotypeObservation):
    raise NotImplementedError

# test PhenotypeObservation data class
def test_phenotype_observation():
    raise NotImplementedError

# test ConsistentReproductionOfPhenotype class
def test_ConsistentReproductionOfPhenotype_class():
    raise NotImplementedError

