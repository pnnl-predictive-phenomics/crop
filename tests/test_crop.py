import pytest
import pandas as pd
import optlang as op


import crop
from crop import PhenotypeObservation
from crop import get_steady_state_dual_constraint
from crop import model_from_stoich_matrix

# fixtures for testing 

# copied from notebook - double check!!

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

### IN PROGRESS
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


# TODO: Run this
def test_get_steady_state_dual_constraint(model, phenotype_observations, growth_objective):
    actual_constraint = get_steady_state_dual_constraint(model, phenotype_observations, growth_objective)
    expected_constraint = {}
    assert actual_constraint == expected_constraint 


# tests needed

# dual flux variable constraint (r_nogrowth)

# steady state flux constraint (Sv=0)
    
# flux variable constraint (v_growth)
    
# reaction constraints (z = 1 for import reactions?)
