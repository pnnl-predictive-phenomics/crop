# code formulation
from dataclasses import dataclass
import cvxpy as cp
import numpy as np
import pandas as pd
import optlang as op
import cobra as cb 
from cobra.core.reaction import Reaction
from cobra.core.metabolite import Metabolite
from optlang import Model, Variable, Constraint, Objective
import json


@dataclass
class PhenotypeObservation:
    """Specifies media conditions, gene knockouts, reaction knockouts, and growth/no growth observation."""
    medium:dict[str,float]
    gene_knockouts:dict[str,bool]  # TODO
    reaction_knockouts:dict[str,bool]
    growth_phenotype:bool

# TEST
def model_from_stoich_matrix(S:pd.DataFrame, name:str, obj:dict[str,int], lower_flux_bound:dict[str,float], upper_flux_bound:dict[str,float])->cb.core.model.Model:
    """Creates a cobra model from a stoichiometric matrix."""
    model = cb.core.model.Model(name)
    model.add_reactions([Reaction(rxn_id, lower_bound=lower_flux_bound[rxn_id], upper_bound=upper_flux_bound[rxn_id], obj=obj[rxn_id]) for rxn_id in S.columns])
    for rxn in model.reactions:
        rxn.add_metabolites({Metabolite(met_id):stoichiometry for met_id,stoichiometry in S[rxn.id].to_dict().items() if stoichiometry !=0})
    return model

# TODO: finish and test
class ConsistentReproductionOfPhenotype:
    """Adding and removing reactions based on phenotype observations."""

    def __init__(self, model:cb.core.model.Model, phenotypes:dict[str,PhenotypeObservation]):
        self.model = model
        self.phenotypes = phenotypes

    def create_problem(self):
        """Create an OptLang representation of the CROP problem."""
        pass

# TODO: finish and test
def generate_flux_bounds_from_growth_observation(phenotype_observation:PhenotypeObservation):
    # go through media and get constraints
    # go through knockouts
    pass

# TODO: finish and test
def generate_flux_dual_bounds_from_nogrowth_observation(phenotype_observation:PhenotypeObservation):
    #
    # go through knockouts
    pass

# TEST
def get_steady_state_dual_constraint(model, phenotype_observations, growth_objective, metabolite_dual, flux_dual):
    """Gets the steady state dual (Gibbs-like) constraints for the dual problem"""
    constraint = {}
    for observation_id, observation in phenotype_observations.items():
        if not observation.growth_phenotype:
            for reaction in model.reactions:
                constraint[f"Gibbs_{reaction.id}"] = Constraint(
                    sum(stoichiometry*metabolite_dual[observation_id][met_id] 
                        for met_id,stoichiometry in reaction.metabolites
                       )-flux_dual[observation_id][reaction.id] + growth_objective[reaction.id], 
                    lb=0, ub=0)
    return constraint

# TEST
def generate_decision_variables(model, phenotype_observations):
    # flux dual, metabolite dual, and flux variables
    flux_dual = {}
    flux = {}
    metabolite_dual = {}
    for observation_id, observation in phenotype_observations.items():
        if observation.growth_phenotype: # primal problem when growth is observed
            flux[observation_id] = {rxn.id:op.Variable(f"v_{observation_id}_{rxn.id}") for rxn in model.reactions}
        else: # dual problem when no growth is observed
            flux_dual[observation_id] = {rxn.id:op.Variable(f"r_{observation_id}_{rxn.id}") for rxn in model.reactions}
            metabolite_dual[observation_id] = {met.id:op.Variable(f"m_{observation_id}_{met.id}") for met in model.metabolites}