import cobra
from pprint import pprint
from cobra.io import write_sbml_model, validate_sbml_model
from cobra.util.array import create_stoichiometric_matrix
import pandas as pd


pd.options.display.float_format = '{:.0f}'.format

mets = ['A', 'B', 'C', 'D', 'E', 'F', 'A_ext', 'C_ext', 'D_ext', 'E_ext',
        'F_ext']

internal_mets = [met for met in mets if 'ext' not in met]
rxns = ['R_{}'.format(i) for i in range(1, 11)]

data = {'R_1': pd.Series({'A': 1, 'A_ext': -1}),
        'R_2': pd.Series({'A': -1, 'B': 1}),
        'R_3': pd.Series({'A': -1, 'C': 1}),
        'R_4': pd.Series({'B': -1, 'D': 2, 'E': -1}),
        'R_5': pd.Series({'E': 1, 'E_ext': -1}),
        'R_6': pd.Series({'B': -2, 'C': 1, 'F': 1}),
        'R_7': pd.Series({'C': -1, 'D': 1}),
        'R_8': pd.Series({'D': -1, 'D_ext': 1}),
        'R_9': pd.Series({'F': -1, 'F_ext': 1}),
        'R_10': pd.Series({'C': -1, 'C_ext': 1})}

model = cobra.Model('example_model')


def add_importer(source, target, name):

    rxn = cobra.Reaction(f'im_{name}')
    rxn.add_metabolites({
        source: -1.0,
        target: 1.0,
    })
    return rxn


def add_metabolite(name):
    metabolite = cobra.Metabolite(f'{name}_c', compartment='c')
    ext_metabolite = cobra.Metabolite(f'{name}_e', compartment='e')
    rxn = add_importer(ext_metabolite, metabolite, name)

    model.add_metabolites([metabolite, ext_metabolite])
    model.add_boundary(ext_metabolite, type="exchange")
    model.add_reactions([rxn])
    return [metabolite, ext_metabolite]


A, A_e = add_metabolite('A')
B, B_e = add_metabolite('B')
C, C_e = add_metabolite('C')
D, D_e = add_metabolite('D')
E, E_e = add_metabolite('E')
F, F_e = add_metabolite('F')


model.add_boundary(model.metabolites.get_by_id("F_e"), type="sink")


r1 = cobra.Reaction('R_A_to_B')
r1.add_metabolites({
    A: -1.0,
    B: 1.0,
})

r2 = cobra.Reaction('R_BE_to_D')
r2.add_metabolites({
    B: -1.0,
    E: -1.0,
    D: 2.0,
})

r3 = cobra.Reaction('R_A_to_C')
r3.add_metabolites({
    A: -1.0,
    C: 1.0,
})

r4 = cobra.Reaction('R_BC_to_F')
r4.add_metabolites({
    B: -2.0,
    C: 1.0,
    F: 1
})

r5 = cobra.Reaction('R_C_to_D')
r5.add_metabolites({
    C: -1.0,
    D: 1.0,
})

r6 = cobra.Reaction('R_D_to_D_ex')
r6.add_metabolites({
    D: -1.0,
    D_e: 1.0,
})

r7 = cobra.Reaction('R_F_to_F_ex')
r7.add_metabolites({
    F: -1.0,
    F_e: 1.0,
})
model.add_reactions([r1, r2, r3, r4, r5, r6, r7])
model.objective = model.reactions.get_by_id('EX_D_e')

write_sbml_model(model, filename='example_model.xml')
report = validate_sbml_model(filename='example_model.xml')
pprint(report)

for rxn in model.reactions:
    print(rxn)
quit()


stoich_matrix = create_stoichiometric_matrix(model, 'DataFrame')

order = ['im_A', 'R_A_to_B', 'R_A_to_C', 'R_BE_to_D', 'im_E', 'R_BC_to_F',
         'R_C_to_D', 'im_D', 'SK_F_e', 'im_C', ]
fullS = pd.DataFrame(data, columns=rxns, index=mets, dtype='int64').fillna(0)
fullS.columns = order


stoich_matrix = stoich_matrix[order]

