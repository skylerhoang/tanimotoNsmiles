import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

# Parse the SMILES strings and convert them into molecular structures
mol_1 = Chem.MolFromSmiles('SMILES_STRING_1')
mol_2 = Chem.MolFromSmiles('SMILES_STRING_2')

# Calculate the Tanimoto similarity between the two compounds
tanimoto_sim = AllChem.GetTanimotoSimilarity(mol_1, mol_2)

print('Tanimoto similarity:', tanimoto_sim)