import os
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt

# function to calculate the Tanimoto similarity of two compounds
def tanimoto_similarity(mol1, mol2):
    # compute the fingerprint for each compound
    fp1 = AllChem.GetMorganFingerprint(mol1, 2)
    fp2 = AllChem.GetMorganFingerprint(mol2, 2)
    # return the Tanimoto similarity
    return rdkit.DataStructs.TanimotoSimilarity(fp1, fp2)

# list to store the compounds
compounds = []

# read in the SDF files from the given folder
folder = 'path/to/folder'
for file in os.listdir(folder):
    if file.endswith('.sdf'):
        # read in the compound from the SDF file
        compound = Chem.SDMolSupplier(os.path.join(folder, file))[0]
        # add the compound to the list
        compounds.append(compound)

# compute the Tanimoto similarity matrix
n = len(compounds)
matrix = [[0 for _ in range(n)] for _ in range(n)]
for i in range(n):
    for j in range(n):
        matrix[i][j] = tanimoto_similarity(compounds[i], compounds[j])

# generate the dendrogram
fig = plt.figure(figsize=(8, 8))
ax1 = fig.add_axes([0.09, 0.1, 0.2, 0.6])
Y = hierarchy.linkage(matrix, method='ward')
Z1 = hierarchy.dendrogram(Y, orientation='right')
ax1.set_xticks([])
ax1.set_yticks([])

# generate the heatmap
ax2 = fig.add_axes([0.3, 0.71, 0.6, 0.2])
ax2.matshow(matrix, aspect='auto', origin='lower', cmap=plt.cm.YlGnBu)
plt.colorbar(ax=ax2)

# show the plot
plt.show()
