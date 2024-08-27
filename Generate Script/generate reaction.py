from rdkit import Chem
from rdkit.Chem import AllChem, Draw

# Define a simple reaction using SMARTS
reaction = AllChem.ReactionFromSmarts('[OH:1][C:2]=[O:3].[CH3:4][I:5]>>[O:1][C:2]=[O:3][CH3:4]')

# Create a new reaction with explicit hydrogens
new_reaction = AllChem.ChemicalReaction()

"""for mol in reaction.GetReactants():
    mol = Chem.AddHs(mol)  # Add explicit hydrogens
    new_reaction.AddReactantTemplate(mol)

for mol in reaction.GetProducts():
    mol = Chem.AddHs(mol)  # Add explicit hydrogens
    new_reaction.AddProductTemplate(mol)"""

# Draw the reaction with explicit hydrogens
reaction_img = Draw.ReactionToImage(new_reaction)
reaction_img.show()

# Save the image to a file
reaction_img.save('reaction_with_hydrogens.png')
