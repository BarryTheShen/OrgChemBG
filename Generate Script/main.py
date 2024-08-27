from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import os

# File paths
input_file = 'compound_list.txt'  # Input file containing SMILES strings, one per line
output_directory = 'output_files'  # Directory to save MOL and PNG files

# Create the output directory if it doesn't exist
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Function to count carbon atoms in a molecule
def count_carbons(mol):
    return sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')

# Function to remove hydrogens from a molecule
def remove_hydrogens(mol):
    return Chem.RemoveHs(mol)

# Read SMILES strings from the input file
with open(input_file, 'r') as f:
    smiles_list = f.readlines()

# Process each SMILES string
for i, smiles in enumerate(smiles_list):
    smiles = smiles.strip()
    try:
        # Convert SMILES to RDKit molecule object
        mol = Chem.MolFromSmiles(smiles)

        if mol is None:
            print(f"Error processing SMILES: {smiles}")
            continue



        # Count carbon atoms to decide image type
        carbon_count = count_carbons(mol)
        is_full_model = carbon_count < 3

        # Generate MOL file name
        mol_name = f"compound_{i + 1}.mol"
        output_mol_path = os.path.join(output_directory, mol_name)

        # Write to MOL file
        with open(output_mol_path, 'w') as mol_file:
            mol_file.write(Chem.MolToMolBlock(mol))

        if is_full_model:
            # Full Lewis structure for molecules with fewer than 3 carbons
            # Add hydrogens explicitly for accurate Lewis structure
            mol_with_hs = Chem.AddHs(mol)

            # Generate 2D coordinates
            AllChem.Compute2DCoords(mol_with_hs)
            options = Draw.MolDrawOptions()
            options.bondLineWidth = 3
            options.addAtomIndices = False  # Do not add atom indices
            options.includeAtomNumbers = False  # Do not show atom numbers
            options.explicitMethyl = True  # Explicitly show methyl groups (CH3)
            options.atomLabelFontSize = 12.0  # Font size for atom labels (adjust as needed)

            # Explicitly show carbon atoms with labels
            for atom in mol_with_hs.GetAtoms():
                if atom.GetSymbol() == 'C':
                    atom.SetProp('atomLabel', 'C')

            # Generate PNG image
            img = Draw.MolToImage(mol_with_hs, size=(600, 600), options=options)
            png_name = f"compound_{i + 1}.png"
            output_png_path = os.path.join(output_directory, png_name)

            # Save the image
            img.save(output_png_path)

        else:
            # Generate 2D coordinates
            AllChem.Compute2DCoords(mol)
            # Bond line structure for molecules with 3 or more carbons
            mol_skeleton = remove_hydrogens(mol)
            options = Draw.MolDrawOptions()
            options.bondLineWidth = 3
            options.addAtomIndices = False
            options.includeAtomNumbers = False
            options.explicitMethyl = False  # Do not explicitly show methyl groups
          #  options.

            # Generate PNG image of the bond line structure
            img = Draw.MolToImage(mol_skeleton, size=(600, 600), kekulize=True, wedgeBonds=True, options=options)
            png_name = f"compound_{i + 1}.png"
            output_png_path = os.path.join(output_directory, png_name)

            # Save the image
            img.save(output_png_path)

        print(f"Successfully processed {smiles} -> {mol_name}, {png_name}")

    except Exception as e:
        print(f"Error processing {smiles}: {e}")