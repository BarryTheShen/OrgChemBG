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

        # Generate 2D coordinates
        AllChem.Compute2DCoords(mol)

        # Generate MOL file name
        mol_name = f"compound_{i + 1}.mol"
        output_mol_path = os.path.join(output_directory, mol_name)

        # Write to MOL file
        with open(output_mol_path, 'w') as mol_file:
            mol_file.write(Chem.MolToMolBlock(mol))

        # Generate PNG image of carbon skeleton in 2D
        img = Draw.MolToImage(mol, size=(300, 300), kekulize=True, wedgeBonds=True)
        png_name = f"compound_{i + 1}.png"
        output_png_path = os.path.join(output_directory, png_name)

        # Save the image
        img.save(output_png_path)

        print(f"Successfully processed {smiles} -> {mol_name}, {png_name}")

    except Exception as e:
        print(f"Error processing {smiles}: {e}")