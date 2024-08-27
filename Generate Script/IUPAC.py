import pubchempy as pcp
import os



# Define the file names
input_file = "compound_list.txt"
output_file = "iupac_names.txt"

# Check if the input file exists
if not os.path.exists(input_file):
    print(f"File '{input_file}' not found.")
else:
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Strip any whitespace characters and ignore empty lines
            smiles = line.strip()
            if smiles:
                try:
                    # Search PubChem for the compound using SMILES
                    compound = pcp.get_compounds(smiles, 'smiles')

                    if compound:
                        # Get the IUPAC name
                        iupac_name = compound[0].iupac_name
                        if iupac_name:
                            print(f"Successfully converted {iupac_name}")
                            outfile.write(f"{iupac_name}\n")
                        else:
                            outfile.write("IUPAC name not found\n")
                    else:
                        outfile.write("Compound not found\n")
                except Exception as e:
                    outfile.write(f"Error: {str(e)}\n")

    print(f"IUPAC names have been saved to '{output_file}'.")
