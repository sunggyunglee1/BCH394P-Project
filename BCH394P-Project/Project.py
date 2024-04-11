# Standard library packages
import io

# Import bioservices module, to run remote UniProt queries
from bioservices import UniProt

# Import Pandas, so we can use dataframes
import pandas as pd

# Import PDB, so we can look at structures
import Bio.PDB.PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch

# Make a link to the UniProt webservice
service = UniProt()

# Build a query string - looking for EC 3.X.X.X
query = "ec:3"

# Define a list of columns we want to retrive - print(service._valid_columns)
columnlist = "accession,ec,xref_pdb,ft_act_site"

# Run the remote search
result = service.search(query, frmt="tsv", columns=columnlist, limit=20)
#print(result)

# Convert the last search result into a dataframe in Pandas
df = pd.read_table(io.StringIO(result))

# New df with NaN lines removed
df_new = df.dropna()
#print(df_new)

# use regular expression to pull out the active site values
extracted_numbers = df_new['Active site'].str.extractall(r'ACT_SITE\s+(\d+)')

# Group active site back into each row
#print(extracted_numbers)
grouped_numbers = extracted_numbers[0].groupby(level=0).apply(list)
#print(grouped_numbers)

# Rewrite active site column to just be numbers
df_new.loc[:, 'Active site'] = grouped_numbers

parser = Bio.PDB.PDBParser(QUIET=True) # QUIET=True avoids comments on errors in the pdb.
#parser = PDBParser() - OG one 
datadir = "C:\\Users\\ameli\\Documents\\Bioinformatics\\Final Project\\PDB Files\\"
for index, row in df_new.iterrows():
    #print(row['PDB'])
    pdb = row["PDB"]
    pdb_ids = pdb.split(';')
    for index, pdb_id in enumerate(pdb_ids):
        if index == 0: #so that it will only do the first pdb file in the list
            path = datadir + pdb_id + ".pdb"
            structures = parser.get_structure(pdb_id, path)
            structure = structures[0] # 'structures' may contain several proteins in this case only one.
            print(pdb_id)
    
            #Put for loop here
            sites = row["Active site"]
            for AS in sites:
            
                target_atom = structure['A'][int(AS)]['CA'] # I think this one is ok

                # Create a NeighborSearch object
                atoms  = Bio.PDB.Selection.unfold_entities(structure, 'A')
                ns = Bio.PDB.NeighborSearch(atoms)
                
                close_residues = ns.search(target_atom.coord, 3, level="R")
                
                # Print information about the nearby atoms
                print(pdb_id, AS, close_residues)
        else:
            pass
#RESOURCES
#https://www.biostars.org/p/269579/ - the VIP
#https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ
#https://bioservices.readthedocs.io/en/main/references.html#module-bioservices.uniprot