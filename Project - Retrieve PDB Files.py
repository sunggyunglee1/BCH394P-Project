# Import bioservices module, to run remote UniProt queries
from bioservices import UniProt

# Import so that we can download PDB files
import urllib.request
import sys
import os
import socket

# Make a link to the UniProt webservice
service = UniProt()

# Build a query string - looking for EC 3.X.X.X
query = "ec:3"

# Define a list of columns we want to retrive - print(service._valid_columns)
PDB_columnlist = "xref_pdb"

# Run the remote search
PDB_result = service.search(query, frmt="tsv", columns=PDB_columnlist, limit=10)

# Strip extra characters and column header
PDB_result_strip = PDB_result.replace('\r', '').replace('\n', '').replace('PDB','')

# Make PDB ID list
PDB_ID_list = PDB_result_strip.split(';')

# Define directory to save PDB files
datadir = r"C:\Users\ameli\Documents\Bioinformatics\Final Project\PDB Files"

# Define a function to download PDB files
def download_pdb(PDB_ID_list, datadir, downloadurl = "https://files.rcsb.org/download/"):
    pdbfn = PDB_ID_list + ".pdb"
    url = downloadurl + pdbfn
    outfnm = os.path.join(datadir, pdbfn)
    socket.setdefaulttimeout(60)
    try:
        print(f"Downloading PDB file {PDB_ID_list} from {url} to {outfnm}...")
        urllib.request.urlretrieve(url, outfnm)
        print(f"PDB file {PDB_ID_list} downloaded successfully.")
        return outfnm
    except Exception as err:
        print(f"Error occurred while downloading PDB file {PDB_ID_list}: {err}", file=sys.stderr)
        return None

# Iterate over all PDB IDs in the list  
for id in PDB_ID_list:
    download_pdb(id, datadir)
    

