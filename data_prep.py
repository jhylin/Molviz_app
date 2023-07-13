# Example SMILES
# Canonical SMILES of calcitriol - CC(CCCC(C)(C)O)C1CCC2C1(CCCC2=CC=C3CC(CC(C3=C)O)O)C
# Canonical SMILES of ubiquinol - CC1=C(C(=C(C(=C1O)OC)OC)O)CC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C

# smiles = """CC(CCCC(C)(C)O)C1CCC2C1(CCCC2=CC=C3CC(CC(C3=C)O)O)C, 
#             CC1=C(C(=C(C(=C1O)OC)OC)O)CC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)CCC=C(C)C,
#             C(CC(=O)N)C(C(=O)O)N
# """
# sm = StringIO(smiles)
# df = pd.read_csv(sm, names = ["SMILES", "Names"])
# df["mol"] = df.SMILES.apply(Chem.MolFromSmiles)



# Import pandas
import pandas as pd
import polars as pl

#import ipywidgets as ipy
#from IPython.display import HTML, IFrame, display

# Add RDKit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.rdmolfiles import SmilesWriter, SmilesMolSupplier
from io import StringIO
import datamol as dm


# ***Specify data source***
# Reading from .smi file
# # Read in a simple list of SMILES from .smi file
# suppl = SmilesMolSupplier("cefe.smi")
# suppl

# # Convert a list of molecules into a dataframe
# mols = dm.to_df(suppl)
# mols

# # Generate a RDKit molecule column
# mols["mol"] = mols.smiles.apply(Chem.MolFromSmiles)
# # Show full dataframe - smiles & mol columns
# mols

# # Display first molecule in dataframe
# mols.iloc[0]["mol"]

# Function:
    # def image(filename):
    #     # Read in a simple list of SMILES from .smi file
    #     suppl = SmilesMolSupplier(input.filename)
    #     # Convert a list of molecules into a dataframe
    #     mols = dm.to_df(suppl)
    #     # Generate a RDKit molecule column
    #     mols["mol"] = mols.smiles.apply(Chem.MolFromSmiles)
    #     # Display first molecule in dataframe
    #     image = mols.iloc[0]["mol"]

    # return image



df = pl.read_csv("df_ai.csv")
#df.head()
df = df.to_pandas()
#df.head()
#type(df)
df["mol"] = df.Smiles.apply(Chem.MolFromSmiles)
df
mols = df["mol"]
mols