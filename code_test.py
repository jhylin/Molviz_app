# --Import libraries
# Import Pandas
import pandas as pd
import polars as pl

# Add RDKit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import PandasTools
from rdkit.Chem.Draw import rdMolDraw2D, MolsToGridImage
#from rdkit.Chem.rdmolfiles import SmilesWriter, SmilesMolSupplier

# **IPythonConsole for Jupyter notebook environment only**
from rdkit.Chem.Draw import IPythonConsole
# Set below to false to show PNG 
IPythonConsole.ipython_useSVG=False
# Set below to True to return SVG
#IPythonConsole.ipython_useSVG=True

import io
from PIL import Image
#from IPython.display import Image
from io import StringIO
import datamol as dm


# --Some code ideas:
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


# --Code test:
# Reading data from .csv file
df = pl.read_csv("df_ai.csv")
#df.head()
df = df.to_pandas()
#df.head()
#type(df)
# Best to use copy of df to avoid changing original dataset object
df = df.copy()


# Generate RDKit molecules as a column from dataframe
df["mol"] = df.Smiles.apply(Chem.MolFromSmiles)
df
mols = df["mol"]
mols
# Pandas series
type(mols)
mols = list(mols)
# Pandas list
type(mols)
# df.set_index(["Name"])
# df

# --Testing MolsToGridImage - gives IPython.core.display.Image object
# which means this'll likely only work for Jupyter notebook environment only
# image = Draw.MolsToGridImage(mols, molsPerRow=4, returnPNG=True)
# image


# --Testing MolsToImage - saving molecules as PNG file & open PNG file directly
# **Shorter code for simple mols to image function only**
# img_test = Draw.MolsToImage(mols)
# img_test
# img_test.save("antiinf.png")
# Potentially replacing below code with the PyShiny example
# of using file path to open PNG image
# image_new = Image.open("anti-inf.png")
# image_new.show()


# --Testing MolToFile 
# **Write a function to allow input of index number & file name to save 2 or more molecules in 1 file**
# MolToFile() can only save a single compound (specify index position) as PNG file

# Draft function v.1 - saving specified compound as PNG file:
# index = index position number of each compound
# file_name = name of PNG file to be saved

# def select_molecules(index, file_name):
#     for mol in df["mol"]:
#         image = Draw.MolToFile(mols[index], f"{file_name}.png")
#         return image
    
# select_molecules(0, "test")
# select_molecules(1, "test1")

# img = Image.open(f"{file_name}".png)
# blank_image = Image.new("RGB", (600, 600))

# #Draw.MolToFile(mols[2], "anti.png")
# Draw.MolToFile(mols[0], "af1.png")
# Draw.MolToFile(mols[1], "af2.png")

# # --Using PIL/Pillow to manipulate images
# img1 = Image.open("af1.png")
# img2 = Image.open("af2.png")

# blank_image = Image.new("RGB", (600, 300))

# blank_image.paste(img1, (0, 0))
# blank_image.paste(img2, (300, 0))
# blank_image.save("merged.png")


# --RDKit Cairo molecule drawer - saving molecules as PNG image file
# **Longer code but with other functions e.g. saving PNG data as string and others**
# Compounds stacked on top of each other in PNG initially (change the frame size parameters)
# Saving 2D compound image as PNG - sample frame size: 500,180,200,180
# Code below:
# drawer = rdMolDraw2D.MolDraw2DCairo(2000,2000,300,300) 
# drawer.drawOptions().useBWAtomPalette()
# drawer.DrawMolecules(mols)
# drawer.FinishDrawing()
# drawer.WriteDrawingText('anti-inf.png')

# Open the PNG file to show image
# image_test = Image.open("anti-inf.png")
# image_test.show()


# --RDKit SVG molecule drawer - produces a long string of SVG data
# drawer = rdMolDraw2D.MolDraw2DSVG(2000,2000,300,300) 
# drawer.drawOptions().useBWAtomPalette()
# drawer.DrawMolecules(mols)
# drawer.FinishDrawing()
# drawer.GetDrawingText()


# --Example of opening byte array data
# import io
# f = io.BytesIO(received_data)
# im = Image.open(f)


# --Trialled FrameToGridImage() - produces IPython.core.display.Image object
# rdkit.Chem.PandasTools.FrameToGridImage(frame, column='ROMol', legendsCol=None, **kwargs)
# df
# Chem.PandasTools.FrameToGridImage(df, column = "mol")


# TODO: Test lasso highlighting in Datamol and substructure highlighting in RDKit