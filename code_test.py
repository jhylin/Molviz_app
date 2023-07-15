# --Import libraries
# Import pandas
import pandas as pd
import polars as pl

#import ipywidgets as ipy
#from IPython.display import HTML, IFrame, display

# Add RDKit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D, MolsToGridImage
#from rdkit.Chem.rdmolfiles import SmilesWriter, SmilesMolSupplier

from rdkit.Chem.Draw import IPythonConsole
# Set below to False to return PNG 
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

# Generate RDKit molecules
df["mol"] = df.Smiles.apply(Chem.MolFromSmiles)
df
mols = df["mol"]
mols
type(mols)
mols = list(mols)
type(mols)

# --Testing MolsToGridImage - gives IPython.core.display.Image object
# image = Draw.MolsToGridImage(mols, molsPerRow=4, returnPNG=True)
# image

# --Testing MolsToImage - saving molecules as PNG file & open PNG directly
# **Shorter code for simple mols to image function only**
img_test = Draw.MolsToImage(mols)
img_test
img_test.save("antiinf.png")
# Potentially replacing below code with the PyShiny example
# of using file path to open PNG image
image_new = Image.open("anti-inf.png")
image_new.show()

# --TODO:
# Try ui.output_image and @render.image from PyShiny
# to show PNG image from file path - code example bookmarked
# which hopefully will show 2D image of compounds in PyShiny


# --Testing MolToFile
# Only can save a single compound (specify index position) as PNG
# Draw.MolToFile(mols[1], "anti.png")


# --RDKit Cairo molecule drawer - saving molecules as PNG image file
# **Longer code but with other functions e.g. saving PNG data as string and others**
# All compound stacked on top of each other in PNG (?because >50 compounds)
# Saving 2D compound image as PNG - default frame size: 500,180,200,180
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


