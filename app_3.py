import pandas as pd
import polars as pl
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
#from rdkit.Chem.rdmolfiles import SmilesWriter, SmilesMolSupplier
from rdkit.Chem.Draw import rdMolDraw2D
# Set below to false to show PNG 
IPythonConsole.ipython_useSVG=False

import io

from PIL import Image
# from IPython.display import Image

import datamol as dm
from shiny import App, render, ui

app_ui = ui.page_fluid(
    ui.output_image("image"),
)


def server(input, output, session):
    @output
    @render.image

    # @render.table
    # def table():
    #     #df = pd.read_csv("df_ai.csv")
    #     df = pl.read_csv("df_ai.csv")
    #     df = df.to_pandas()
    #     df["mol"] = df["Smiles"].apply(lambda x: dm.to_mol(x))
    #     mols = df["mol"]
    #     mols = list(mols)
    #     # image = dm.viz.to_image(mols)
    #     Draw.MolsToGridImage(mols)
    #     return table

    def image():
        #df = pd.read_csv("df_ai.csv")
        df = pl.read_csv("df_ai.csv")
        df = df.to_pandas()
        df["mol"] = df["Smiles"].apply(lambda x: dm.to_mol(x))
        mols = df["mol"]
        mols = list(mols)

        #img = Draw.MolsToGridImage(mols) 

        # --Testing MolDraw2DCairo
        # Saving 2D compound image as PNG
        # drawer = rdMolDraw2D.MolDraw2DCairo(500,180,200,180)
        # drawer.drawOptions().useBWAtomPalette()
        # drawer.DrawMolecules(mols)
        # drawer.FinishDrawing()
        # drawer.WriteDrawingText('anti-inf.png')

        # Open the PNG file to show image
        # image_test = Image.open("anti-inf.png")
        # image_test.show()

        # --Testing MolsToImage
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

        return image


app = App(app_ui, server)