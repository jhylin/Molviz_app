# Import libraries---
import pandas as pd
import polars as pl
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
#from rdkit.Chem.rdmolfiles import SmilesWriter, SmilesMolSupplier
from rdkit.Chem.Draw import rdMolDraw2D

import datamol as dm

from pathlib import Path
import io
from PIL import Image
# IPython console similar to Jupyter notebook
# from IPython.display import Image

#from shiny import App, render, ui
from shiny import App, Inputs, Outputs, Session, render, ui
from shiny.types import ImgData


# Input--- **work-in-progress**
app_ui = ui.page_fluid(
    ui.output_image("image")
)

# Output---
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
        # Pandas library renders df columns in strange data types
        # So using Polars instead then convert to Pandas df
        df = pl.read_csv("df_ai.csv")
        df = df.to_pandas()
        df["mol"] = df["Smiles"].apply(lambda x: dm.to_mol(x))
        mols = df["mol"]
        mols = list(mols)

        # MolsToGridImage() works in Jupyter notebook only (returns IPython.core.display.Image object)
        #img = Draw.MolsToGridImage(mols) 

        # --Testing MolDraw2DCairo
        # Saving 2D compound image as PNG
        # drawer = rdMolDraw2D.MolDraw2DCairo(2000,2000,300,300)
        # Change to black & white
        # drawer.drawOptions().useBWAtomPalette()
        # drawer.DrawMolecules(mols)
        # drawer.FinishDrawing()
        # drawer.WriteDrawingText('anti-inf.png')

        # Open the PNG file to show image (for MacOS - opens file using Preview)
        # image_test = Image.open("anti-inf.png")
        # image_test.show()

        # --Testing MolsToImage
        #Image of molecules all aligned horizontally (not in grid/table format)
        # img_test = Draw.MolsToImage(mols, fitImage = True)
        # img_test
        # img_test.save("antiinf.png")
        # image_new = Image.open("antiinf.png")
        # image_new.show()

        # --TODO:
        # For image of list of molecules: ?improve image resolution (blurry at the moment)
        
        # Trial single molecule image, but write function to add e.g. 2 or 3 molecules as PNG file
        # then show image of selection of molecules (linking to input selection?)

        # **Function to allow input selection to save 2 or more molecules in 1 file**
        # Can only save a single compound (specify index position) as PNG file
        #def select_molecule(i):
        # --Testing MolToFile 
        Draw.MolToFile(mols[0], "af1.png")
        Draw.MolToFile(mols[1], "af2.png")

        # Consider using Image.paste() to place two pngs together in one image for comparison
        


        # --Use local file path to import PNG image file

        # Showing a list of molecules
        # dir = Path(__file__).resolve().parent
        # img: ImgData = {"src": str(dir / "anti-inf.png"), "width": "2350px", "height": "2350px"}
        # return img

        # Showing a single molecule
        dir = Path(__file__).resolve().parent
        img: ImgData = {"src": str(dir / "af1.png")} 
        return img 

        #return image


app = App(app_ui, server)