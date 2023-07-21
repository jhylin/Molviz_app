# App showing 2D images of small molecules of interests
# Rough plan of adding structure alignment/substructure highlighting once input determined

# Import libraries---
import pandas as pd
import polars as pl
from rdkit import Chem
from rdkit.Chem import Draw
#from rdkit.Chem.Draw import IPythonConsole
#from rdkit.Chem.rdmolfiles import SmilesWriter, SmilesMolSupplier
from rdkit.Chem.Draw import rdMolDraw2D

import datamol as dm

from pathlib import Path
#import io
from PIL import Image
# IPython console similar to Jupyter notebook
# from IPython.display import Image

#from shiny import App, render, ui
from shiny import App, Inputs, Outputs, Session, render, ui, reactive
from shiny.types import ImgData

# Data source---
# Pandas library renders df columns in strange data types
# So using Polars instead then convert to Pandas df
df = pl.read_csv("df_ai.csv")
df = df.to_pandas()
df["mol"] = df["Smiles"].apply(lambda x: dm.to_mol(x))
mols = df["mol"]
mols = list(mols)


# Input--- **work-in-progress**
app_ui = ui.page_fluid(
    ui.input_numeric("mol", "Molecule index:", 0, min=0, max=143),
    ui.input_text("filename", "File name of PNG file:"),
    ui.input_action_button("btn", "Confirm", class_="btn-success"),
    ui.output_image("image")
)

# Output--- **work-in-progress**
def server(input, output, session):
    @output
    # Not requiring to render text as only images are needed as final output
    # If using both render.text & render.image at the same time,
    # only PNG file of specified molecule will be saved with no image shown!
    #@render.text
    @render.image
    def image():
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
        # For image of list of molecules: ?improve image resolution (blurry at the moment)

        # --Testing MolsToImage
        #Image of molecules all aligned horizontally (not in grid/table format)
        # img_test = Draw.MolsToImage(mols, fitImage = True)
        # img_test
        # img_test.save("antiinf.png")
        # image_new = Image.open("antiinf.png")
        # image_new.show()

        # Place action button here to take a reactive dependency
        input.btn()

        # MolToFile() - only saving a single compound as PNG file (via index position)
        # Draft function v.1 - saving specified compound as PNG file:
        # index = index position number of each compound
        # file_name = name of PNG file to be saved

        # **this part works with PNG file saved as input number**
        # Use reactive.isolate to not take a reactive dependency 
        # when entering/specifying compound and name of PNG file saved
        with reactive.isolate():
            Draw.MolToFile(mols[input.mol()], f"{input.filename()}.png")

        
        # --Testing MolToFile 
        # Draw.MolToFile(mols[0], "af1.png")
        # Draw.MolToFile(mols[1], "af2.png")
        # Draw.MolToFile(mols[2], "af3.png")
        # Draw.MolToFile(mols[3], "af4.png")


        #img = Image.open(f"{input.filename()}.png")

        # --Using PIL/Pillow to manipulate images
        # img1 = Image.open("af1.png")
        # img2 = Image.open("af2.png")
        # img3 = Image.open("af3.png")
        # img4 = Image.open("af4.png")

        # Create a blank image template - (width, height)
        #blank_image = Image.new("RGB", (600, 600))

        # Paste img1 to img4 together in a grid (top left, right & bottom left, right)
        # blank_image.paste(img1, (0, 0))
        # blank_image.paste(img2, (300, 0))
        # blank_image.paste(img3, (0, 300))
        # blank_image.paste(img4, (300, 300))

        # Save combined img1 & img2 as a new PNG file
        #blank_image.save("merged.png")


        # --Use local file path to import PNG image file
        # Showing a list of molecules
        # dir = Path(__file__).resolve().parent
        # img: ImgData = {"src": str(dir / "anti-inf.png"), "width": "2350px", "height": "2350px"}
        # return img

        # TODO: **Not showing up in app yet!** ?to add reactive effect
        # Showing image from the saved PNG file
        
        dir = Path(__file__).resolve().parent
        img: ImgData = {"src": str(dir / f"{input.filename()}.png")}
        return img
                

        #return image


app = App(app_ui, server)