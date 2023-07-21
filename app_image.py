# App showing 2D images of small molecules of interests
# Rough plan of adding structure alignment/substructure highlighting once input determined

# Import libraries---
import pandas as pd
import polars as pl
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D

import datamol as dm

from pathlib import Path
from PIL import Image

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
# TODO: input types ready, need to set up page style e.g. add columns etc.
app_ui = ui.page_fluid(
    ui.input_numeric("mol", "Molecule index:", 0, min=0, max=143),
    ui.input_text("filename", "File name of PNG file:"),
    ui.input_action_button("btn", "Confirm", class_="btn-success"),
    ui.output_image("image")
)

# Output--- **work-in-progress**
# Below code works for inputting one molecule, saved as PNG file with image shown in app
def server(input, output, session):
    @output
    @render.image
    def image():
        # Place action button here to take a reactive dependency
        input.btn()

        # Using MolToFile() - only saving a single compound as PNG file (via index position)

        # Use reactive.isolate to not take a reactive dependency,
        # when entering/specifying compound and name of PNG file saved
        # only execute the PNG file generating & saving for the specified compound
        # when pressing the "Confirm" action button
        with reactive.isolate():
            Draw.MolToFile(mols[input.mol()], f"{input.filename()}.png")


        # TODO: ?can somehow formulate to take in 2 PNG images (or more e.g. 4)

        # --Testing MolToFile - 4 PNG files
        # Draw.MolToFile(mols[0], "af1.png")
        # Draw.MolToFile(mols[1], "af2.png")
        # Draw.MolToFile(mols[2], "af3.png")
        # Draw.MolToFile(mols[3], "af4.png")

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


        # Show image of saved PNG file of specified compound from input
        dir = Path(__file__).resolve().parent
        img: ImgData = {"src": str(dir / f"{input.filename()}.png")}
        return img
                


app = App(app_ui, server)