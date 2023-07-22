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
df = pl.read_csv("df_ai_cleaned.csv")
df = df.to_pandas()
df["mol"] = df["standard_smiles"].apply(lambda x: dm.to_mol(x))
mols = df["mol"]
mols = list(mols)


# Input--- **work-in-progress**
# TODO: input types ready, need to set up page style e.g. add columns etc.
app_ui = ui.page_fluid(
    ui.h4("Compound input selections"),
    ui.row(
        ui.column(3, ui.input_numeric("mol", "Index number of compound:", 0, min=0, max=143)),
        ui.column(3, ui.input_text("filename", "PNG file name for compound:")),
        ),
        ui.input_action_button("btn", "Confirm", class_="btn-success"),
        ui.output_image("image"),
    ui.row(
        ui.column(3, ui.input_numeric("mol1", "Specify index number of compound:", 0, min=0, max=143)),
        ui.column(3, ui.input_text("filename1", "PNG file name for compound:")),
        ),
        ui.input_action_button("btn1", "Confirm", class_="btn-success"),
        ui.output_image("image1")
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

        # Show image of saved PNG file of specified compound from input
        dir = Path(__file__).resolve().parent
        img: ImgData = {"src": str(dir / f"{input.filename()}.png")}
        return img
    
    @output
    @render.image
    def image1():
        # Place action button here to take a reactive dependency
        input.btn1()

        # Use reactive.isolate to not take a reactive dependency,
        # when entering/specifying compound and name of PNG file saved
        # only execute the PNG file generating & saving for the specified compound
        # when pressing the "Confirm" action button
        with reactive.isolate():
            Draw.MolToFile(mols[input.mol1()], f"{input.filename1()}.png")

        # Show image of saved PNG file of specified compound from input
        dir = Path(__file__).resolve().parent
        img1: ImgData = {"src": str(dir / f"{input.filename1()}.png")}
        return img1
    
    
    # def merge_image():


    #     img = Image.open(f"{input.filename()}.png")

    #     # --Using PIL/Pillow to manipulate images
    #     # img1 = Image.open("af1.png")
    #     # img2 = Image.open("af2.png")
    #     # img3 = Image.open("af3.png")
    #     # img4 = Image.open("af4.png")

    #     # Create a blank image template - (width, height)
    #     blank_image = Image.new("RGB", (600, 600))


    #     blank_image.paste(img, (0, 0))

    #     # Paste img1 to img4 together in a grid (top left, right & bottom left, right)
    #     # blank_image.paste(img1, (0, 0))
    #     # blank_image.paste(img2, (300, 0))
    #     # blank_image.paste(img3, (0, 300))
    #     # blank_image.paste(img4, (300, 300))


    #     blank_image.save(f"{input.merge_filename()}.png")

    #     # Save combined img1 & img2 as a new PNG file
    #     #blank_image.save("merged.png")


        
                


app = App(app_ui, server)