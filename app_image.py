# Shiny for Python app for viewing and saving 2D images of small molecules of interests
# Showing 2D images of molecules via MolToFile()
# TODO: Adding substructure highlighting via MolToImage()

# Import libraries---
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

import datamol as dm

from pathlib import Path
from PIL import Image

from shiny import App, Inputs, Outputs, Session, render, ui, reactive
from shiny.types import ImgData


# Data source---
# Using pre-processed and standardised data 
df = pd.read_csv("df_ai_cleaned.csv")
# Avoid any changes to original dataset object by using .copy()
df = df.copy()
df["mol"] = df["standard_smiles"].apply(lambda x: dm.to_mol(x))
mols = df["mol"]
mols = list(mols)


# Input--- 
# TODO: formatting page style
app_ui = ui.page_fluid(
    ui.h4("Compound input selections"),
    ui.row(
        ui.column(
            3, 
            ui.input_numeric("mol1", "Select index number of 1st compound:", 0, min=0, max=143),
            ui.input_text("filename1", "File name for compound:"),
            ui.input_action_button("btn1", "Confirm", class_="btn"),
            ui.output_image("image1")
        ),
        ui.column(
            3, 
            ui.input_numeric("mol2", "Select index number of 2nd compound:", 0, min=0, max=143),
            ui.input_text("filename2", "File name for compound:"),
            ui.input_action_button("btn2", "Confirm", class_="btn"),
            ui.output_image("image2"),
        ),
        ui.column(
            3,
            ui.input_numeric("mol3", "Select index number of 3rd compound:", 0, min=0, max=143),
            ui.input_text("filename3", "File name for compound:"),
            ui.input_action_button("btn3", "Confirm", class_="btn"),
            ui.output_image("image3"),
        ),
        ui.column(
            3,
            ui.input_numeric("mol4", "Select index number of 4th compound:", 0, min=0, max=143),
            ui.input_text("filename4", "File name for compound:"),
            ui.input_action_button("btn4", "Confirm", class_="btn"),
            ui.output_image("image4"),
        )
    ),
    ui.row(
        ui.column(3, ui.input_text("merge_filename", "File name for merged images:")),
    ),
        ui.input_action_button("btn_merge", "Confirm", class_="btn"),
        ui.output_image("merge_image"),
    # Potentially adding MolToImage() to highlight atoms and bonds
    # ui.row(
    #     ui.column(3, ui.input_text("highlight_mol", "Enter:")),
    # ),
    #     ui.input_action_button("btn_highlight", "Send", class_="btn"),
    #     ui.output_image("highlight")
)

# Old app_ui layout:
# app_ui = ui.page_fluid(
#     ui.h4("Compound input selections"),
#     ui.row(
#         ui.column(3, ui.input_numeric("mol1", "Index number of compound:", 0, min=0, max=143)),
#         ui.column(3, ui.input_text("filename1", "File name for compound:")),
#         ),
#         ui.input_action_button("btn1", "Confirm", class_="btn"),
#         ui.output_image("image1"),
#     ui.row(
#         ui.column(3, ui.input_numeric("mol2", "Specify index number of compound:", 0, min=0, max=143)),
#         ui.column(3, ui.input_text("filename2", "File name for compound:")),
#         ),
#         ui.input_action_button("btn2", "Confirm", class_="btn"),
#         ui.output_image("image2"),
#     ui.row(
#         ui.column(3, ui.input_text("merge_filename", "File name for merged images:")),
#         ),
#         ui.input_action_button("btn_merge", "Confirm", class_="btn"),
#         ui.output_image("merge_image")
# )
        

# Server output--- 
def server(input, output, session):
    @output
    @render.image
    # Function to show 1st selected PNG image
    def image1():

        # Place action button here to take a reactive dependency
        input.btn1()

        # Use reactive.isolate to not take a reactive dependency,
        # but only when entering/specifying compound & name of PNG file saved
        # then execute the generating & saving of PNG file for the specified compound
        # after pressing "Confirm" action button
        with reactive.isolate():
            Draw.MolToFile(mols[input.mol1()], f"{input.filename1()}.png")

        # Show saved PNG file from selected compound
        dir = Path(__file__).resolve().parent
        img1: ImgData = {"src": str(dir / f"{input.filename1()}.png")}
        return img1
    
    @output
    @render.image
    # Function to show 2nd selected PNG image
    def image2():

        # 2nd action button
        input.btn2()

        # Use reactive.isolate() to isolate MolToFile() code until confirming with action button
        with reactive.isolate():
            Draw.MolToFile(mols[input.mol2()], f"{input.filename2()}.png")

        # Show saved PNG file from selected compound
        dir = Path(__file__).resolve().parent
        img2: ImgData = {"src": str(dir / f"{input.filename2()}.png")}
        return img2
    
    @output
    @render.image
    # Function to show 3rd selected PNG image
    def image3():

        # 3rd action button
        input.btn3()

        # Use reactive.isolate() to isolate MolToFile() code until confirming with action button
        with reactive.isolate():
            Draw.MolToFile(mols[input.mol3()], f"{input.filename3()}.png")

        # Show saved PNG file from selected compound
        dir = Path(__file__).resolve().parent
        img3: ImgData = {"src": str(dir / f"{input.filename3()}.png")}
        return img3
    
    @output
    @render.image
    # Function to show 4th selected PNG image
    def image4():

        # 4th action button
        input.btn4()

        # Use reactive.isolate() to isolate MolToFile() code until confirming with action button
        with reactive.isolate():
            Draw.MolToFile(mols[input.mol4()], f"{input.filename4()}.png")

        # Show saved PNG file from selected compound
        dir = Path(__file__).resolve().parent
        img4: ImgData = {"src": str(dir / f"{input.filename4()}.png")}
        return img4
    
    @output
    @render.image
    # Function to paste selected PNG files into a table of 4 images
    def merge_image():

        input.btn_merge()

        with reactive.isolate():
            # Using PIL/Pillow to manipulate images
            img1 = Image.open(f"{input.filename1()}.png")
            img2 = Image.open(f"{input.filename2()}.png")
            img3 = Image.open(f"{input.filename3()}.png")
            img4 = Image.open(f"{input.filename4()}.png")

            # Create a blank image template - (width, height)
            blank_image = Image.new("RGB", (600, 600))

            # Paste img1 to img4 together in a grid (top left, right, bottom left, right)
            blank_image.paste(img1, (0, 0))
            blank_image.paste(img2, (300, 0))
            blank_image.paste(img3, (0, 300))
            blank_image.paste(img4, (300, 300))

            # Save combined/merged image as a new PNG file
            blank_image.save(f"{input.merge_filename()}.png")
        
        # Show merged image of selected PNG images
        dir = Path(__file__).resolve().parent
        img_merge: ImgData = {"src": str(dir / f"{input.merge_filename()}.png")}
        return img_merge
    

    # TODO: Potentially adding MolToImage() to highlight bonds and atoms in compounds
    # @output
    # @render.image

    # Testing MolToImage() - code below worked in code_test.py
    # again this works for one molecule at a time via index position
    # from matplotlib.colors import ColorConverter 
    # img = Draw.MolToImage(mols[1], highlightAtoms=[1,2,3], highlightBonds = [1,2], highlightColor=ColorConverter().to_rgb("aqua")) 
    # img.save("molecule.png")
    

app = App(app_ui, server)