# Shiny for Python app for viewing and saving 2D images of small molecules of interests
# 2D images of molecules saved as PNG files via MolToFile() 
# Atoms & bonds highlighting via MolToImage() 
# Select atom & bond highlighting (on or off) with or without atom indices


# TODO: 
# To add shinyswatch theme
# To add introduction of app - how to use, functions to view, save and look up compounds

# Import libraries---
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
import datamol as dm
from itables.shiny import DT
from pathlib import Path
from PIL import Image
from matplotlib.colors import ColorConverter 
from shiny import App, Inputs, Outputs, Session, render, ui, reactive
from shiny.types import ImgData



# Data source---
# Using pre-processed and standardised data 
df = pd.read_csv("df_ai_cleaned.csv")
# Avoid any changes to original dataset object by using .copy()
df = df.copy()
#df = df.reset_index()
df["mol"] = df["standard_smiles"].apply(lambda x: dm.to_mol(x))
mols = df["mol"]
mols = list(mols)


# Input--- 
# Added ui.navset_tab_card & ui.nav() to keep each of the 4 cpds in separate tabs! 
# - refer to "Orbit simulation" example
# Changed the merged PNG image from ui.row() to ui.column() to be in the right hand side of the space
# - refer to "Orbit simulation" example
# - to allow quick references of compound index numbers in the same app
# Removed multiple input_numeric fields for atom & bond highlightings 
# - using direct number inputs in input text areas then string to integer conversion
# Added interactive dataframe beneath PNG image processing area
# Added option to show atom indices in PNG images
# Changed tabs
# - 1 for non-indexed molecules & 1 for indexed molecules 
# - both with option to turn on/off highlighting

app_ui = ui.page_fluid(
    ui.h4("Compound input selections"),
    ui.row(
        ui.column(
            4,
            ui.navset_tab_card(
                ui.nav(
                    "Unindexed",
                    ui.input_numeric("mol1", "Select compound via index number:", 0, min = 0, max = 143),
                    ui.input_select("image_style1", "Choose substructure highlights:", 
                                    {"no_high": "Without highlights",
                                     "high": "With highlights"}),
                    ui.input_text_area("atom1", "Enter atom number to highlight:", placeholder = "e.g. 0, 1, 2, 3..."),
                    ui.input_text_area("bond1", "Enter bond number to highlight:", placeholder = "e.g. 0, 1, 2, 3..."),
                    ui.input_text("filename1", "File name for PNG image:"),
                    # Add caption explaining everytime confirm button is clicked, 
                    # a new PNG file will be saved with whatever file name was inputted above
                    # Make desired selections before clicking on the confirm button with diff file names for each
                    ui.input_action_button("btn1", "Save and view", class_= "btn"),
                    ui.output_image("image1")
                ),
                ui.nav(
                    "Indexed",
                    ui.input_numeric("mol2", "Select compound via index number:", 0, min = 0, max = 143),
                    ui.input_select("image_style2", "Choose substructure highlights:", 
                                    {"no_high": "Without highlights ", 
                                     "high": "With highlights"}),
                    ui.input_text_area("atom2", "Enter atom number to highlight", placeholder = "e.g. 0, 1, 2, 3..."),
                    ui.input_text_area("bond2", "Enter bond number to highlight:", placeholder = "e.g. 0, 1, 2, 3..."),
                    ui.input_text("filename2", "File name for PNG image:"),
                    ui.input_action_button("btn2", "Save and view", class_= "btn"),
                    ui.output_image("image2")
                ),
            )
        ),
        ui.column(
            6, 
            ui.input_text("merge_filename", "File name for merged PNG images:"),
            ui.input_action_button("btn_merge", "Confirm", class_="btn"),
            ui.output_image("merge_image"),
            ),
        ui.row(
            ui.page_fluid(ui.HTML(DT(df)))
        )
    ),
)


# Server output--- 
def server(input, output, session):

    @output
    @render.image
    @reactive.Calc
    # Function to show unindexed compound as PNG image
    def image1():

        input.btn1()

        with reactive.isolate():

            if input.image_style1() == "high":
                    # Highlight atoms & bonds
                    img = Draw.MolToImage(mols[input.mol1()], 
                                          highlightAtoms = [int(n) for n in input.atom1().split(",")],
                                          highlightBonds = [int(n) for n in input.bond1().split(",")], 
                                          highlightColor = ColorConverter().to_rgb("aqua")
                                          ) 
                    # Save image
                    img.save(f"{input.filename1()}.png")
                    # Show saved PNG file from selected compound
                    dir = Path(__file__).resolve().parent
                    img_1: ImgData = {"src": str(dir / f"{input.filename1()}.png")}
                    return img_1
                
            elif input.image_style1() == "no_high":
                    Draw.MolToFile(mols[input.mol1()], f"{input.filename1()}.png")
                    dir = Path(__file__).resolve().parent
                    img_1: ImgData = {"src": str(dir / f"{input.filename1()}.png")}
                    return img_1
       
            else:
                return None
    


    @output
    @render.image
    @reactive.Calc
    # Function to show indexed compound as PNG image
    def image2():

        input.btn2()

        with reactive.isolate():

            if input.image_style2() == "no_high":
                for atom in mols[input.mol2()].GetAtoms():
                    atom.SetProp('atomLabel', str(atom.GetIdx()))
                Draw.MolToFile(mols[input.mol2()], f"{input.filename2()}.png")
                dir = Path(__file__).resolve().parent
                img_2: ImgData = {"src": str(dir / f"{input.filename2()}.png")}
                return img_2
            
            elif input.image_style2() == "high":
                for atom in mols[input.mol2()].GetAtoms():
                    atom.SetProp('atomLabel', str(atom.GetIdx()))
                # Highlight atoms & bonds
                img = Draw.MolToImage(mols[input.mol2()], 
                                        highlightAtoms = [int(n) for n in input.atom2().split(",")],
                                        highlightBonds = [int(n) for n in input.bond2().split(",")], 
                                        highlightColor = ColorConverter().to_rgb("aqua")
                                        ) 
                # Save image
                img.save(f"{input.filename2()}.png")
                # Show saved PNG file from selected compound
                dir = Path(__file__).resolve().parent
                img_2: ImgData = {"src": str(dir / f"{input.filename2()}.png")}
                return img_2
            
            else:
                return None
    

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
    

app = App(app_ui, server)