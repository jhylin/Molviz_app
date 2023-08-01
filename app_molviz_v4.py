# Shiny for Python app for viewing and saving 2D images of small molecules of interests
# 2D images of molecules saved as PNG files via MolToFile() 
# Atoms & bonds highlighting via MolToImage() 
# Select atom & bond highlighting (on or off) with or without atom indices


# TODO: 
# To add introduction of app 
# - how to use, features to view, save and look up compounds
# - split sections: Panel to add compound, panel to view merged image, df to search below

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
import shinyswatch



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
# ***Added features***
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
# - 2 for non-indexed molecules & 2 for indexed molecules 
# - both with option to turn on/off highlighting
# - to allow quick comparison between unindexed & indexed molecules side-by-side
# Added shinyswatch theme solar

app_ui = ui.page_fluid(
    shinyswatch.theme.solar(),
    ui.h4("Small molecule visualisation web application"),
    ui.row(
        ui.p(
             """
             This is a Shiny for Python app for viewing and saving 2D images of small molecules of interests. 
             The main backbone libraries used are RDKit and Datamol.
             Two types of PNG images are available: with atom index (atoms labelled with numbers), or without (native skeletal forms).
             2D images of molecules can be saved as PNG files via MolToFile(). 
             Atoms & bonds highlighting are done via MolToImage().
             An option to save each of the four images as a merged version in a grid format
             (1 - top left, 2 - top right, 3 - bottom left & 4 - bottom right).
             A final feature is the compound name search function in the interactive dataframe at the bottom of the app.
             """
        ),
        ui.column(
            4,
            ui.navset_tab_card(
                ui.nav(
                    "Unindexed - 1",
                    ui.input_numeric("mol1", "Select compound via index number:", 0, min = 0, max = 143),
                    ui.input_select("image_style1", "Choose substructure highlights:", 
                                    {"no_high": "Without highlights",
                                     "high": "With highlights"}),
                    ui.input_text_area("atom1", "Enter atom number to highlight:", placeholder = "e.g. 0, 1, 2, 3..."),
                    ui.input_text_area("bond1", "Enter bond number to highlight:", placeholder = "e.g. 0, 1, 2, 3..."),
                    ui.input_text("filename1", "File name for PNG image:"),
                    ui.input_action_button("btn1", "Save and view", class_= "btn-success"),
                    ui.output_image("image1")
                ),
                ui.nav(
                    "Indexed - 2",
                    ui.input_numeric("mol2", "Select compound via index number:", 0, min = 0, max = 143),
                    ui.input_select("image_style2", "Choose substructure highlights:", 
                                    {"no_high": "Without highlights", 
                                     "high": "With highlights"}),
                    ui.input_text_area("atom2", "Enter atom number to highlight:", placeholder = "e.g. 0, 1, 2, 3..."),
                    ui.input_text_area("bond2", "Enter bond number to highlight:", placeholder = "e.g. 0, 1, 2, 3..."),
                    ui.input_text("filename2", "File name for PNG image:"),
                    ui.input_action_button("btn2", "Save and view", class_= "btn-success"),
                    ui.output_image("image2")
                ),
                ui.nav(
                    "Unindexed - 3",
                    ui.input_numeric("mol3", "Select compound via index number:", 0, min = 0, max = 143),
                    ui.input_select("image_style3", "Choose substructure highlights:", 
                                    {"no_high": "Without highlights",
                                     "high": "With highlights"}),
                    ui.input_text_area("atom3", "Enter atom number to highlight:", placeholder = "e.g. 0, 1, 2, 3..."),
                    ui.input_text_area("bond3", "Enter bond number to highlight:", placeholder = "e.g. 0, 1, 2, 3..."),
                    ui.input_text("filename3", "File name for PNG image:"),
                    ui.input_action_button("btn3", "Save and view", class_= "btn-success"),
                    ui.output_image("image3")
                ),
                ui.nav(
                    "Indexed - 4",
                    ui.input_numeric("mol4", "Select compound via index number:", 0, min = 0, max = 143),
                    ui.input_select("image_style4", "Choose substructure highlights:", 
                                    {"no_high": "Without highlights", 
                                     "high": "With highlights"}),
                    ui.input_text_area("atom4", "Enter atom number to highlight:", placeholder = "e.g. 0, 1, 2, 3..."),
                    ui.input_text_area("bond4", "Enter bond number to highlight:", placeholder = "e.g. 0, 1, 2, 3..."),
                    ui.input_text("filename4", "File name for PNG image:"),
                    ui.input_action_button("btn4", "Save and view", class_= "btn-success"),
                    ui.output_image("image4")
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
    # Function to show unindexed compound as PNG image (file 1)
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
    # Function to show indexed compound as PNG image (file 2)
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
    @reactive.Calc
    # Function to show unindexed compound as PNG image (file 3)
    def image3():

        input.btn3()

        with reactive.isolate():

            if input.image_style3() == "high":
                    # Highlight atoms & bonds
                    img = Draw.MolToImage(mols[input.mol3()], 
                                          highlightAtoms = [int(n) for n in input.atom3().split(",")],
                                          highlightBonds = [int(n) for n in input.bond3().split(",")], 
                                          highlightColor = ColorConverter().to_rgb("aqua")
                                          ) 
                    # Save image
                    img.save(f"{input.filename3()}.png")
                    # Show saved PNG file from selected compound
                    dir = Path(__file__).resolve().parent
                    img_3: ImgData = {"src": str(dir / f"{input.filename3()}.png")}
                    return img_3
                
            elif input.image_style3() == "no_high":
                    Draw.MolToFile(mols[input.mol3()], f"{input.filename3()}.png")
                    dir = Path(__file__).resolve().parent
                    img_3: ImgData = {"src": str(dir / f"{input.filename3()}.png")}
                    return img_3
       
            else:
                return None


    @output
    @render.image
    @reactive.Calc
    # Function to show indexed compound as PNG image (file 4)
    def image4():

        input.btn4()

        with reactive.isolate():

            if input.image_style4() == "no_high":
                for atom in mols[input.mol4()].GetAtoms():
                    atom.SetProp('atomLabel', str(atom.GetIdx()))
                Draw.MolToFile(mols[input.mol4()], f"{input.filename4()}.png")
                dir = Path(__file__).resolve().parent
                img_4: ImgData = {"src": str(dir / f"{input.filename4()}.png")}
                return img_4
            
            elif input.image_style4() == "high":
                for atom in mols[input.mol4()].GetAtoms():
                    atom.SetProp('atomLabel', str(atom.GetIdx()))
                # Highlight atoms & bonds
                img = Draw.MolToImage(mols[input.mol4()], 
                                        highlightAtoms = [int(n) for n in input.atom4().split(",")],
                                        highlightBonds = [int(n) for n in input.bond4().split(",")], 
                                        highlightColor = ColorConverter().to_rgb("aqua")
                                        ) 
                # Save image
                img.save(f"{input.filename4()}.png")
                # Show saved PNG file from selected compound
                dir = Path(__file__).resolve().parent
                img_4: ImgData = {"src": str(dir / f"{input.filename4()}.png")}
                return img_4
            
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