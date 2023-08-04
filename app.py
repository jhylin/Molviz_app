# Shiny for Python app for viewing and saving 2D images of small molecules of interests

# Current version: Improved image resolution by using rdMolDraw2D module, with highlighting option

# old version: 2D images of molecules saved as PNG files via MolToFile() 
# old version: Atoms & bonds highlighting via MolToImage() 
# old version: Select atom & bond highlighting (on or off) with or without atom indices
# old version: Changed to "atomNote" for atom labelling



# Import libraries---
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
import datamol as dm
from itables.shiny import DT
from pathlib import Path
from PIL import Image
import cairosvg
#from matplotlib.colors import ColorConverter 
from shiny import App, Inputs, Outputs, Session, render, ui, reactive
from shiny.types import ImgData
import shinyswatch


# Data source---
# Using cleaned, pre-processed and standardised data 
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
    shinyswatch.theme.cyborg(),
    {"class": "col-lg-20 py-4 mx-auto text-left"},
    ui.div(
        {"style": "font-weight: bold;"},
        ui.h4("Molecule visualiser - Molviz"),
    ),
    ui.row(
        ui.column(
            12,
            ui.div(
                {"class": "app-col"},
                ui.p(
                    """
                    This is an application built by using Shiny for Python web application framework. It provides features for viewing and saving 2D images of small molecules of interests. 
                    The main libraries used are RDKit and Datamol.
                    The data source is based on compound data in a dataframe, which included molecular representations such as SMILES.
                    For demonstration purpose, the example provided here is a set of anti-infectives sourced from an older version of ChEMBL.
                    Users can generate similar app by using tailored compound data for structure-activity relationships or other drug discovery-related work.
                    For code and license:
                    """,
                ui.a("here", href="https://github.com/jhylin/Molviz_app"),
                    """
                    or visit https://github.com/jhylin/Molviz_app.
                    """
                ),
            ),
        )
    ),
    ui.row(
        ui.column(
            5,
            ui.navset_tab_card(
                ui.nav(
                    "1",
                    ui.input_numeric("mol1", "Select compound via index number from dataframe:", 0, min = 0, max = 143),
                    ui.input_select("image_style1", "With or without substructure highlights:", 
                                    {"no_high": "Without highlights",
                                     "high": "With highlights"}),
                    ui.input_text_area("atom1", "Enter atom numbers to highlight:", placeholder = "e.g. 0, 1, 2, 3..."),
                    ui.input_text_area("bond1", "Enter bond numbers to highlight:", placeholder = "e.g. 0, 1, 2, 3..."),
                    ui.input_text("filename1", "File name for PNG image:"),
                    ui.input_action_button("btn1", "Save and view", class_= "btn-success"),
                    ui.output_image("image1")
                ),
                ui.nav(
                    "2",
                    ui.input_numeric("mol2", "Select compound via index number from dataframe:", 0, min = 0, max = 143),
                    ui.input_select("image_style2", "With or without substructure highlights:", 
                                    {"no_high": "Without highlights", 
                                     "high": "With highlights"}),
                    ui.input_text_area("atom2", "Enter atom numbers to highlight:", placeholder = "e.g. 0, 1, 2, 3..."),
                    ui.input_text_area("bond2", "Enter bond numbers to highlight:", placeholder = "e.g. 0, 1, 2, 3..."),
                    ui.input_text("filename2", "File name for PNG image:"),
                    ui.input_action_button("btn2", "Save and view", class_= "btn-success"),
                    ui.output_image("image2")
                ),
                ui.nav(
                    "3",
                    ui.input_numeric("mol3", "Select compound via index number from dataframe:", 0, min = 0, max = 143),
                    ui.input_select("image_style3", "With or without substructure highlights:", 
                                    {"no_high": "Without highlights",
                                     "high": "With highlights"}),
                    ui.input_text_area("atom3", "Enter atom numbers to highlight:", placeholder = "e.g. 0, 1, 2, 3..."),
                    ui.input_text_area("bond3", "Enter bond numbers to highlight:", placeholder = "e.g. 0, 1, 2, 3..."),
                    ui.input_text("filename3", "File name for PNG image:"),
                    ui.input_action_button("btn3", "Save and view", class_= "btn-success"),
                    ui.output_image("image3")
                ),
                ui.nav(
                    "4",
                    ui.input_numeric("mol4", "Select compound via index number from dataframe:", 0, min = 0, max = 143),
                    ui.input_select("image_style4", "With or without substructure highlights:", 
                                    {"no_high": "Without highlights", 
                                     "high": "With highlights"}),
                    ui.input_text_area("atom4", "Enter atom numbers to highlight:", placeholder = "e.g. 0, 1, 2, 3..."),
                    ui.input_text_area("bond4", "Enter bond numbers to highlight:", placeholder = "e.g. 0, 1, 2, 3..."),
                    ui.input_text("filename4", "File name for PNG image:"),
                    ui.input_action_button("btn4", "Save and view", class_= "btn-success"),
                    ui.output_image("image4")
                ),
            )
        ),
        ui.column(
            6, 
            ui.tags.ul(
            {"class": "col-lg-20 py-2 mx-auto text-left"},
            ui.div(
                ui.tags.li(
                """
                Save each compound individually as PNG files, with atom indices shown"""
            ), 
                ui.tags.li(
                """
                Save merged image, which contains 4 images in this order: 
                top left (1), right (2), bottom left (3) and right (4) once file name is entered and confirmed
                - numbering corresponds to 4 tabs on the left"""
            ), 

            )

            ),
            ui.input_text("merge_filename", "File name for merged PNG images:"),
            ui.input_action_button("btn_merge", "Confirm", class_="btn-success"),
            ui.output_image("merge_image"),
            ),
            ui.div(
                {"class": "col-lg-15 py-4 mx-auto text-left"},
                {"style": "font-weight: bold;"},
                ui.h5("Interactive data table"),
            ), 
            ui.page_fluid(
                ui.div(
                {"class": "col-lg-15 py-1 mx-auto text-left"},
                ui.tags.ul(
                ui.tags.li(
                    """
                    Search by compound name from the data table below
                    """
                ),
                ui.tags.li(
                    """
                    Browse each page by clicking on the page number buttons at the bottom
                    """
                ),
                ui.tags.li(
                    """
                    Note: searching by index number may not work as well as by name (tend to include mixed results)
                    """
                ))), 
                ui.HTML(DT(df))
                ),
            ),
)


# Server output--- 
def server(input, output, session):

    @output
    @render.image
    @reactive.Calc
    # Added reactive.event() to avoid triggering initial image when 1st opening the app 
    # reaction will only occur after hitting action button
    @reactive.event(input.btn1)
    # Function to save & show 1st compound as PNG image (file 1)
    def image1():

        input.btn1()

        with reactive.isolate():

            if input.image_style1() == "high":
                d = rdMolDraw2D.MolDraw2DSVG(300, 300)
                d.drawOptions().addAtomIndices = True
                # Potentially fills up white spaces too much if using bond indices
                #d.drawOptions().addBondIndices = False
                d.DrawMolecule(mols[input.mol1()], 
                               highlightAtoms = [int(n) for n in input.atom1().split(",")], 
                               highlightBonds = [int(n) for n in input.bond1().split(",")]
                               )
                d.FinishDrawing()
                cairosvg.svg2png(bytestring = d.GetDrawingText().encode(), 
                                 # dpi = dots per inch
                                 dpi = 50000,
                                 scale = 10, 
                                 write_to = f"{input.filename1()}.png",
                                 output_height = 400,
                                 output_width = 400
                                 )
                dir = Path(__file__).resolve().parent
                img_1: ImgData = {"src": str(dir / f"{input.filename1()}.png")}
                return img_1
                
            elif input.image_style1() == "no_high":
                d = rdMolDraw2D.MolDraw2DSVG(300, 300)
                d.drawOptions().addAtomIndices = True
                d.DrawMolecule(mols[input.mol1()])
                d.FinishDrawing()
                cairosvg.svg2png(bytestring = d.GetDrawingText().encode(), 
                                 # dpi = dots per inch
                                 dpi = 50000,
                                 scale = 10, 
                                 write_to = f"{input.filename1()}.png",
                                 output_height = 400,
                                 output_width = 400
                                 )

                dir = Path(__file__).resolve().parent
                img_1: ImgData = {"src": str(dir / f"{input.filename1()}.png")}
                return img_1
       
            else:
                return None
    


    @output
    @render.image
    @reactive.Calc
    @reactive.event(input.btn2)
    # Function to save & show 2nd compound as PNG image (file 2)
    def image2():

        input.btn2()

        with reactive.isolate():

            if input.image_style2() == "no_high":
                d = rdMolDraw2D.MolDraw2DSVG(300, 300)
                d.drawOptions().addAtomIndices = True
                d.DrawMolecule(mols[input.mol2()])
                d.FinishDrawing()
                cairosvg.svg2png(bytestring = d.GetDrawingText().encode(), 
                                 dpi = 30000,
                                 scale = 5, 
                                 write_to = f"{input.filename2()}.png",
                                 output_height = 400,
                                 output_width = 400
                                 )
                dir = Path(__file__).resolve().parent
                img_2: ImgData = {"src": str(dir / f"{input.filename2()}.png")}
                return img_2
            
            elif input.image_style2() == "high":
                d = rdMolDraw2D.MolDraw2DSVG(300, 300)
                d.drawOptions().addAtomIndices = True
                #d.drawOptions().addBondIndices = False
                d.DrawMolecule(mols[input.mol2()], 
                               highlightAtoms = [int(n) for n in input.atom2().split(",")], 
                               highlightBonds = [int(n) for n in input.bond2().split(",")]
                               )
                d.FinishDrawing()
                cairosvg.svg2png(bytestring = d.GetDrawingText().encode(), 
                                 dpi = 30000,
                                 scale = 5, 
                                 write_to = f"{input.filename2()}.png",
                                 output_height = 400,
                                 output_width = 400
                                 )
                dir = Path(__file__).resolve().parent
                img_2: ImgData = {"src": str(dir / f"{input.filename2()}.png")}
                return img_2
            
            else:
                return None
            
    @output
    @render.image
    @reactive.Calc
    @reactive.event(input.btn3)
    # Function to save & show 3rd compound as PNG image (file 3)
    def image3():

        input.btn3()

        with reactive.isolate():

            if input.image_style3() == "high":
                d = rdMolDraw2D.MolDraw2DSVG(300, 300)
                d.drawOptions().addAtomIndices = True
                #d.drawOptions().addBondIndices = False
                d.DrawMolecule(mols[input.mol3()], 
                               highlightAtoms = [int(n) for n in input.atom3().split(",")], 
                               highlightBonds = [int(n) for n in input.bond3().split(",")]
                               )
                d.FinishDrawing()
                cairosvg.svg2png(bytestring = d.GetDrawingText().encode(), 
                                 dpi = 30000,
                                 scale = 5, 
                                 write_to = f"{input.filename3()}.png",
                                 output_height = 400,
                                 output_width = 400
                                 )
                dir = Path(__file__).resolve().parent
                img_3: ImgData = {"src": str(dir / f"{input.filename3()}.png")}
                return img_3

            elif input.image_style3() == "no_high":
                d = rdMolDraw2D.MolDraw2DSVG(300, 300)
                d.drawOptions().addAtomIndices = True
                d.DrawMolecule(mols[input.mol3()])
                d.FinishDrawing()
                cairosvg.svg2png(bytestring = d.GetDrawingText().encode(), 
                                 dpi = 30000,
                                 scale = 5, 
                                 write_to = f"{input.filename3()}.png",
                                 output_height = 400,
                                 output_width = 400
                                 )
                dir = Path(__file__).resolve().parent
                img_3: ImgData = {"src": str(dir / f"{input.filename3()}.png")}
                return img_3

            else:
                return None


    @output
    @render.image
    @reactive.Calc
    @reactive.event(input.btn4)
    # Function to save & show 4th compound as PNG image (file 4)
    def image4():

        input.btn4()

        with reactive.isolate():

            if input.image_style4() == "no_high":
                d = rdMolDraw2D.MolDraw2DSVG(300, 300)
                d.drawOptions().addAtomIndices = True
                d.DrawMolecule(mols[input.mol4()])
                d.FinishDrawing()
                cairosvg.svg2png(bytestring = d.GetDrawingText().encode(), 
                                 dpi = 30000,
                                 scale = 5, 
                                 write_to = f"{input.filename4()}.png",
                                 output_height = 400,
                                 output_width = 400
                                 )
                dir = Path(__file__).resolve().parent
                img_4: ImgData = {"src": str(dir / f"{input.filename4()}.png")}
                return img_4

            
            elif input.image_style4() == "high":
                d = rdMolDraw2D.MolDraw2DSVG(300, 300)
                d.drawOptions().addAtomIndices = True
                #d.drawOptions().addBondIndices = False
                d.DrawMolecule(mols[input.mol4()], 
                               highlightAtoms = [int(n) for n in input.atom4().split(",")], 
                               highlightBonds = [int(n) for n in input.bond4().split(",")]
                               )
                d.FinishDrawing()
                cairosvg.svg2png(bytestring = d.GetDrawingText().encode(), 
                                 dpi = 30000,
                                 scale = 5, 
                                 write_to = f"{input.filename4()}.png",
                                 output_height = 400,
                                 output_width = 400
                                 )
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
            blank_image = Image.new("RGB", (800, 800))

            # Paste img1 to img4 together in a grid (top left, right, bottom left, right)
            blank_image.paste(img1, (0, 0))
            blank_image.paste(img2, (400, 0))
            blank_image.paste(img3, (0, 400))
            blank_image.paste(img4, (400, 400))

            # Save combined/merged image as a new PNG file
            blank_image.save(f"{input.merge_filename()}.png")
        
        # Show merged image of selected PNG images
        dir = Path(__file__).resolve().parent
        img_merge: ImgData = {"src": str(dir / f"{input.merge_filename()}.png")}
        return img_merge
    

app = App(app_ui, server)