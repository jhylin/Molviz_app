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

app_ui = ui.page_fluid(
    ui.h4("Compound input selections"),
    ui.row(
        ui.column(
            4,
            ui.navset_tab_card(
                ui.nav(
                    "1",
                    ui.input_numeric("mol1", "Select 1st compound via index number:", 0, min = 0, max = 143),
                    ui.input_select("image_style", "Choices for compound image:", 
                                    {"no_high_no_idx": "No highlights and no index",
                                     "high_no_idx": "Highlights without index",
                                     "no_high_idx": "No highlights with index", 
                                     "high_idx": "Highlights with index"}),
                    ui.input_text_area("atom1", "Enter atom number to highlight", placeholder = "e.g. 0, 1, 2, 3..."),
                    ui.input_text_area("bond1", "Enter bond number to highlight:", placeholder = "e.g. 0, 1, 2, 3..."),
                    #ui.input_checkbox("on", "Show index"),
                    ui.input_text("filename1_1", "File name for compound image 1:"),
                    ui.input_text("filename1_2", "File name for compound image 2:"),
                    # Add caption explaining everytime confirm button is clicked, 
                    # a new PNG file will be saved with whatever file name was inputted above
                    # Make desired selections before clicking on the confirm button with diff file names for each
                    ui.input_action_button("btn1", "Save", class_= "btn"),
                    ui.output_image("image1")
                ),
                ui.nav(
                    "2",
                    ui.input_numeric("mol2", "Select 2nd compound via index number:", 0, min = 0, max = 143),
                    ui.input_text_area("atom2", "Enter atom number to highlight", placeholder = "e.g. 0, 1, 2, 3..."),
                    ui.input_text_area("bond2", "Enter bond number to highlight:", placeholder = "e.g. 0, 1, 2, 3..."),
                    ui.input_select("image_style2", "Choices for compound image:", {"no_highlight": "Without highlights", "highlight": "With highlights"}),
                    ui.input_checkbox("on2", "Show index"),
                    ui.input_text("filename2", "File name for compound:"),
                    ui.input_action_button("btn2", "Confirm", class_="btn"),
                    ui.output_image("image2")
                ),
                ui.nav(
                    "3",
                    ui.input_numeric("mol3", "Select 3rd compound via index number:", 0, min = 0, max = 143),
                    ui.input_text_area("atom3", "Enter atom number to highlight", placeholder = "e.g. 0, 1, 2, 3..."),
                    ui.input_text_area("bond3", "Enter bond number to highlight:", placeholder = "e.g. 0, 1, 2, 3..."),
                    ui.input_select("image_style3", "Choices for compound image:", {"no_highlight": "Without highlights", "highlight": "With highlights"}),
                    ui.input_checkbox("on3", "Show index"),
                    ui.input_text("filename3", "File name for compound:"),
                    ui.input_action_button("btn3", "Confirm", class_="btn"),
                    ui.output_image("image3")
                ),
                ui.nav(
                    "4",
                    ui.input_numeric("mol4", "Select 4th compound via index number:", 0, min = 0, max = 143),
                    ui.input_text_area("atom4", "Enter atom number to highlight", placeholder = "e.g. 0, 1, 2, 3..."),
                    ui.input_text_area("bond4", "Enter bond number to highlight:", placeholder = "e.g. 0, 1, 2, 3..."),
                    ui.input_select("image_style4", "Choices for compound image:", {"no_highlight": "Without highlights", "highlight": "With highlights"}),
                    ui.input_checkbox("on4", "Show index"),
                    ui.input_text("filename4", "File name for compound:"),
                    ui.input_action_button("btn4", "Confirm", class_="btn"),
                    ui.output_image("image4")
                )
            )
        ),
        ui.column(
            3, 
            ui.input_text("merge_filename", "File name for merged images:"),
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

    # One-off effect for each image option - top down only, for saving PNG image only 
    # (viewing available after saving process)
    # To flip back from index to no index is currently not possible to view in app 
    # (refer to saved PNG images instead)
    @output
    @render.image
    @reactive.Calc
    #@reactive.event(input.btn1)
    # Function to show 1st compound as PNG image
    def image1():

        input.btn1()

        with reactive.isolate():

            #for mol in mols:

                if input.image_style() == "high_no_idx":
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
                
                elif input.image_style() == "no_high_no_idx":
                    Draw.MolToFile(mols[input.mol1()], f"{input.filename1_1()}.png")
                    dir = Path(__file__).resolve().parent
                    img_1: ImgData = {"src": str(dir / f"{input.filename1_1()}.png")}
                    return img_1
                
                elif input.image_style() == "no_high_idx":
                    for atom in mols[input.mol1()].GetAtoms():
                        atom.SetProp('atomLabel', str(atom.GetIdx()))
                    Draw.MolToFile(mols[input.mol1()], f"{input.filename1_2()}.png")
                    dir = Path(__file__).resolve().parent
                    img_1: ImgData = {"src": str(dir / f"{input.filename1_2()}.png")}
                    return img_1
                
                elif input.image_style() == "high_idx":
                    for atom in mols[input.mol1()].GetAtoms():
                        atom.SetProp('atomLabel', str(atom.GetIdx()))
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
            
                else:
                    return None


    @output
    @render.image
    @reactive.Calc
    # Function to show 2nd compound as PNG image
    def image2():

        input.btn2()

        with reactive.isolate():

            if input.on2() and input.image_style2() == "no_highlight":
                for atom in mols[input.mol2()].GetAtoms():
                    atom.SetProp('atomLabel', str(atom.GetIdx()))
                Draw.MolToFile(mols[input.mol2()], f"{input.filename2()}.png")
                dir = Path(__file__).resolve().parent
                img_2: ImgData = {"src": str(dir / f"{input.filename2()}.png")}
                return img_2
            
            elif input.on2() and input.image_style2() == "highlight":
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
            
            elif input.image_style2() == "no_highlight":
                Draw.MolToFile(mols[input.mol2()], f"{input.filename2()}.png")
                dir = Path(__file__).resolve().parent
                img_2: ImgData = {"src": str(dir / f"{input.filename2()}.png")}
                return img_2
           
            elif input.image_style2() == "highlight":
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
    # Function to show 3rd compound as PNG image
    def image3():

        # 3rd action button
        input.btn3()

        with reactive.isolate():

            if input.on3() and input.image_style3() == "no_highlight":
                for atom in mols[input.mol3()].GetAtoms():
                    atom.SetProp('atomLabel', str(atom.GetIdx()))
                Draw.MolToFile(mols[input.mol3()], f"{input.filename3()}.png")
                dir = Path(__file__).resolve().parent
                img_3: ImgData = {"src": str(dir / f"{input.filename3()}.png")}
                return img_3
            
            elif input.on3() and input.image_style3() == "highlight":
                for atom in mols[input.mol3()].GetAtoms():
                    atom.SetProp('atomLabel', str(atom.GetIdx()))
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
            
            elif input.image_style3() == "no_highlight":
                Draw.MolToFile(mols[input.mol3()], f"{input.filename3()}.png")
                dir = Path(__file__).resolve().parent
                img_3: ImgData = {"src": str(dir / f"{input.filename3()}.png")}
                return img_3
           
            elif input.image_style3() == "highlight":
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
            
            else:
                return None
    
    @output
    @render.image
    @reactive.Calc
    # Function to show 4th compound as PNG image
    def image4():

        # 4th action button
        input.btn4()

        with reactive.isolate():

            if input.on4() and input.image_style4() == "no_highlight":
                for atom in mols[input.mol4()].GetAtoms():
                    atom.SetProp('atomLabel', str(atom.GetIdx()))
                Draw.MolToFile(mols[input.mol4()], f"{input.filename4()}.png")
                dir = Path(__file__).resolve().parent
                img_4: ImgData = {"src": str(dir / f"{input.filename4()}.png")}
                return img_4
            
            elif input.on4() and input.image_style4() == "highlight":
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
            
            elif input.image_style4() == "no_highlight":
                Draw.MolToFile(mols[input.mol4()], f"{input.filename4()}.png")
                dir = Path(__file__).resolve().parent
                img_4: ImgData = {"src": str(dir / f"{input.filename4()}.png")}
                return img_4
           
            elif input.image_style4() == "highlight":
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