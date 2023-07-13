# ***Import all libraries or packages needed***
# Import shiny ui, app
from shiny import *
# Import shinywidgets
from shinywidgets import output_widget, render_widget
# Import shinyswatch to add themes
import shinyswatch

# Import pandas
import pandas as pd

#import ipywidgets as ipy
#from IPython.display import HTML, IFrame, display

# Add RDKit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.rdmolfiles import SmilesWriter, SmilesMolSupplier
from io import StringIO
import datamol as dm


# ***Specify data source***



# User interface---
# Add inputs & outputs
app_ui = ui.page_fluid(
        # Add theme
        shinyswatch.theme.superhero(),
        # Add heading
        ui.h3("Headings"),
        # Input
        ui.row(
            # Specify input
            ui.input_select(
                "filename",
                "Choose a file:", 
                {"cefe.smi": "cefepime", 
                 "cpd3.smi": "cpd3"},
            ),

            #ui.input_action_button("go", "Go!", class_="btn-success"),
            ),
        # Output
        ui.output_ui("image")
        #output_widget("image")
    )


# Server---
# Add code within my_widget function within the server function
def server(input, output, session):
    @output
    @render_widget
    @render.ui
    #@reactive.event(input.go, ignore_none=False)

    def image(filename):
        # Read in a simple list of SMILES from .smi file
        suppl = SmilesMolSupplier(filename)
        # Convert a list of molecules into a dataframe
        mols = dm.to_df(suppl)
        # Generate a RDKit molecule column
        mols["mol"] = mols.smiles.apply(Chem.MolFromSmiles)
        # Display first molecule in dataframe
        image = mols.iloc[0]["mol"]

    return image
        
    
        
# Combine UI & server into Shiny app
app = App(app_ui, server)