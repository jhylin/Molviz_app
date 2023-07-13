# ***Import all libraries or packages needed***
# Import shiny ui, app
from shiny import *
# Import shinywidgets
from shinywidgets import output_widget, render_widget
# Import shinyswatch to add themes
import shinyswatch

# Import pandas
import pandas as pd
import polars as pl

#import ipywidgets as ipy
#from IPython.display import HTML, IFrame, display

# Add RDKit
from rdkit import Chem
from rdkit.Chem import Draw
#from rdkit.Chem.rdmolfiles import SmilesWriter, SmilesMolSupplier
from io import StringIO
import datamol as dm


# ***Specify data source***
df = pl.read_csv("df_ai.csv")
#df.head()
df = df.to_pandas()
#df.head()
#type(df)
df["mol"] = df.Smiles.apply(Chem.MolFromSmiles)
#df.head()


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
                "Name",
                "Choose a compound:", 
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


        
    
        
# Combine UI & server into Shiny app
app = App(app_ui, server)