# App reading and showing dataframe as an interactive table
# for viewing and searching data by using keywords

# Import libraries---
import pandas as pd
import polars as pl
from rdkit import Chem
from rdkit.Chem import Draw
import datamol as dm
from shiny import App, render, ui
from itables.shiny import DT


# Data source---
df = pl.read_csv("df_ai.csv")
df = df.to_pandas()
#df["mol"] = df["Smiles"].apply(lambda x: dm.to_mol(x))


# Input---
app_ui = ui.page_fluid(
    #ui.output_table("table"),
    ui.page_fluid(ui.HTML(DT(df)))
)


# Server not required so set it to none
app = App(app_ui, server = None)