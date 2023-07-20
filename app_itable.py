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
df = pd.read_csv("df_ai.csv")
# Keeping index column
df = df.reset_index()


# Input---
app_ui = ui.page_fluid(
    ui.page_fluid(ui.HTML(DT(df)))
)


# Server not required so set it to none
app = App(app_ui, server = None)