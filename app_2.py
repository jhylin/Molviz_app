import pandas as pd
import polars as pl
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.rdmolfiles import SmilesWriter, SmilesMolSupplier
import datamol as dm
from shiny import App, render, ui

app_ui = ui.page_fluid(
    ui.output_table("table"),
)


def server(input, output, session):
    @output
    @render.table
    def table():
        #infile = Path(__file__).parent / "test.csv"
        #df = pd.read_csv("df_ai.csv")
        df = pl.read_csv("df_ai.csv")
        df = df.to_pandas()
        df["mol"] = df["Smiles"].apply(lambda x: dm.to_mol(x))
        return df


app = App(app_ui, server)