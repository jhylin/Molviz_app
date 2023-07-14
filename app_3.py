import pandas as pd
import polars as pl
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
#from rdkit.Chem.rdmolfiles import SmilesWriter, SmilesMolSupplier
from rdkit.Chem.Draw import rdMolDraw2D
import io
from PIL import Image

import datamol as dm
from shiny import App, render, ui

app_ui = ui.page_fluid(
    ui.output_image("table"),
)


def server(input, output, session):
    @output
    @render.image
    # def table():
    #     #df = pd.read_csv("df_ai.csv")
    #     df = pl.read_csv("df_ai.csv")
    #     df = df.to_pandas()
    #     df["mol"] = df["Smiles"].apply(lambda x: dm.to_mol(x))
    #     mols = df["mol"]
    #     mols = list(mols)
    #     # image = dm.viz.to_image(mols)
    #     Draw.MolsToGridImage(mols)

    #     return table

    def image():
        mol = Chem.MolFromSmiles("CC(CCCC(C)(C)O)C1CCC2C1(CCCC2=CC=C3CC(CC(C3=C)O)O)C")
        Draw.MolsToGridImage(mol)
        return image


app = App(app_ui, server)