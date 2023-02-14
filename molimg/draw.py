import logging
import io
from pathlib import Path

import pandas as pd
from PIL import Image
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from molimg.helpers import remove_extensions

logger = logging.getLogger(__name__)


def df_columns_of_smiles_to_pngs(
    df: pd.DataFrame, columns: list[str], save_folder: str
):
    """Saves smiles in multiple columns to png formats.

    Args:
        df: dataframe with smiles columns
        columns: list of column names with smiles strings
        save_folder: Folder to save column folders in. Column folders will be the
            same name as the column string.
    """
    # check prior to running to avoid clean up
    for column in columns:
        if column not in df.columns:
            raise KeyError(f"{column} not in the dataframe")

    save_folder = Path(remove_extensions(save_folder))
    save_folder.mkdir(exist_ok=True)
    for column in columns:
        df_column_of_smiles_to_pngs(df, column, str(save_folder / column))


def df_column_of_smiles_to_pngs(df: pd.DataFrame, column: str, save_folder: str):
    """Saves all smiles in df column to png format

    Args:
        df: dataframe with smiles column
        column: column in dataframe with smiles strings to save
        save_folder: folder to save images to
    """
    logger.debug(f"Converting smiles in {column} to images")
    if column not in df.columns:
        raise KeyError(f"{column} not in the dataframe")
    save_folder = Path(remove_extensions(save_folder))
    save_folder.mkdir(exist_ok=True)
    for index, row in df.iterrows():
        try:
            if not pd.isna(row[column]):
                smiles_to_png(row[column], f"{str(save_folder / row[column])}.png")
        # To get Boost.Python.ArgumentError. This is a C++ error and can only
        # be caught with Excpetion
        except Exception as e:
            logger.warning(f"For {row[column]}, {e}")
            continue


def smiles_to_png(smiles: str, filename: str, width: int = 200, height: int = 200):
    """Creates png of smiles string"""
    logger.debug(f"Converting {smiles} to image")
    mol = Chem.MolFromSmiles(smiles)
    drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
    drawer.SetFontSize(6.0)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    png = drawer.GetDrawingText()
    im = Image.open(io.BytesIO(png))
    im.save(remove_extensions(filename) + ".png", "PNG")
