import logging
import io
from pathlib import Path

import pandas as pd
from PIL import Image
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

# import logging details
import molimg
logger = logging.getLogger(__name__)

def df_column_of_smiles_to_pngs(df: pd.DataFrame, column: str, save_folder: str):
    """Saves all smiles in df column to png format

    Args:
        df: dataframe with smiles column
        column: column in dataframe with smiles strings to save
        save_folder: folder to save images to
    """
    logger.debug(f"Converting smiles in {column} to images")
    save_folder = Path(save_folder)
    save_folder.mkdir(exist_ok=True)
    for index, row in df.iterrows():
        try:
            smiles_to_png(row[column], f"{str(save_folder / row[column])}.png")
        except Exception as e:
            logger.warning(e)
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
    im.save(filename, "PNG")
