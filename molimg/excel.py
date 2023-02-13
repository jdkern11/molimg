from typing import List

import pandas as pd
import xlsxwriter

def write(df: pd.DataFrame, smiles_columns: List[str], filename: str):
    """Writes dataframe to excel file after converting smiles in smiles_columns to images

    Images are embedded in separate columns in the excel file.

    Args:
        df: 
            Dataframe to save data in
        smiles_columns: 
            List of columns in the dataframe that contain the smiles
            strings. These smiles will be converted to images and the images will
            be embedded in the excel file in a new column.
        filename: 
            Name of the saved excel files.
    """


