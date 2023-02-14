import shutil
import uuid
import logging
from pathlib import Path

import pandas as pd
import xlsxwriter

from molimg.helpers import remove_extensions
from molimg.draw import df_columns_of_smiles_to_pngs

logger = logging.getLogger(__name__)


def write(df: pd.DataFrame, smiles_columns: list[str], filename: str):
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
            Name of the saved excel files. Extensions will be removed
            (filename.extension) and filename will be saved as "filename.xlsx".
    """
    logger.debug(f"Converting smiles in {smiles_columns} to images")
    filename = remove_extensions(filename)
    folder_name = str(uuid.uuid4())
    # make sure folder doesn't already exist
    while Path(folder_name).exists():
        folder_name = str(uuid.uuid4())

    df_columns_of_smiles_to_pngs(df, smiles_columns, folder_name)

    workbook = xlsxwriter.Workbook(f"{filename}.xlsx")
    worksheet = workbook.add_worksheet()

    r, c = 0, 0
    for column in df.columns:
        worksheet.write(r, c, column)
        for index, row in df.iterrows():
            r += 1
            if pd.isna(row[column]):
                worksheet.write_blank(r, c, None)
            else:
                worksheet.write(r, c, row[column])

        # image column should come next if column is in smiles column
        if column in smiles_columns:
            r = 0
            c += 1
            worksheet.write(r, c, f"{column}_image")
            worksheet.set_column_pixels(c, c, 205)
            for index, row in df.iterrows():
                r += 1
                filepath = Path(folder_name) / column / f"{row[column]}.png"
                if filepath.exists():
                    worksheet.insert_image(
                        r, c, str(filepath), {"x_offset": 2.5, "y_offset": 2.5}
                    )
                    worksheet.set_row_pixels(r, 205)
        c += 1
        r = 0

    workbook.close()
    shutil.rmtree(folder_name)
