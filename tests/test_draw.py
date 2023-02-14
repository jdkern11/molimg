from pathlib import Path

import pandas as pd
import pytest

from molimg import draw


@pytest.fixture
def df():
    folder = Path(__file__).resolve().parent
    return pd.read_csv(str(folder / "data" / "example_data.csv"))


def test_smiles_to_png(tmp_path, df):
    for index, row in df.iterrows():
        if row["smiles1"] == "UrBi" or pd.isna(row["smiles1"]):
            # To get Boost.Python.ArgumentError. This is a C++ error and can only
            # be caught with Excpetion
            with pytest.raises(Exception):
                draw.smiles_to_png(
                    row["smiles1"], str(tmp_path / f"{row['smiles1']}.png")
                )
        else:
            draw.smiles_to_png(row["smiles1"], str(tmp_path / f"{row['smiles1']}.png"))

    for index, row in df.iterrows():
        if row["smiles1"] != "UrBi" and not pd.isna(row["smiles1"]):
            assert (tmp_path / f"{row['smiles1']}.png").exists()
        else:
            assert not (tmp_path / f"{row['smiles1']}.png").exists()


def test_column_of_smiles_to_pngs(tmp_path, df):
    draw.df_column_of_smiles_to_pngs(df, "smiles1", tmp_path)
    for index, row in df.iterrows():
        if row["smiles1"] != "UrBi" and not pd.isna(row["smiles1"]):
            assert (tmp_path / f"{row['smiles1']}.png").exists()
        else:
            assert not (tmp_path / f"{row['smiles1']}.png").exists()


def test_columns_of_smiles_to_pngs(tmp_path, df):
    draw.df_columns_of_smiles_to_pngs(df, ["smiles1", "smiles2"], tmp_path)
    for index, row in df.iterrows():
        for column in ["smiles1", "smiles2"]:
            if row[column] != "UrBi" and not pd.isna(row[column]):
                assert (tmp_path / column / f"{row[column]}.png").exists()
            else:
                assert not (tmp_path / column / f"{row[column]}.png").exists()
