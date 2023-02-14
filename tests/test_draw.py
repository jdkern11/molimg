from pathlib import Path

import pandas as pd
import pytest

from molimg import draw

@pytest.fixture
def df():
    return pd.read_csv(str(Path('data') / 'example_data.csv'))

def test_smiles_to_png(tmp_path, df):
    for index, row in df.iterrows():
        if row['smiles1'] == 'UrBi':
            # To get Boost.Python.ArgumentError. This is a C++ error and can only
            # be caught with Excpetion
            with pytest.raises(Exception):
                draw.smiles_to_png(row['smiles1'], str(tmp_path / f"{row['smiles1']}.png"))
        else:
            draw.smiles_to_png(row['smiles1'], str(tmp_path / f"{row['smiles1']}.png"))

    for index, row in df.iterrows():
        if row['smiles1'] != 'UrBi':
            assert (tmp_path / f"{row['smiles1']}.png").exists()


def test_column_of_smiles_to_pngs(tmp_path, df):
    draw.df_column_of_smiles_to_pngs(df, 'smiles1', tmp_path)
    for index, row in df.iterrows():
        if row['smiles1'] != 'UrBi':
            assert (tmp_path / f"{row['smiles1']}.png").exists()
        else:
            assert not (tmp_path / f"{row['smiles1']}.png").exists()

