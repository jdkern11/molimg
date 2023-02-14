from pathlib import Path
import pandas as pd
import pytest

from molimg import excel


@pytest.fixture
def df():
    folder = Path(__file__).resolve().parent
    return pd.read_csv(str(folder / "data" / "example_data.csv"))


def test_write(tmp_path, df):
    columns = [col for col in df.columns if "smiles" in col]
    filename = "testing.xlsx"
    excel.write(df, columns, str(tmp_path / filename))
    assert (tmp_path / filename).exists()
