from pathlib import Path

import pandas as pd
import pytest

@pytest.fixture
def data():
    return pd.read_csv(str(Path('data') / 'example_data.csv'))


