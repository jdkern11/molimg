from pathlib import Path


def remove_extensions(filename: str):
    return str(Path(filename).with_suffix(""))
