from pathlib import Path

from PIL import Image
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

def smiles_to_png(smiles: str, filename: str):
    """Creates png of smiles string"""
    mol = Chem.MolFromSmiles(smiles)
	drawer = rdMolDraw2D.MolDraw2DCairo(200, 200)
    drawer.SetFontSize(6.0)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    png = drawer.GetDrawingText()
    im = Image.open(io.BytesIO(png))
    im.save(name, "PNG")
