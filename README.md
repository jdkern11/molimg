# Molecular Imager
Do you ever wish you could easily embed the images of your smiles strings into 
an excel sheet? Wish no more! molimg is here to do just that!

Take the following data in a csv:

![image of example data](https://raw.githubusercontent.com/jdkern11/molimg/main/images/example_csv.png)

and molimg will convert it like so:

![image of example data](https://raw.githubusercontent.com/jdkern11/molimg/main/images/example_csv_with_images.png)

## Usage
Assuming I am writing a script in the folder that contains the example data,
then I would write the following code:

```Python
import pandas as pd
from molimg import excel

df = pd.read_csv('example_data.csv')
smiles_columns = ['smiles1', 'smiles2']
excel.write(
    df=df, 
    smiles_columns=smiles_columns, 
    filename='example_data_with_images.xlsx'
)
```

The order the columns appear in df.columns is how the columns will be saved in
the new excel sheet. The new smiles columns with images will always appear to the right
of the data they originate from. 
