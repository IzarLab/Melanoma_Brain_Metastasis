from os.path import basename
import datatable as dt
from pathlib import Path
import pandas as pd


def get_marker_positive_cells(path, marker, cutoff):
    # Generate a list with all the files on the folder
    files = Path(path).rglob("*.csv")

    dict = {}
    # Iterate over each file
    for f in files:

        # Read in data
        data = dt.fread(f, columns={marker:marker, ...:None}).to_pandas()

        # Get positive cells
        dict[basename(f)] = len(data[data[marker]>=cutoff])

    # Create a dataframe and output
    pos = pd.Series(dict)

    return pos