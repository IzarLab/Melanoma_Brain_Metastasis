from os.path import basename
import datatable as dt
from pathlib import Path
import pandas as pd

def get_marker_intensities(path, marker, cutoff):
    # Generata a list with all the files on the folder
    files = Path(path).rglob("*.csv")

    dict = {}
    # Iterate over each file
    for f in files:

        # Read in data
        data = dt.fread(f, columns={marker:marker, ...:None}).to_pandas()

        # Get average intensity
        dict[basename(f)] = data[data[marker] >= cutoff][marker].mean()

    # Create an output Series
    out  = pd.Series(dict)

    return out