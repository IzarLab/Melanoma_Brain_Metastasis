from os.path import basename
import pandas as pd
import datatable as dt
from pathlib import Path

def read_markers_labeled(path, marker, metadata):
    data = list()

    # Generate a list with all the files on the folder
    files = Path(path).rglob("*.csv")

    # Iterate over the results folder and load marker(s) with labels
    for f in files:
        frame = dt.fread(f, columns={marker:marker, ...:None}).to_pandas()

        # Add labels
        for m in metadata:
            frame[m] = metadata.loc[basename(f)][m]

        data.append(frame)

    # Concatenate dataframes
    data = pd.concat(data)
    data = data.reset_index(drop=True)

    return data