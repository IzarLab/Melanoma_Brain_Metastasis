import numpy as np
import datatable as dt
import pandas as pd
from sklearn.metrics import pairwise_distances
from pathlib import Path
from os.path import basename


def measure_marker_positive_mean_distance(path, marker, cutoff):
    # Generata a list with all the files on the folder
    files = Path(path).rglob("*.csv")

    # Iterate over all the files
    di = dict()
    for f in files:

        # Read in data using datatable
        data = dt.fread(f, columns={marker:marker, "X_centroid":"X", "Y_centroid":"Y", ... : None}).to_pandas()

        # Call marker positive cells by applying the cutoff value
        data = data[data[marker]>=cutoff]

        # If less than 2 cells are present a distance can't be calculated
        if len(data)>1:
            # Calculate pairwise distances between cells
            dists = pairwise_distances(data[["X", "Y"]], metric="euclidean", n_jobs=-1)

            # Mean distance of marker positive cells
            dists = np.mean(dists)/2
        else:
            dists = np.nan

        di[basename(f)] = dists

    # Convert to Series
    dists = pd.Series(di)

    return dists




