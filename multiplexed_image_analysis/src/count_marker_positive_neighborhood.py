import numpy as np
import pandas as pd
from pathlib import Path
import datatable as dt
from sklearn.metrics import pairwise_distances
from os.path import basename

def count_marker_positive_neighborhood(path, marker, cutoff, distance=30):
    # Generata a list with all the files on the folder
    files = Path(path).rglob("*.csv")

    # Iterate over all files
    di = dict()
    for f in files:

        # Read in data using datatable
        data = dt.fread(f, columns={marker:marker, "X_centroid":"X", "Y_centroid":"Y", ... : None}).to_pandas()

        # Call marker positive cells by applying the cutoff value
        data = data[data[marker] >= cutoff]

        # If less than 2 cells are present a distance can't be calculated
        if len(data)>1:

            # Calculate pairwise distances between cells
            dists = pairwise_distances(data[["X", "Y"]], metric="euclidean", n_jobs=-1)

            # Get marker positive cells in neighborhood
            neighborhood = np.where(dists <= distance, dists, 0)

            # Count cells  in neighborhood
            neighborhood = np.count_nonzero(neighborhood, axis=0)

            # Get average number of marker positive cells in neighborhood
            neighborhood = neighborhood.mean()

        else:
            neighborhood = np.nan

        # Aggregate the results together
        di[basename(f)] = neighborhood

    # Convert dictionary to Series
    neighborhood = pd.Series(di)

    return neighborhood


