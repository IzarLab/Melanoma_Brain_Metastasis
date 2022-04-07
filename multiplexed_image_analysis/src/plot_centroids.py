import datatable as dt
import pandas as pd
import plotly.express as px

def plot_centroids(path, marker, cutoff=0, show=False, out=False):
    # Read in data
    data = dt.fread(path, columns={"CellID":"CellID", marker:marker, "X_centroid":"X", "Y_centroid":"Y", ...:None}).to_pandas()

    # Get positive cells for marker based on cutoff
    data[f"{marker}_positive"] = data[marker] >= cutoff

    # Make a plotly plot
    fig = px.scatter(data, x="X", y="Y", color=f"{marker}_positive", width=1400, height=1000, text=marker, opacity=0.5)
    fig.update_layout(yaxis=dict(autorange="reversed"), template="plotly_dark")

    # Show figure
    if show:
        fig.show()

    # Save figure
    if out:
        fig.write_html(out)