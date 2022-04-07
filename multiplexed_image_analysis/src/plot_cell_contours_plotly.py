import plotly.express as px
from skimage import io, measure
from sklearn.preprocessing import MinMaxScaler


def plot_cell_contours_plotly(im, mask, data, out, marker, cutoff, channel=0, sx=(0,-1), sy=(0,-1)):
    """
    This function creates an interactive plot of a masked marker together with the cell contours.
    It can color the cell contours based on a given cutoff value.
    """

    # Reset index so the CellID is the index
    data = data.set_index("CellID")

    # Cell
    flag = data[marker].apply(lambda x: marker if x >= cutoff else "")
    color = data[marker].apply(lambda x: "cyan" if x >= cutoff else "red")

    # Select image channel
    im = im[channel]

    # Scale image
    im = MinMaxScaler(feature_range=(0,1), clip=True).fit_transform(im)

    # Sub-setting the image if required
    im = im[sy[0]:sy[1], sx[0]:sx[1]]
    mask = mask[sy[0]:sy[1], sx[0]:sx[1]]

    # Get range of cells
    min = mask[np.nonzero(mask)].min()
    max = mask[np.nonzero(mask)].max()

    # Plot masked marker
    fig = px.imshow(im*mask.astype(bool), color_continuous_scale="Viridis")

    # Plot contours each cell at a time
    for i in range (min, max+1):
        try:
            # Find contours
            y, x = measure.find_contours(mask==i, 0.8)[0].T

            #
            fig.add_trace(go.Scatter(x=x, y=y,
                                     mode='lines',
                                     showlegend=False,
                                     opacity=0.4,
                                     text=flag[i],
                                     line=dict(color=color[i])))
        except IndexError:
            print(f"Can not find contour on mask {i}")

    # Save image
    fig.write_html(out)