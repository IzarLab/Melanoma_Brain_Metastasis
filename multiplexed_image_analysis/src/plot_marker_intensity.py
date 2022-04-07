import matplotlib.pyplot as plt
import seaborn as sns

def plot_marker_intensity(data, marker, cutoff=None, out=None):
    # Subsample data to 10^5 if bigger than 10^5
    data = data.sample(10000, random_state=42) if len(data)>=10000 else data

    # Create canvas
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))

    # Plot
    sns.histplot(data=data, x=marker, bins=200, ax=ax)
    sns.rugplot(data=data, x=marker, height=-.01, clip_on=False, lw=1, alpha=.1, ax=ax)

    # If cutoff provided draw a vertical line
    if cutoff is not None:
        ax.axvline(cutoff, 0, 1, color="r")

    ax.set_title(marker, fontsize=20)
    ax.set_xlabel("Intensity")
    sns.despine()
    plt.tight_layout()

    # Save figure
    if out is not None:
        out = [out] if isinstance(out, str) else out
        for f in out:
            fig.savefig(f, dpi=300)