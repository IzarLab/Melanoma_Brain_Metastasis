from os.path import abspath
import matplotlib.pyplot as plt
import seaborn as sns

def plot_sample_distributions(data, x, marker, hue, palette, out):
    # Create canvas
    fig, ax = plt.subplots(1, 1, figsize=(15, 10))

    # Plot
    sns.scatterplot(data=data, x=x, y=marker, hue=hue, palette=palette, linewidth=0, ax=ax)
    ax.set_title(marker, fontsize=20)
    ax.tick_params(axis="x", labelrotation=90)
    plt.tight_layout()

    # Save figure
    if out is not None:
        out = [out] if isinstance(out, str) else out
        for f in out:
            fig.savefig(f, dpi=300)