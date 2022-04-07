import matplotlib.pyplot as plt
import seaborn as sns


def plot_distributions(data, marker, hue, palette, out=None):
    # Create canvas
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))

    # Plot
    sns.violinplot(data=data, x=hue, y=marker, hue=hue, palette=palette, ax=ax)
    ax.set_title(f"{marker} Intensity distribution", fontsize=20)
    sns.despine(left=True)
    plt.tight_layout()

    # Save figure
    if out is not None:
        out = [out] if isinstance(out, str) else out
        for f in out:
            fig.savefig(f, dpi=300)


def plot_distribution_boxen(data, marker, hue, palette, out=None):
    # Create canvas
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))

    # Plot
    sns.boxenplot(data=data, x=hue, y=marker, hue=hue, palette=palette, ax=ax)
    ax.set_title(marker, fontsize=20)
    sns.despine(left=True)
    plt.tight_layout()

    # Save figure
    if out is not None:
        out = [out] if isinstance(out, str) else out
        for f in out:
            fig.savefig(f, dpi=300)