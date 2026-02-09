import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def palette_plot(
    colors,
    ax=None,
    order='sort',
    width=0.25,
    height=0,
    vertical=False,
    legend_on_right=True,
    display_ticks=False,
):

    color_series = pd.Series(colors)

    if order == 'sort':
        order = sorted(color_series.index)

    if vertical:
        data = pd.Series(1, index=order)
        if ax is None:
            height = height or 0.3 * len(colors)
            _, ax = plt.subplots(figsize=(height, width))
        data.plot(
            kind='bar', color=[color_series[label] for label in data.index], width=1, ax=ax
        )
        ax.set_yticks([])
    else:
        data = pd.Series(1, index=order[::-1])
        if ax is None:
            height = height or 0.3 * len(colors)
            _, ax = plt.subplots(figsize=(width, height))
        data.plot(
            kind='barh', color=[color_series[label] for label in data.index], width=1, ax=ax
        )
        ax.set_xticks([])
        if legend_on_right:
            ax.yaxis.tick_right()

    sns.despine(offset={'left': -2}, ax=ax)
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)

    if not display_ticks:
        ax.tick_params(length=0)

    return ax

def plot_color_annotations(
    color_series,
    ax=None,  
    missing_color='#ffffff',
    shift=0,
    remove_ticks=True,
    remove_borders=True,
):

    if ax is None:
        _, ax = plt.subplots(figsize=(max(len(color_series) / 15.0, 6), 0.5))

    num_items = len(color_series)

    x_positions = np.arange(num_items) - shift
    y_positions = pd.Series([1] * num_items, index=color_series.index)

    with sns.axes_style("white"):
        ax.bar(
            x_positions,
            y_positions,
            color=color_series.fillna(missing_color),
            width=1,
            align='edge',
            edgecolor=color_series.fillna(missing_color),
        )

    ax.set_ylim(0, 1)
    ax.set_xlim(0, num_items)

    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.xaxis.label.set_visible(False)
    ax.set_ylabel(color_series.name, rotation=0, labelpad=10, va='center', ha='right')

    if remove_ticks:
        ax.tick_params(length=0)

    if remove_borders:
        for spine in ['bottom', 'top', 'left', 'right']:
            ax.spines[spine].set_visible(False)

    return ax

def plot_color_annotations_palette(
    val_vector,
    palette,
    ax=None,
    nan_color='#ffffff',
    hide_ticks=True,
    hide_borders=True,
    **kwargs,
):
    return plot_color_annotations(
        val_vector.map(palette),
        ax=ax,
        missing_color=nan_color,
        remove_ticks=hide_ticks,
        remove_borders=hide_borders,
        **kwargs,
    )

def plot_signature_kde(genesets_dict, exp_tpm, nrows=2, ncols=4, figsize=(20, 10)):
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=False, sharey=False)
    axes = axes.flatten()

    for idx, clus in enumerate(genesets_dict.keys()):
        ax = axes[idx]
        num_genes = len(genesets_dict[clus])
        for gene in genesets_dict[clus]:
            if gene in exp_tpm.columns:
                sns.kdeplot(np.log2(exp_tpm[gene] + 1), ax=ax, fill=True)
        ax.set_title(f'Signature {clus}, {num_genes} genes')
        ax.set_xlabel('log2(TPM + 1)')

    plt.tight_layout()
    plt.show()