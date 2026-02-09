import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import umap


def plot_umap(features, umap_kwargs=None, color_vec=None, color_pallet=None, scale=False, title: str = None):
    if scale:
        raise NotImplementedError
    umap_kwargs = {} if umap_kwargs is None else umap_kwargs
    default_umap_kwargs = {"random_state": 42,
                           "n_neighbors": 5,
                           }
    default_umap_kwargs.update(umap_kwargs)
    color_pallet = sns.color_palette() if color_pallet is None else color_pallet

    color = None
    if color_vec is not None:
        color_map = {c: i for i, c in enumerate(pd.Series(color_vec).unique())}
        color = [color_pallet[x] for x in color_vec.map(color_map)]

    reducer = umap.UMAP(**default_umap_kwargs)
    embedding = reducer.fit_transform(features)
    plt.scatter(embedding[:, 0], embedding[:, 1], c=color)
    plt.gca().set_aspect('equal', 'datalim')
    if title is not None:
        plt.title(f'{title}', fontsize=24)

