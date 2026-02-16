import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import copy
from src.utils import intersect_df, assign_quantiles, intersect_series_with_index
from src.utils import generate_color_palette

from lifelines import KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts
from lifelines.statistics import logrank_test, multivariate_logrank_test

_UNIT_FACTORS = {
    ('Days',   'Weeks'):  7.0,
    ('Days',   'Months'): 365.25 / 12,
    ('Days',   'Years'):  365.25,
    ('Months', 'Days'):   30.0,
    ('Months', 'Weeks'):  30.0 / 7.0,
    ('Months', 'Years'):  12.0,
}

class SurvivalInfo:
    """Holds aligned time-to-event data for downstream KM analysis."""

    def __init__(
        self,
        times: pd.Series,
        events: pd.Series,
        units='Days',
        survival_type='PFS',
        max_time=None,
    ):
        clean_times = times.dropna()
        clean_events = events.dropna()
        self.times, self.events = intersect_df([clean_times, clean_events])

        self.times = self.times.astype(float)
        self.events = self.events.astype(float).astype(bool)

        if max_time:
            self.events = (self.times < max_time) & self.events
            self.times = self.times.clip(upper=max_time)

        self.units = units
        self.survival_type = survival_type
        self.index = self.times.index
        self.data = pd.DataFrame({'times': self.times, 'events': self.events})

def create_survival_annotation(
    times: pd.Series,
    events: pd.Series,
    max_time=None,
    in_units='Days',
    out_units='Months',
    survival_type='PFS',
) -> SurvivalInfo:
    converted = times.copy()

    if in_units != out_units:
        divisor = _UNIT_FACTORS.get((in_units, out_units))
        converted = converted / divisor

    return SurvivalInfo(converted, events, out_units, survival_type, max_time)

def _align_groups(survival_index, data):
    """Return aligned group series filtered to non-null, sorted categories."""
    aligned = intersect_series_with_index(survival_index, data)
    cats = sorted(aligned.dropna().unique())
    return aligned[aligned.isin(cats)], cats


def _compute_pvalue(survival, aligned_groups, categories):
    """Compute log-rank p-value for 2 groups or multivariate for >2."""
    n_cats = len(categories)
    if n_cats == 2:
        mask_a = aligned_groups == categories[0]
        mask_b = aligned_groups == categories[1]
        idx_a = aligned_groups[mask_a].index
        idx_b = aligned_groups[mask_b].index
        return logrank_test(
            survival.times[idx_a],
            survival.times[idx_b],
            event_observed_A=survival.events[idx_a],
            event_observed_B=survival.events[idx_b],
        ).p_value
    elif n_cats > 2:
        return multivariate_logrank_test(
            survival.times[aligned_groups.index],
            aligned_groups,
            survival.events[aligned_groups.index],
        ).p_value
    return None


def _fit_and_draw_group(kmf, survival, group_idx, color, ax, ci_show):
    kmf.fit(
        survival.times[group_idx],
        survival.events[group_idx],
        label='',
    )
    kmf.plot_survival_function(
        ax=ax,
        ci_show=ci_show,
        show_censors=True,
        c=color,
    )
    return survival.times[group_idx].max()


# main plotting function
def plot_kaplan_meier(
    data: pd.Series,
    survival: SurvivalInfo,
    title='',
    palette=None,
    pvalue=True,
    ax=None,
    figsize=(4, 4.5),
    p_digits=3,
    cmap=plt.cm.rainbow,
    max_time=None,
    legend='in',
    title_n_samples=False,
    ci_show=False,
    add_at_risk=True,
    **kwargs,
):
    aligned_groups, order = _align_groups(survival.index, data)

    track_max = max_time is None
    if track_max:
        max_time = 0

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    colors = (
        copy.copy(palette) if palette is not None
        else generate_color_palette(pd.Series(order), cmap=cmap)
    )

    kmf = KaplanMeierFitter()
    fitted_models = []

    for grp in order:
        members = aligned_groups[aligned_groups == grp]
        if members.empty:
            continue

        group_max = _fit_and_draw_group(
            kmf, survival, members.index, colors[grp], ax, ci_show,
        )

        # overlay the legend label via the last plotted line
        ax.lines[-1].set_label(str(grp))

        if track_max:
            max_time = max(max_time, group_max)

        snapshot = copy.copy(kmf)
        snapshot._label = grp
        fitted_models.append(snapshot)

    parts = [title]
    if title_n_samples:
        parts.append(f'N={len(aligned_groups)}')
    title_str = ''.join(parts)

    if pvalue:
        p = _compute_pvalue(survival, aligned_groups, order)
        if p is not None:
            sep = '\n' if title_str else ''
            title_str += f'{sep}p={p:.{p_digits}g}'

    ax.set_title(title_str, loc='left')

    if legend == 'in':
        ax.legend(loc='best')
    elif legend == 'out':
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    # at-risk table
    if add_at_risk and fitted_models:
        add_at_risk_counts(*fitted_models, ax=ax, rows_to_show=['At risk'])

    ax.set_xlim(0, max_time)

    return ax


def calculate_logrank_stats(
    data: pd.Series,
    survival: SurvivalInfo,
):
    aligned_groups, categories = _align_groups(survival.index, data)

    counts = {cat: int((aligned_groups == cat).sum()) for cat in categories}

    p_value = None
    if len(categories) == 2:
        idx_a = aligned_groups[aligned_groups == categories[0]].index
        idx_b = aligned_groups[aligned_groups == categories[1]].index
        if len(idx_a) > 0 and len(idx_b) > 0:
            p_value = logrank_test(
                survival.times[idx_a],
                survival.times[idx_b],
                event_observed_A=survival.events[idx_a],
                event_observed_B=survival.events[idx_b],
            ).p_value
        else:
            print('> 2 categories, skipped')

    return {
        'N': len(aligned_groups),
        'categories': categories,
        'counts': counts,
        'p_value': p_value,
    }


def plot_kaplan_meier_quantiles(
    data: pd.Series,
    survival: SurvivalInfo,
    q=(0.5,),
    cmap=plt.cm.Greens,
    palette=None,
    show_pvalue=True,
    **kwargs,
):
    aligned_data = intersect_series_with_index(survival.index, data)
    quantile_data = assign_quantiles(aligned_data, q)

    if palette is None:
        palette = generate_color_palette(quantile_data, cmap=cmap, min_v=0.4)

    kwargs['pvalue'] = show_pvalue
    return plot_kaplan_meier(quantile_data, survival, palette=palette, **kwargs)


def logrank_on_quantiles(
    data: pd.Series,
    survival: SurvivalInfo,
    q=(0.5,),
):
    aligned_data = intersect_series_with_index(survival.index, data)
    quantile_data = assign_quantiles(aligned_data, q)
    groups = sorted(quantile_data.dropna().unique())

    p_value = None
    if len(groups) == 2:
        idx_a = quantile_data[quantile_data == groups[0]].index
        idx_b = quantile_data[quantile_data == groups[1]].index

        if len(idx_a) > 0 and len(idx_b) > 0:
            p_value = logrank_test(
                survival.times[idx_a],
                survival.times[idx_b],
                event_observed_A=survival.events[idx_a],
                event_observed_B=survival.events[idx_b],
            ).p_value

    return len(quantile_data), p_value