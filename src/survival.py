import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import warnings
import copy
from src.utils import intersect_df, assign_quantiles, intersect_series_with_index
from src.utils import generate_color_palette

from lifelines.plotting import add_at_risk_counts
from lifelines.statistics import logrank_test
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test

class SurvivalInfo:
    def __init__(
        self,
        times: pd.Series,
        events: pd.Series,
        units='Days',
        survival_type='PFS',
        max_time=None,
    ):
        self.times, self.events = intersect_df(
            [times.dropna(), events.dropna()]
        )

        self.times = self.times.astype(float)
        self.events = self.events.astype(float).astype(bool)

        if max_time:
            self.events = np.logical_and(self.times < max_time, self.events)
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
    times_copy = times.copy()

    UNIT_CONVERSIONS = {
        ('Days', 'Weeks'): 7.0,
        ('Days', 'Months'): 365.25 / 12,
        ('Days', 'Years'): 365.25,
        ('Months', 'Days'): 30.0,
        ('Months', 'Weeks'): 30.0 / 7.0,
        ('Months', 'Years'): 12.0,
    }

    if in_units != out_units:
        factor = UNIT_CONVERSIONS.get((in_units, out_units))
        times = times / factor
    
    return SurvivalInfo(times, events, out_units, survival_type, max_time)

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
    kmf = KaplanMeierFitter()

    auto_max_time = False
    if max_time is None:
        max_time = 0
        auto_max_time = True

    aligned_groups = intersect_series_with_index(survival.index, data)
    order = list(sorted(aligned_groups.dropna().unique()))
    aligned_groups = aligned_groups[aligned_groups.isin(order)]

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    if palette is None:
        color_scheme = generate_color_palette(pd.Series(order), cmap=cmap)
    else:
        color_scheme = copy.copy(palette)

    fitted_models = []
    for group in order:
        subset = aligned_groups[aligned_groups == group]
        if len(subset):
            kmf.fit(
                survival.times[subset.index], survival.events[subset.index], label=''
            )
            kmf.plot_survival_function(
                ax=ax,
                ci_show=ci_show,
                show_censors=True,
                c=color_scheme[group],
                label=str(group),
            )

            if auto_max_time:
                max_time = max(max_time, survival.times[subset.index].max())

            kmf._label = group
            fitted_models.append(copy.copy(kmf))
    if title_n_samples:
        title += f'N={len(aligned_groups)}'

    if pvalue:
        if len(title):
            title += '\n'
        if len(order) == 2:
            group1 = aligned_groups[aligned_groups == order[0]]
            group2 = aligned_groups[aligned_groups == order[1]]

            p_value = logrank_test(
                survival.times[group1.index],
                survival.times[group2.index],
                event_observed_A=survival.events[group1.index],
                event_observed_B=survival.events[group2.index],
            ).p_value
            title += f'p={p_value:.{p_digits}g}'
        elif len(order) > 2:
            p_value = multivariate_logrank_test(
                survival.times[aligned_groups.index],
                aligned_groups,
                survival.events[aligned_groups.index]
            ).p_value
            title += f'p={p_value:.{p_digits}g}'

    ax.set_title(title, loc = 'left')   

    if legend == 'in':
        ax.legend(loc='best')
    elif legend == 'out':
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    if add_at_risk:
        add_at_risk_counts(*fitted_models, ax=ax, rows_to_show=['At risk'])

    if max_time is not None or auto_max_time:
        ax.set_xlim(0, max_time)
    
    return ax

def calculate_logrank_stats(
    data: pd.Series,
    survival: SurvivalInfo,
):
    aligned_groups = intersect_series_with_index(survival.index, data)
    
    categories = list(sorted(aligned_groups.dropna().unique()))

    aligned_groups = aligned_groups[aligned_groups.isin(categories)]
    
    result = {
        'N': len(aligned_groups),
        'categories': categories,
        'counts': {cat: sum(aligned_groups == cat) for cat in categories},
        'p_value': None
    }
    
    # Calculate p-value if we have exactly 2 categories
    if len(categories) == 2:
        group1 = aligned_groups[aligned_groups == categories[0]]
        group2 = aligned_groups[aligned_groups == categories[1]]
        
        if len(group1) > 0 and len(group2) > 0:
            result['p_value'] = logrank_test(
                survival.times[group1.index],
                survival.times[group2.index],
                event_observed_A=survival.events[group1.index],
                event_observed_B=survival.events[group2.index],
            ).p_value

        else:
            print('> 2 categories, skipped')
    
    return result

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
    order = list(sorted(quantile_data.dropna().unique()))

    p_value = None
    if len(order) == 2:
        group1 = quantile_data[quantile_data == order[0]]
        group2 = quantile_data[quantile_data == order[1]]

        if not group1.empty and not group2.empty:
            p_value = logrank_test(
                survival.times[group1.index],
                survival.times[group2.index],
                event_observed_A=survival.events[group1.index],
                event_observed_B=survival.events[group2.index],
            ).p_value

    return len(quantile_data), p_value