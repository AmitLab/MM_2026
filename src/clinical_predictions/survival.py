from lifelines.statistics import logrank_test
from lifelines import KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import warnings
import copy
from src.utils import assign_quantiles, generate_color_palette

def intersect_samples(series_list):
    if not series_list:
        return series_list
    common_idx = series_list[0].index
    for s in series_list[1:]:
        common_idx = common_idx.intersection(s.index)
    return [s.loc[common_idx] for s in series_list]

class SurvivalInfo:
    def __init__(
        self,
        times: pd.Series,
        events: pd.Series,
        units='Days',
        survival_type='OS',
        max_time=None,
    ):
        self.times, self.events = intersect_samples(
            [times.dropna(), events.dropna()]
        )
        if len(self.times) == 0:
            raise Exception('Different id\'s provided')

        self.times = self.times.astype(float)
        try:
            self.events = self.events.astype(float).astype(bool)
        except Exception as e:
            raise Exception('Bad censorship ' + str(e))

        if max_time:
            self.events = np.logical_and(self.times < max_time, self.events)
            self.times = self.times.clip(upper=max_time)

        self.units = units
        self.survival_type = survival_type
        self.index = self.times.index

        self.data = pd.DataFrame({'times': self.times, 'events': self.events})

    def plot(self, ax=None, label='Samples', **kwargs):
        """
        Simple KM plot of the all samples
        """
        kwargs['palette'] = kwargs.get('palette', {label: 'black'})
        plot_kaplan_meier(
            pd.Series(label, index=self.index), self, ax=ax, pvalue=False, **kwargs
        )

    def intersect_index(self, data: pd.Series) -> pd.Series:
        """
        Returns a subset of data with common indexes
        """
        data_common = data.loc[self.index.intersection(data.index)]
        discarded_samples = len(data) - len(data_common)
        if discarded_samples:
            warnings.warn(
                f'{discarded_samples} out of {len(data)} discarded due to no survival annotation'
            )
        return data_common

    def event_at_time(self, time: float) -> pd.Series:
        """
        Returns a series with Event/Censored/No event for each sample at a specific time
        """
        return pd.concat(
            [
                self.data[self.data.times <= time].events.map(
                    {True: 'Event', False: 'Censored'}
                ),
                pd.Series('No event', index=self.data[self.data.times > time].index),
            ]
        )

def create_survival_annotation(
    times: pd.Series,
    events: pd.Series,
    max_time=None,
    in_units='Days',
    out_units='Months',
    survival_type='OS',
) -> SurvivalInfo:
    """
    Prepare SurvivalInfo with survival annotation suitable for plot_kaplan_meier()
    :param times: Series with survival in days
    :param events: Series with Death events as True/False or 1/0
    :param max_time: Limit survival time and events (in out units)
    :param in_units: Should be Days, others not implemented
    :param out_units: To convert days into Weeks/Months/Years
    :param survival_type:
    :return: pd.DataFrame, index - patients, columns - ['times', 'events']
    """
    times_copy = times.copy()

    conversion_factors = {
        ('Days', 'Weeks'): 7.0,
        ('Days', 'Months'): 365.25 / 12,
        ('Days', 'Years'): 365.25,
        ('Months', 'Days'): 30.0,
        ('Months', 'Weeks'): 30.0 / 7.0,
        ('Months', 'Years'): 12.0,
    }

    if in_units != out_units:
        try:
            factor = conversion_factors[(in_units, out_units)]
            times_copy /= factor
        except KeyError:
            raise NotImplementedError(
                f'{in_units} to {out_units} conversion not supported. Try constructing SurvivalInfo by yourself'
            )

    return SurvivalInfo(times_copy, events, out_units, survival_type, max_time)

def plot_kaplan_meier(
    groups: pd.Series,
    survival: SurvivalInfo,
    loglogs=False,
    title='',
    palette=None,
    pvalue=True,
    ax=None,
    figsize=(4, 4.5),
    p_digits=3,
    order=None,
    cmap=plt.cm.rainbow,
    max_time=None,
    legend='in',
    ci_show=False,
    title_n_samples=False,
    add_at_risk=True,
    weightings=None,
    **kwargs,
):
    kmf = KaplanMeierFitter()

    emt = False
    if max_time is None:
        max_time = 0
        emt = True

    groups_common = survival.intersect_index(groups)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    if order is None:
        order = list(sorted(groups_common.dropna().unique()))

    groups_common = groups_common[groups_common.isin(order)]

    if palette is None:
        color_palette = generate_color_palette(pd.Series(order), cmap=cmap)
    else:
        color_palette = copy.copy(palette)

    kmfs = []
    for group_name in order:
        samples = groups_common[groups_common == group_name]
        if len(samples):
            kmf.fit(
                survival.times[samples.index], survival.events[samples.index], label=''
            )
            if loglogs:
                kmf.plot_loglogs(
                    ax=ax,
                    show_censors=True,
                    c=color_palette[group_name],
                    label=str(group_name),
                )
            else:
                kmf.plot_survival_function(
                    ax=ax,
                    ci_show=ci_show,
                    show_censors=True,
                    c=color_palette[group_name],
                    label=str(group_name),
                )

            if emt:
                max_time = max(max_time, survival.times[samples.index].max())

            kmf._label = group_name
            kmfs.append(copy.copy(kmf))

    if len(title):
        title += ', '
    if title_n_samples:
        title += f'N={len(groups_common)}'

    if not loglogs and pvalue:
        if len(title):
            title += '\n'
        if len(order) == 2:
            group1 = groups_common[groups_common == order[0]]
            group2 = groups_common[groups_common == order[1]]

            pval = logrank_test(
                survival.times[group1.index],
                survival.times[group2.index],
                event_observed_A=survival.events[group1.index],
                event_observed_B=survival.events[group2.index],
                weightings=weightings,
            ).p_value
            title += f'p={pval:.{p_digits}g}'

    ax.set_title(title)
    if legend == 'in':
        ax.legend(loc='best')
    elif legend == 'out':
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    if add_at_risk:
        add_at_risk_counts(*kmfs, ax=ax, rows_to_show=['At risk'])

    return ax

def plot_kaplan_meier_quantiles(
    data: pd.Series,
    survival: SurvivalInfo,
    q=(0.5,),
    cmap=plt.cm.Oranges,
    min_v=0.4,
    palette=None,
    show_pvalue=True,  # New parameter
    **kwargs,
):
    data_common = survival.intersect_index(data)
    quantile_data = assign_quantiles(data_common, q)

    if 'cmap' in kwargs:
        cmap = kwargs['cmap']

    if 'min_v' in kwargs:
        min_v = kwargs['min_v']

    if palette is None:
        palette = generate_color_palette(quantile_data, cmap=cmap, min_v=min_v, sort=False)

    # Set pvalue parameter based on show_pvalue
    kwargs['pvalue'] = show_pvalue

    return plot_kaplan_meier(quantile_data, survival, palette=palette, **kwargs)

def calculate_kaplan_meier_quantiles(
    data: pd.Series,
    survival: SurvivalInfo,
    q=(0.5,),
    weightings=None,
):
    data_common = survival.intersect_index(data)
    quantile_data = assign_quantiles(data_common, q)

    kmf = KaplanMeierFitter()

    order = list(sorted(quantile_data.dropna().unique()))
    quantile_data = quantile_data[quantile_data.isin(order)]

    kmfs = []
    for group_name in order:
        samples = quantile_data[quantile_data == group_name]
        if len(samples):
            kmf.fit(
                survival.times[samples.index], survival.events[samples.index], label=''
            )
            kmf._label = group_name
            kmfs.append(copy.copy(kmf))

    pval = None
    if len(order) == 2:
        group1 = quantile_data[quantile_data == order[0]]
        group2 = quantile_data[quantile_data == order[1]]

        if not group1.empty and not group2.empty:
            pval = logrank_test(
                survival.times[group1.index],
                survival.times[group2.index],
                event_observed_A=survival.events[group1.index],
                event_observed_B=survival.events[group2.index],
                weightings=weightings,
            ).p_value

    return len(quantile_data), pval
