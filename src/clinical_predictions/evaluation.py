from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.base import ClassifierMixin
from sklearn.metrics import classification_report, accuracy_score, f1_score


@dataclass
class ExperimentData:
    X_train: pd.DataFrame
    X_test: pd.DataFrame
    y_train: pd.Series
    y_test: pd.Series
    y_pred: pd.Series
    model: ClassifierMixin
    report: dict


def train_and_eval_model(X_train, X_test, y_train, y_test, model, extra_for_report=None,
                         baseline="RANDOM_CHOICE", base_line_n_test: int = 10) -> ExperimentData:
    model.fit(X_train, y_train)
    y_pred = pd.Series(model.predict(X_test))
    report = {}
    if extra_for_report is not None:
        report["extra"] = extra_for_report

    if baseline == "RANDOM_CHOICE":
        if base_line_n_test < 1:
            raise ValueError('baseline number of iterations must be larger than 1')
        running_accuracy = 0
        running_f1_score = 0
        for _ in range(base_line_n_test):
            probas = y_train.value_counts(normalize=True)
            y_pred_baseline = np.random.choice(len(probas), len(y_test), p=[probas.loc[0], probas.loc[1]])
            running_accuracy += accuracy_score(y_true=y_test, y_pred=y_pred_baseline)
            running_f1_score += f1_score(y_true=y_test, y_pred=y_pred_baseline)

        report["baseline"] = classification_report(y_true=y_test, y_pred=y_pred_baseline, output_dict=True)
        report["baseline"]['1']["accuracy"] = running_accuracy / base_line_n_test
        report["baseline"]['1']["f1-score"] = running_f1_score / base_line_n_test

    report.update(classification_report(y_true=y_test, y_pred=y_pred, output_dict=True))
    report['1']["accuracy"] = accuracy_score(y_true=y_test, y_pred=y_pred)

    return ExperimentData(X_train, X_test, y_train, y_test, y_pred, model, report)


# results visualization
def _add_suffix_to_index(df, sufix):
    df.index = [str(ind) + f"_{sufix}" for ind in df.index]
    return df


def generate_datasets_summerization(monitor):
    df_y_train = pd.concat([exp_data.y_train.value_counts().rename(f"{treatment}")
                            for treatment, exp_data in monitor.items()], axis=1)
    df_y_train = _add_suffix_to_index(df_y_train, "y_train")

    df_y_test = pd.concat([exp_data.y_test.value_counts().rename(f"{treatment}")
                           for treatment, exp_data in monitor.items()], axis=1)
    df_y_test = _add_suffix_to_index(df_y_test, "y_test")

    df_y_pred = pd.concat([exp_data.y_pred.value_counts().rename(f"{treatment}")
                           for treatment, exp_data in monitor.items()], axis=1)
    df_y_pred = _add_suffix_to_index(df_y_pred, "y_pred")

    df_report = pd.concat([pd.Series(exp_data.report['1']).rename(treatment)
                           for treatment, exp_data in monitor.items()], axis=1)
    df_report = _add_suffix_to_index(df_report, "report")

    df_report_baseline = pd.concat([pd.Series(exp_data.report["baseline"]['1']).rename(treatment)
                                    for treatment, exp_data in monitor.items()], axis=1)
    df_report_baseline = _add_suffix_to_index(df_report_baseline, "baseline")

    series_model = pd.Series({treatment: str(exp_data.model) for treatment, exp_data in monitor.items()}).rename(
        "model")

    extras = pd.Series({treatment: str(exp_data.report["extra"]) if "extra" in exp_data.report else ""
                        for treatment, exp_data in monitor.items()}).rename("extra")

    return pd.concat([df_y_train.T, df_y_test.T, df_y_pred.T, df_report.T, df_report_baseline.T, series_model, extras],
                     axis=1)
