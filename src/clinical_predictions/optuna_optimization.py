from functools import partial
from typing import List, Optional

import numpy as np
import optuna
import sklearn
from scipy.stats import hmean
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SequentialFeatureSelector
from sklearn.model_selection import ShuffleSplit, cross_validate
from xgboost import XGBClassifier

from clinical_predictions.utils import balanced_subsample


def classifiaction_cv_objective(trial, X_train, y_train, use_feature_selection: bool = False,
                                try_balance_with_subsample: bool = False,
                                classifier_names: Optional[List] = None, precision_alpha: float = 0.8):
    classifier_names = classifier_names if classifier_names is not None else ['RandomForest',
                                                                              # 'XGBoost',
                                                                              'LogisticRegression']  # ['LogisticRegression', 'SVC','RandomForest', 'XGBoost', 'LogisticRegression'])
    classifier_name = trial.suggest_categorical('classifier', classifier_names)
    if classifier_name == 'SVC':
        svc_c = trial.suggest_float('svc_c', 1e-1, 1e3, log=True)
        model = sklearn.svm.SVC(C=svc_c, gamma='auto', probability=True)
    elif classifier_name == 'RandomForest':
        rf_max_depth = trial.suggest_int('rf_max_depth', 2, 3, log=True)
        rf_n_estimators = trial.suggest_int('rf_n_estimators', 2, 64, log=True)
        rf_max_samples = trial.suggest_categorical("rf_max_samples", [0.8])  # [None, 0.8]
        model = RandomForestClassifier(max_depth=rf_max_depth, n_estimators=rf_n_estimators, max_samples=rf_max_samples)
    elif classifier_name == 'LogisticRegression':
        logistic_regression_c = trial.suggest_float('logistic_regression_c', 1e-4, 1e2, log=True)
        logr_penalty = trial.suggest_categorical('logr_penalty', ["l1", "l2"])
        class_weight = trial.suggest_categorical('class_weight', [None, "balanced"])
        model = sklearn.linear_model.LogisticRegression(C=logistic_regression_c, penalty=logr_penalty,
                                                        solver='liblinear', class_weight=class_weight)
    else:
        param = {
            'booster': trial.suggest_categorical('xgb_booster', ['gbtree', 'dart']),
            'n_estimators': trial.suggest_int('xgb_n_estimators', 2, 32, log=True),
            'max_depth': trial.suggest_int('xgb_max_depth', 1, 3),
            # 'alpha': trial.suggest_float('alpha', 1e-8, 1.0, log=True)
        }
        model = XGBClassifier(**param)

    if try_balance_with_subsample:
        subsample = trial.suggest_categorical('subsample', [True, False])
        if subsample:
            subsample_indexs = balanced_subsample(y_train)
            X_train = X_train.loc[subsample_indexs]
            y_train = y_train.loc[subsample_indexs]

    cv = ShuffleSplit(n_splits=8, test_size=0.3, random_state=0)

    if use_feature_selection:
        fix_feature_selection = trial.suggest_categorical('fix_feature_selection', [True, False])
        alpha_features_to_select = 'auto'
        if fix_feature_selection:
            alpha_features_to_select = trial.float("alpha_features_to_select", 0.1, 0.9)
        sfs = SequentialFeatureSelector(model, scoring='f1_weighted', n_features_to_select=alpha_features_to_select,
                                        cv=cv, n_jobs=4)
        X_train = sfs.fit_transform(X_train, y_train)

    scores = cross_validate(model, X_train, y_train, cv=cv,
                            scoring=['accuracy', 'precision', 'f1', 'f1_weighted', 'f1_macro'])

    mean_f1_prec_score = hmean(scores['test_f1_weighted'] + scores['test_precision'] * precision_alpha)
    # mean_f1_prec_score = np.mean(scores['test_f1_weighted'] + scores['test_precision'] * precision_alpha)
    final_score = mean_f1_prec_score

    trial.set_user_attr(key="best_booster", value=model)
    trial.set_user_attr(key="scores", value=scores)
    return final_score


def regression_cv_objective(trial, X_train, y_train):
    regressor_names = ['RandomForestRegressor', 'SVR']  # LogisticRegression
    regressor_name = trial.suggest_categorical('regressor', regressor_names)
    if regressor_name == 'SVR':
        svc_c = trial.suggest_float('svc_c', 1e-1, 1e3, log=True)
        model = sklearn.svm.SVC(C=svc_c, gamma='auto')
    elif regressor_name == 'RandomForestRegressor':
        rf_max_depth = trial.suggest_int('rf_max_depth', 2, 3, log=True)
        rf_n_estimators = trial.suggest_int('rf_n_estimators', 2, 64, log=True)
        model = RandomForestClassifier(max_depth=rf_max_depth, n_estimators=rf_n_estimators)
    else:
        raise ValueError("model name error")
    cv = ShuffleSplit(n_splits=8, test_size=0.3, random_state=0)
    scores = cross_validate(model, X_train, y_train, cv=cv, scoring=['r2'])
    final_score = np.mean(scores['test_r2'])
    trial.set_user_attr(key="best_booster", value=model)
    trial.set_user_attr(key="scores", value=scores)
    return final_score


def save_best_model_callback(study, trial):
    if study.best_trial.number == trial.number:
        study.set_user_attr(key="best_booster", value=trial.user_attrs["best_booster"])


def get_best_model_with_optuna(X_train, y_train, n_trials=30, use_regression_score=False, **extra_params_for_hp_search):
    if not use_regression_score:
        study = optuna.create_study(direction='maximize')
        study.optimize(
            partial(classifiaction_cv_objective, X_train=X_train, y_train=y_train, **extra_params_for_hp_search),
            n_trials=n_trials, callbacks=[save_best_model_callback])
    else:
        study = optuna.create_study(direction='maximiza')
        study.optimize(
            partial(regression_cv_objective, X_train=X_train, y_train=y_train, **extra_params_for_hp_search),
            n_trials=n_trials, callbacks=[save_best_model_callback])
    best_model = study.user_attrs["best_booster"]
    best_trail = study.best_trial
    return best_model, best_trail
