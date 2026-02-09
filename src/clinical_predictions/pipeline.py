from functools import partial

import numpy as np
import optuna
from sklearn.model_selection import ShuffleSplit, cross_val_score
from xgboost import XGBClassifier


def objective(trial, X_train, y_train):
    param = {
        'silent': 1,
        'objective': 'binary:logistic',
        'booster': trial.suggest_categorical('booster', ['gbtree', 'gblinear', 'dart']),
        'n_estimators': trial.suggest_int('n_estimators', 1, 8),
        'max_depth': trial.suggest_int('n_estimators', 1, 8),
        # 'lambda': trial.suggest_float('lambda', 1e-8, 1.0, log=True),
        # 'alpha': trial.suggest_float('alpha', 1e-8, 1.0, log=True)
    }

    model = XGBClassifier(*param)
    cv = ShuffleSplit(n_splits=5, test_size=0.3, random_state=0)
    scores = cross_val_score(model, X_train, y_train, cv=cv)
    mean_f1_score = np.mean(scores)
    trial.set_user_attr(key="best_booster", value=model)
    return mean_f1_score


def callback(study, trial):
    if study.best_trial.number == trial.number:
        study.set_user_attr(key="best_booster", value=trial.user_attrs["best_booster"])


def get_best_model_with_optuna(X_train, y_train):
    study = optuna.create_study(direction='maximize')
    study.optimize(partial(objective, X_train=X_train, y_train=y_train), n_trials=10, callbacks=[callback])
    best_model = study.user_attrs["best_booster"]
    return best_model
