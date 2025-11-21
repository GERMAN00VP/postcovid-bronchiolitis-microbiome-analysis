import numpy as np
import pandas as pd
from sklearn.linear_model import LassoCV
from sklearn.metrics import roc_auc_score

from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

from bootstrap_laso.clr_preprocess import PreprocessModule   # IMPORTA TU PREPROCESS



class BootstrapLassoRunner:
    """
    Ejecuta bootstrap -> preprocess -> LassoCV -> AUC en OOB sample.
    Ideal para repetirlo muchas veces o para permutación posterior.
    """

    def __init__(self, 
                 clr_variables=None,
                 standardize_variables=None,
                 n_splits_cv=5,
                 random_state=None):
        
        self.clr_variables = clr_variables
        self.standardize_variables = standardize_variables
        self.n_splits_cv = n_splits_cv
        self.random_state = random_state

    def run_once(self, X, y):
        """
        Hace un bootstrap estratificado adaptado para clases desbalanceadas + entrenamiento + evaluación en OOB.
        Devuelve AUC, betas del modelo, lambda seleccionado, n_features seleccionadas y coef de features.
        """

        X_boot_list = []
        y_boot_list = []

        rng = np.random.default_rng(self.random_state)  # reproducible

        # -------- 1) Stratified bootstrap con OOB seguro ----------
        for cls in y.unique():
            X_cls = X[y == cls]
            y_cls = y[y == cls]
            n_total = len(X_cls)
            n_bootstrap = max(n_total - 1, 1)  # deja al menos 1 fuera
            indices = rng.integers(0, n_total, size=n_bootstrap)
            X_res = X_cls.iloc[indices]
            y_res = y_cls.iloc[indices]
            X_boot_list.append(X_res)
            y_boot_list.append(y_res)

        X_boot = pd.concat(X_boot_list)
        y_boot = pd.concat(y_boot_list)

        # -------- 2) Out-of-bag samples ----------
        in_boot_idx = X_boot.index
        oob_idx = X.index.difference(in_boot_idx)

        X_oob = X.loc[oob_idx]
        y_oob = y.loc[oob_idx]

        # Guardar info de OOB
        n_oob = len(y_oob)
        oob_class_counts = y_oob.value_counts().to_dict()

        # chequeo clases
        if y_oob.nunique() < 2:
            return None

        # -------- 3) Preprocess y Lasso ----------
        prep = PreprocessModule(
            clr_variables=self.clr_variables,
            standardize_variables=self.standardize_variables
        )
        X_boot_prep = prep.fit_transform(X_boot)
        X_oob_prep = prep.transform(X_oob)

        model = LassoCV(cv=self.n_splits_cv, random_state=self.random_state)
        model.fit(X_boot_prep, y_boot)

        preds = model.predict(X_oob_prep)
        auc = roc_auc_score(y_oob, preds)

        # -------- 4) Número de features seleccionadas --------
        coef_series = pd.Series(model.coef_, index=X.columns)
        n_features_selected = (coef_series != 0).sum()

        out = {
            "auc": auc,
            "model": model,
            "lambda": model.alpha_,
            "coef": coef_series,
            "n_features_selected": n_features_selected,
            "Features": model.coef_,
            "intercept": model.intercept_,
            "oob_index": oob_idx,
            "n_oob": n_oob,
            "oob_class_counts": oob_class_counts
            
        }

        return out




    def run_many(self, X, y, n_boot=500):
        """
        Corre run_once() muchas veces y devuelve lista de resultados.
        """

        results = []

        for i in range(n_boot):
            r = self.run_once(X, y)
            if r is not None:
                results.append(r)

        return results
