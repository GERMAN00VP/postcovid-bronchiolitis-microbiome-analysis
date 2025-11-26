import numpy as np
import pandas as pd
from scipy import stats

import copy
from concurrent.futures import ProcessPoolExecutor, as_completed
from scipy.stats import wasserstein_distance

from sklearn.linear_model import LassoCV, RidgeCV, ElasticNetCV
from sklearn.metrics import roc_auc_score

from bootstrap_laso.clr_preprocess import PreprocessModule

# -------------------------------------------------------------------
# Wrapper necesario para ejecutar run_once con multiprocessing
# -------------------------------------------------------------------
def run_once_worker(args):
    """
    Wrapper pickleable para permitir paralelización real.
    """
    obj, X, y, seed = args
    obj.random_state = seed  # reproducibilidad por bootstrap
    return obj.run_once(X, y)

# -------------------------------------------------------------------
# Funciones de distancia de Wasserstein
# -------------------------------------------------------------------
def distance_matrix_wasserstein(list_of_distributions):
    # list_of_distributions: lista de arrays (cada array son las muestras para esa distribución)
    N = len(list_of_distributions)
    D = np.zeros((N, N))
    for i in range(N):
        for j in range(i+1, N):
            d = wasserstein_distance(list_of_distributions[i], list_of_distributions[j])
            D[i, j] = D[j, i] = d
    return D

def p_from_distance_matrix(D, real_idx=-1):
    N = D.shape[0]
    perm_idx = [i for i in range(N) if i != real_idx]
    d_real = D[real_idx, perm_idx].mean()
    d_perm = []
    for k in perm_idx:
        others = [i for i in perm_idx if i != k]
        d_perm.append(D[k, others].mean())

    d_perm = np.array(d_perm)
    p = np.mean(d_perm >= d_real)
    return p, d_real, d_perm

def extract_coefs(results, test_dict=None):
    """
    Construye un diccionario donde cada clave contiene una LISTA de distribuciones.
    test_dict[feat] = [dist_real, dist_perm1, dist_perm2, ...]
    """

    # ------------------------------------------------------
    # 1. Primera ejecución: crear contenedores
    # ------------------------------------------------------
    if test_dict is None:
        test_dict = {"features_to_keep": None}

        # ---- Determinar features no-siempre-cero en REAL ----
        all_coefs = [coef_series for coef_series in results['coef']]
        coefs_df = pd.DataFrame(all_coefs)

        non_zero_mask = (coefs_df != 0).any(axis=0)
        features_to_keep = coefs_df.columns[non_zero_mask].tolist()

        test_dict["features_to_keep"] = features_to_keep

        # Crear una lista vacía para cada feature
        for feat in features_to_keep:
            test_dict[feat] = []

        # Crear contenedor para AUC
        test_dict["AUC"] = []

        # Primera distribución real
        test_dict["AUC"].append(results["auc"].values)
        for feat in features_to_keep:
            test_dict[feat].append(coefs_df[feat].values)

        return test_dict

    # ------------------------------------------------------
    # 2. Ejecuciones posteriores (permutaciones)
    # ------------------------------------------------------
    features_to_keep = test_dict["features_to_keep"]

    # Distribución AUC permutada
    test_dict["AUC"].append(results["auc"].values)

    # Distribuciones coef permutadas
    perm_coefs = pd.DataFrame([coef for coef in results["coef"]])

    for feat in features_to_keep:
        test_dict[feat].append(perm_coefs[feat].values)

    return test_dict


# ===================================================================
#                    BOOTSTRAP + LASSO / RIDGE / EN
# ===================================================================
class BootstrapLassoRunner:

    def __init__(self, 
                 clr_variables=None,
                 standardize_variables=None,
                 n_splits_cv=5,
                 random_state=None,
                 penalty="l1",        # "l1", "l2", "elasticnet"
                 alpha=None,          # si no usas CV
                 l1_ratio=0.5,        # solo ElasticNet
                 alphas=None):        # grid de alphas para CV
        
        self.clr_variables = clr_variables
        self.standardize_variables = standardize_variables
        self.n_splits_cv = n_splits_cv
        self.random_state = random_state

        self.penalty = penalty.lower()
        self.alpha = alpha
        self.l1_ratio = l1_ratio
        self.alphas = alphas  # si None → rango por defecto

    # ----------------------------------------------------------------
    def _select_model(self):
        """
        Devuelve el modelo adecuado según el tipo de penalización.
        Con valores por defecto apropiados si alphas=None.
        """

        # ---- LASSO ----
        if self.penalty == "l1":
            return LassoCV(
                cv=self.n_splits_cv,
                random_state=self.random_state,
                alphas=self.alphas if self.alphas is not None 
                                   else np.logspace(-4, 4, 100)
            )

        # ---- RIDGE ----
        elif self.penalty == "l2":
            return RidgeCV(
                alphas=self.alphas if self.alphas is not None 
                                   else np.logspace(-4, 4, 100)
            )

        # ---- ELASTIC NET ----
        elif self.penalty == "elasticnet":
            return ElasticNetCV(
                cv=self.n_splits_cv,
                random_state=self.random_state,
                l1_ratio=self.l1_ratio,
                alphas=self.alphas if self.alphas is not None 
                                   else np.logspace(-4, 4, 100)
            )

        else:
            raise ValueError("penalty debe ser 'l1', 'l2' o 'elasticnet'.")

    # ----------------------------------------------------------------
    def run_once(self, X, y):
        """
        Bootstrap estratificado + preprocess + modelo penalizado + AUC OOB.
        """

        X_boot_list = []
        y_boot_list = []

        rng = np.random.default_rng(self.random_state)

        # -------- 1) Stratified bootstrap ----------
        for cls in y.unique():
            X_cls = X[y == cls]
            y_cls = y[y == cls]
            n_total = len(X_cls)
            n_bootstrap = max(n_total - 1, 1)

            idx = rng.integers(0, n_total, size=n_bootstrap)
            X_boot_list.append(X_cls.iloc[idx])
            y_boot_list.append(y_cls.iloc[idx])

        X_boot = pd.concat(X_boot_list)
        y_boot = pd.concat(y_boot_list)

        # -------- 2) OOB ----------
        oob_idx = X.index.difference(X_boot.index)
        if len(oob_idx) == 0:
            return None

        X_oob = X.loc[oob_idx]
        y_oob = y.loc[oob_idx]

        for cls in y.unique():
            if (y_oob == cls).sum() < 1:
                return None

        # -------- 3) Preprocessing ----------
        prep = PreprocessModule(
            clr_variables=self.clr_variables,
            standardize_variables=self.standardize_variables
        )
        X_boot_prep = prep.fit_transform(X_boot)
        X_oob_prep  = prep.transform(X_oob)

        # -------- 4) Modelo penalizado ----------
        model = self._select_model()
        model.fit(X_boot_prep, y_boot)

        # -------- 5) Predicción y AUC ----------
        preds = model.predict(X_oob_prep)
        auc = roc_auc_score(y_oob, preds)

        # -------- 6) Coeficientes ----------
        coef_series = pd.Series(model.coef_, index=X.columns)

        return {
            "auc": auc,
            "lambda": getattr(model, "alpha_", self.alpha),
            "coef": coef_series,
            "intercept": model.intercept_,
            "oob_class_counts": y_oob.value_counts().to_dict(),
        }

    # ----------------------------------------------------------------
    def run_many(self, X, y, n_boot=500, n_jobs=8):
        """
        Paralelización REAL usando ProcessPoolExecutor.
        Escala con todos los cores.
        """

        results = []
        base_obj = copy.deepcopy(self)

        # Seeds independientes por bootstrap
        seeds = np.random.default_rng(self.random_state).integers(0, 2**32 - 1, size=n_boot)

        with ProcessPoolExecutor(max_workers=n_jobs) as executor:

            futures = [
                executor.submit(run_once_worker, (copy.deepcopy(base_obj), X, y, seed))
                for seed in seeds
            ]

            for fut in as_completed(futures):
                r = fut.result()
                if r is not None:
                    results.append(r)

        return pd.DataFrame(results)

    # ----------------------------------------------------------------

    def empirical_p_values(self, X, y, results_real, n_boot=100, n_permutations=100,
                        random_state=None, n_jobs=8):

        rng = np.random.default_rng(random_state)

        # ----------------------------------------------------------
        # 1. Obtener distribuciones reales
        # ----------------------------------------------------------
        test_dict = extract_coefs(results_real)

        # ----------------------------------------------------------
        # 2. Obtener distribuciones permutadas
        # ----------------------------------------------------------
        for i in range(n_permutations):
            X_perm = X.reset_index(drop=True)
            y_perm = y.sample(frac=1.0, replace=False,
                            random_state=rng.integers(0, 2**32-1)).reset_index(drop=True)

            perm_results = self.run_many(X_perm, y_perm, n_boot=n_boot, n_jobs=n_jobs)
            test_dict = extract_coefs(results=perm_results, test_dict=test_dict)

        # ----------------------------------------------------------
        # 3. Calcular p-values por Wasserstein para cada feature
        # ----------------------------------------------------------
        results_out = []

        features = ["AUC"] + test_dict["features_to_keep"]

        for feat in features:
            # Lista de distribuciones
            dists = test_dict[feat]     # dists[0] es REAL
            D = distance_matrix_wasserstein(dists)
            p, d_real, d_perm = p_from_distance_matrix(D, real_idx=0)

            results_out.append({
                "Feature": feat,
                "p_value": p,
                "d_real": d_real,
                "d_perm_mean": np.mean(d_perm),
                "d_perm_std": np.std(d_perm)
            })

        return pd.DataFrame(results_out)
