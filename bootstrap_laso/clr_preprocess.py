import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.preprocessing import StandardScaler


class CLRTransformer(BaseEstimator, TransformerMixin):
    """
    CLR transformer aplicado solo a columnas específicas.
    """

    def __init__(self, variables=None):
        self.variables = variables

    def fit(self, X, y=None):
        return self

    def _clr_transform(self, df):
        def relab(row):
            total = row.sum()
            row = row.copy()
            row[row != 0] = row[row != 0] / total
            return row

        rel_df = df.apply(relab, axis=1)

        pseudocount = rel_df[rel_df > 0].min().min() / 100

        data = rel_df.to_numpy(float)
        data += pseudocount

        gm = np.exp(np.mean(np.log(data), axis=1, keepdims=True))

        clr_vals = np.log(data / gm)

        return pd.DataFrame(
            clr_vals, columns=df.columns, index=df.index
        )

    def transform(self, X):
        X = X.copy()

        if self.variables is None:
            return X

        clr_df = self._clr_transform(X[self.variables])

        # Replace only the requested columns
        X[self.variables] = clr_df

        return X


class PreprocessModule(BaseEstimator, TransformerMixin):
    """
    Módulo de preprocesamiento:
       1) CLR en columnas seleccionadas
       2) StandardScaler en columnas seleccionadas (después del CLR)
    """

    def __init__(self, clr_variables=None, standardize_variables=None):
        self.clr_variables = clr_variables
        self.standardize_variables = standardize_variables
        self.clr_transformer = None
        self.scaler = None

    def fit(self, X, y=None):
        X_tmp = X.copy()

        # --- CLR -----------------------------------------
        if self.clr_variables is not None:
            self.clr_transformer = CLRTransformer(self.clr_variables)
            X_tmp = self.clr_transformer.transform(X_tmp)

        # --- StandardScaler (después del CLR) -------------
        if self.standardize_variables is not None:
            self.scaler = StandardScaler()
            self.scaler.fit(X_tmp[self.standardize_variables])

        return self

    def transform(self, X):
        X = X.copy()

        # --- CLR primero ---
        if self.clr_transformer is not None:
            X = self.clr_transformer.transform(X)

        # --- StandardScaler después ---
        if self.scaler is not None:
            X[self.standardize_variables] = self.scaler.transform(
                X[self.standardize_variables]
            )

        return X

    def fit_transform(self, X, y=None):
        return self.fit(X, y).transform(X)
