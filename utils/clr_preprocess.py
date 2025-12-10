import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin


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
