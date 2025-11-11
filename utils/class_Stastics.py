import pandas as pd
import anndata as ad
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib
import re
import seaborn as sns
import statsmodels.api as sm
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from venn import venn
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.cluster import KMeans
# SKBIO IMPORTS
from skbio.stats import subsample_counts
from skbio.diversity import alpha_diversity
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa



class Stastics:

    def __init__(self, data=None, metadata=None,adata=None):

        assert type(adata)!=type(None) or type(data)!=type(None) , "You have to provide some data and metadata (or AnnData) in order to create the object"

        if type(adata)==type(None):

            self.adata = ad.AnnData(data,obs=metadata)
        
        else:

            self.adata = adata

        # Removal of the duplicated columns
        self.adata.obs = self.adata.obs.T.drop_duplicates().T


        # Automatic detection of the column type
        self.adata.obs = self.adata.obs.infer_objects()

        # Create an initial description file
        self.do_description()


        # Initialize a dictionary and a datafame to store continous parways comparisons
        self.dict_comp_1_1 = {"Name":["Normal_Data","Test","P-value","Mean_Difference","Hodges_Lehmann_Estimator","N"]}

        # Initialize a dictionary and a datafame to store categorical parways comparisons
        self.dict_chi_sq = {"Name":["Chi2_statistic","P-value","N"]}
        
        # Initialize a dictionary for anova results
        self.dict_anova= {"Name":["F","P-value","N","Post-Hoc (Tukey)"]}

        
    def __remove_unwanted_chars(self,df,unwanted_chars_dict={},revert_dict={},revert=False):

        df1 = df.copy()
        if revert:
            revert_dict.update({"signopos":"+","_espacio_":" ","signoneg":"-","7barra7":"/","_CORCHETE1_":"[","_CORCHETE2_":"]"})

            for wanted_char in revert_dict.keys():
                df1.columns = df1.columns.str.replace(wanted_char,revert_dict[wanted_char],regex=True)

            return df1

        unwanted_chars_dict.update({"\+":"signopos"," ":"_espacio_","-":"signoneg","/":"7barra7","\[":"_CORCHETE1_","\]":"_CORCHETE2_"})
        for unwanted_char in unwanted_chars_dict.keys():
            df1.columns = df1.columns.str.replace(unwanted_char,unwanted_chars_dict[unwanted_char],regex=True)

        return df1
    


    def do_description(self, name="All", subset="All"):
        """
        Generates a descriptive statistics summary for the observations in the AnnData object.

        Args:
            name (str, optional): Name to be saved with in the adata.uns object. Defaults to "All".
            subset (tuple, optional): Defaults to "All". 
                If provided, creates a subset of the observations by a column of the metadata and a value for this column. E.g., (Group, control)

        Returns:
            None: The resulting DataFrame is stored in the adata.uns attribute of the object.
        """

        # Access the observation data from the AnnData object
        df = self.adata.obs

        # If a subset is specified, filter the DataFrame
        if subset != "All":
            samples = self.adata.obs[self.adata.obs[subset[0]] == subset[1]].index
            df = self.adata.obs.loc[samples]

        # Initialize lists to store summary statistics
        variables, condicion, cuenta, means, medians, normality,iqr = [], [], [], [], [], [], []

        # Iterate over each column in the DataFrame
        for col in df.columns:

            if df[col].dtype == float:
                # Handle numeric columns
                cuentas = np.nan
                condicion.append(cuentas)
                cuen = len(df[col].dropna()) 
                iqr.append(f"{round(df[col].quantile(0.25),2)} - {round(df[col].quantile(0.75),2)}")
                cuenta.append(f"{cuen} ({round(cuen / len(df) * 100, 2)})")
                variables.append(col)
                mean = df[col].mean()
                desvest = df[col].describe()["std"]
                means.append(f"{round(mean, 2)} ± {round(desvest, 2)}")
                medians.append(round(df[col].median(), 2))
                try:
                    normality.append(sp.stats.shapiro(df[col].dropna())[1] > 0.05)
                except :
                    normality.append(None)

            else:
                # Handle categorical columns        
                cuentas = df[col].value_counts()
                for cuent in cuentas.index:
                    condicion.append(cuent)
                    iqr.append(np.nan)
                    cuen = cuentas[cuent]
                    cuenta.append(f"{cuen} ({round(cuen / len(df) * 100, 2)})")
                    variables.append(col)
                    means.append(np.nan)
                    medians.append(np.nan)
                    normality.append(np.nan)

        # Create a DataFrame from the collected summary statistics
        df_res = pd.DataFrame(
            [variables, condicion, cuenta, means, medians, normality,iqr],
            index=["Variable", "Class", "Count (%)", "Mean ± Std_dev", "Median", "Normal Data","IQR"]
        ).T

        # Store the resulting DataFrame in the adata.uns attribute
        self.adata.uns[f"Description_{name}"] = df_res


    def filter(self,var_name,subgroup="DROPNAS",in_place=True):

        """Method that can filter the adata object by a variable of data or metadata

        Args:
            var_name (str): A column of metadata or data to filter by.

            subgroup (str|list, optional): List of values from the selected variable to be kept . Defaults to "DROPNAS", in that case removes all the obs that contains np.nan values for the selected variable.

        Raises:
            ValueError: If subgroup isnt DROPNAS or a list of values.
            KeyError: If the var_name isnt a column of data or metadata.
        """

        if var_name in self.adata.obs.columns: # Filter by metadata

            if subgroup=="DROPNAS": # Filter out the missing values by var_name

                indexes = self.adata.obs.dropna(subset=var_name).index

            elif type(subgroup)==list: # Filter out the  values not in subgroup of the var_name column

                indexes= self.adata.obs_names[self.adata.obs[var_name].isin(subgroup)]

            else:

                print(f"Subgroup must be a list or leave empty, not {type(subgroup)}")

                raise ValueError
            


        elif var_name in self.adata.var_names:  # Filter by variables

            if subgroup=="DROPNAS":# Filter out the missing values by var_name

                indexes = self.adata.to_df().dropna(subset=var_name).index

            elif type(subgroup)==list: # Filter out the  values not in subgroup of the var_name column

                indexes= self.adata.obs_names[self.adata.to_df()[var_name].isin(subgroup)]

            else:

                print(f"Subgroup must be a list or leave empty, not {type(subgroup)}")

        else:

            print(f"The var_name '{var_name}' isn`t a column of data or metadata")
            raise KeyError

        if len(indexes)==0:

            print("Any of the observations meet the desired filtering conditions.\nThe filtering won`t be performed.")

        else:

            if in_place:        
                # Do the filtering
                self.adata = self.adata[indexes,:]

                # Regenerate the description with the new data
                self.do_description()

            else:
                # Return a filtered object
                return Stastics(adata = self.adata[indexes,:])
        
        
    def find_var(self, var, adata_df=None):
        """
        Finds and returns the specified variable from the AnnData object.

        Parameters:
        var (str, list, pd.Series, pd.Index, np.ndarray): The variable(s) to be found.
        adata_df (pd.DataFrame, optional): A DataFrame representation of the AnnData object. If None, it will be created.

        Returns:
        pd.Series or pd.DataFrame: The found variable(s) from the AnnData object.

        Raises:
        KeyError: If the specified variable is not found in the AnnData object.
        AssertionError: If the `var` is not of the accepted types.
        """
        
        # Accepted types for the variable
        accepted_types = (list, pd.Series, pd.Index, np.ndarray)  # Use tuple for type checking
        
        # Check if adata_df is not provided or not a DataFrame
        if type(adata_df) != pd.DataFrame:
            # Convert AnnData object to DataFrame if not provided
            adata_df = self.adata.to_df()  # Store DataFrame once
        
        # If var is an accepted type (list, pd.Series, pd.Index, np.ndarray)
        if isinstance(var, accepted_types):
            # Recursively find variables in the list-like structure
            finded_vars = pd.DataFrame([self.find_var(v, adata_df=adata_df) for v in var]).T
            return finded_vars.infer_objects()
        
        # Ensure var is a string
        assert isinstance(var, str), f"The value of var must be str or one of {accepted_types}"
        
        # Check if the variable is in the observation columns
        if var in self.adata.obs.columns:
            variab = self.adata.obs[var]
            return variab.infer_objects()
        
        # Check if the variable is in the variable names
        if var in self.adata.var_names:
            variab = adata_df[var]
            return variab.infer_objects()
        
        # Raise an error if the variable is not found
        raise KeyError(f"The variable '{var}' can't be found in the anndata object")

    def make_dummy(self,var,in_place=False):
        col = self.find_var(var).dropna()
        dict_dummy = {}
        for v in col.unique():
            if in_place:
                self.adata.obs[f"{var}_{v}"] = col==v
            else:
                dict_dummy[f"{var}_{v}"]=col==v

        if not in_place:
            return pd.DataFrame(dict_dummy)
    

    
    def __check_normality(self ,values:pd.Series,condition=None):

        """ A method that checks the normality of data:

        Args:
            values (pd.Series): the values to check.
            condition (bool, optional): Whether there is a condition (check the normality of both conditions). Defaults to False.


        Returns:
            Boolean, true if normal data.

        """

        if type(condition)!="NoneType": # If a condition is provided

            normal = np.all([sp.stats.shapiro(values[condition.index[condition==i]].astype(float))[1]>0.05  for i in condition.unique()])
            homocedastic = sp.stats.levene(values[condition.index[condition==condition.unique()[0]]].astype(float),
                                           values[condition.index[condition==condition.unique()[1]]].astype(float))[1]>0.05
            # Check the normality of every condition
            return normal & homocedastic
        
        # Check the normality of all values
        return sp.stats.shapiro(values)[1]>0.05

           

        
    def calc_alpha_div(self, vars = "", name="Alpha_div",level=False, index= "shannon"):       
        """
        Calculate alpha diversity for the dataset and store the results in the `obs` attribute.

        Parameters:
        -----------
        vars : str or list, optional
            A list of variable names (features) to include in the calculation. If an empty string 
            (default), all variables in `adata.var.index` are used.

        name : str, optional
            The base name for the new column in `adata.obs` where the results will be stored. 
            Defaults to "Alpha_div". If `level` is specified, the name will include the level 
            and index type.

        level : str or bool, optional
            If specified, aggregates data by the given level (a column in `adata.var`). The 
            aggregated data will be grouped and summed based on this level. Defaults to `False`.

        index : str, optional
            The diversity index to calculate. Must be one of ["shannon", "chao1", "simpson"]. 
            Defaults to "shannon".

        Returns:
        --------
        None
            The calculated alpha diversity values are stored in the `obs` attribute of the 
            `adata` object, under the column name specified by `name`.

        Notes:
        ------
        - If `level` is specified, the data is grouped by the specified level in `adata.var`.
        - The column name in `adata.obs` will include the `level` and `index` if applicable.
        - An assertion ensures that the `index` parameter is valid.
        """

        #assert index in ["shannon", "chao1","simpson","pielou_e"]

        if vars=="":
            vars=self.adata.var.index.tolist()
            df = self.adata.to_df()[vars]
        
        if level:
            df = self.adata.to_df().T
            df[level]=self.adata.var[level]
            df =df.groupby(level).sum()
            df = df.rename({"":"Unclassified"}).T
            name = name+"_level_"+level +  ("_"+index if index != "shannon" else "")

        self.adata.obs[name] = alpha_diversity(index,counts= df,ids=self.adata.obs_names)
        

    def calc_beta_div(self,vars="",in_place=True,level=False):

        if vars=="":
            vars=self.adata.var.index.tolist()
            data = self.adata.to_df()[vars]
        
        if level:
            data = self.adata.to_df().T
            data[level]=self.adata.var[level]
            data =data.groupby(level).sum()
            data = data.rename({"":"Unclassified"}).T

        dm = beta_diversity("braycurtis", data, self.adata.obs_names)

        coord = pcoa(dm)

        df_plo = coord.samples[["PC1",'PC2']]

        distancematrix = pd.DataFrame(dm.data,columns=self.adata.obs_names,index=self.adata.obs_names)

        exp_var=round(coord.proportion_explained*100,2)

        if in_place:

            self.adata.obsm["PCoA"]  = df_plo

            self.adata.obsm["Distance_matrix"] = pd.DataFrame(dm.data,columns=self.adata.obs_names,index=self.adata.obs_names)

            self.adata.uns["PCoA_exp_var"]= exp_var

            print("Saving PCoA and Distance matrix in obsm as: PCoA and Distance_matrix")

        else:

            print("Returning PCoA and Distance matrix")

            return df_plo, distancematrix, exp_var, 

    
    def clr_transformation(self, pseudocount=1e-9, name="CLR_DATA"):
        """
        Apply CLR transformation to the raw counts of an AnnData object and store the result in adata.obsm.
        
        Parameters:
        - adata: AnnData object with raw counts in adata.X.
        - pseudocount: Small value added to avoid issues with zeroes (default is 1e-9).
        - name: The name for storing the CLR-transformed data in adata.obsm (default is "CLR_DATA").
        """
        
        # Convert raw counts (adata.X) to a pandas DataFrame for easier manipulation
        df_ab = self.adata.to_df()  # This assumes raw counts are stored in adata.raw.X
        
        # Convert the DataFrame to numpy array
        data_clr = df_ab.to_numpy(float)
        
        # Add a pseudocount to avoid issues with zeroes
        data_clr += pseudocount
        
        # Compute the geometric mean across rows (samples)
        geometric_mean = np.exp(np.mean(np.log(data_clr), axis=1, keepdims=True))
        
        # Apply CLR transformation
        clr_transformed = np.log(data_clr / geometric_mean)

        clr_transformed = pd.DataFrame(clr_transformed, columns=df_ab.columns, index=df_ab.index)
        
        # Store the CLR transformed data in adata.obsm with the given name
        self.adata.obsm[name] = clr_transformed
        
        # Optionally, return the transformed data as a DataFrame for inspection
        return clr_transformed



    def __hodges_lehmann_estimator(self, group1, group2):
        """
        Calculate the Hodges-Lehmann estimator for two independent samples.
        The Hodges-Lehmann estimator is the median of all pairwise differences 
        between the observations in the two groups.

        Parameters:
        group1 (array-like): First sample
        group2 (array-like): Second sample

        Returns:
        float: Hodges-Lehmann estimator
        """

        group1 = np.asarray(group1)
        group2 = np.asarray(group2)
        
        # Compute all pairwise differences
        diffs = []
        for x in group1:
            for y in group2:
                diffs.append(x - y)
        
        # Return the median of differences
        return np.median(diffs)

    
    def comparisons_1_1(self, target, condition: str = None,name="Comp_1_1"):
        """
        Performs T-test or Mann-Whitney U test comparisons between two variables.

        Args:
            target (str): The target variable to compare.
            condition_name (str): The condition name for grouping the target variable.

        Returns:
            pd.DataFrame: Updated results report of the comparisons.
        """
        
        if type(condition)==type(None):
            assert type(target) in [np.ndarray,pd.Series,list] and len(target)==2,"If condition not provided, the target must be an Iterable wuth len(target)==2"
            # Extract the values of the target variable, dropping missing values
            df_test = self.find_var(var=target).dropna().infer_objects().melt()
            comparison_name = ": ".join([target[0],target[1]])
            N=int(df_test.shape[0]/2)
            target="value"
            condition="variable"
            
        else:
            # Extract the values of the target variable, dropping missing values
            df_test = self.find_var(var=[target,condition]).dropna().infer_objects()
            comparison_name = ": ".join([condition,target])
            N=int(df_test.shape[0])

        values = df_test[target]

        assert values.dtype!= object, "Target must be a continious variable"

        # Extract the corresponding condition values
        condition_df = df_test[condition]

        assert len(condition_df.unique())==2, "There must be two levels to compare"


        # Check normality of the values based on the condition
        normal = self.__check_normality(values, condition_df)

        # Split the values into two Series based on the unique conditions
        sep_values = [values[condition_df.index[condition_df == cond]].astype(float) for cond in condition_df.unique()]


        if normal:
            test = "T-test"
            mean_dif = sep_values[0].mean() - sep_values[1].mean()
            median_dif = np.nan
            pval = sp.stats.ttest_ind(sep_values[0], sep_values[1])[1]
        else:
            test = "Mann-Whitney U"
            pval = sp.stats.mannwhitneyu(sep_values[0], sep_values[1])[1]
            mean_dif = np.nan
            median_dif = self.__hodges_lehmann_estimator(sep_values[0], sep_values[1])


        
        if f"{name}_dict" not in self.adata.uns.keys():

            self.adata.uns[f"{name}_dict"] = {}
            
        self.adata.uns[f"{name}_dict"]["Name"] = self.dict_comp_1_1["Name"]


        self.adata.uns[f"{name}_dict"][comparison_name] = [normal, test, pval, mean_dif, median_dif, N]

        # Create a DataFrame from the comparison results dictionary
        results_df = pd.DataFrame(self.adata.uns[f"{name}_dict"]).set_index("Name").T

                      
        if results_df.shape[0]>1:

            results_df["FDR"] = sp.stats.false_discovery_control(results_df["P-value"].tolist(), method="bh")

            results_df["Significant"] = results_df["FDR"]<0.05
        

        # Update the instance attribute with the results DataFrame
        self.adata.uns[name] = results_df

        return results_df


    def __correlations(self, df1, df2, name):
        """
        Calculate the correlation matrix between two dataframes and store the results.

        Parameters:
        df1 (pd.DataFrame): First dataframe with the variables of interest.
        df2 (pd.DataFrame): Second dataframe with the variables of interest.
        name (str): Name to prefix the stored results in the annotated data (adata) object.

        Returns:
        None: The function stores the correlation matrix, p-value matrix, and sample size matrix 
            in the `adata.uns` dictionary with keys based on the provided name.
        """
        
        # Initialize the correlation matrix
        corr_matrix = pd.DataFrame(index=df1.columns, columns=df2.columns)
        pval_matrix = pd.DataFrame(index=df1.columns, columns=df2.columns)
        n_matrix = pd.DataFrame(index=df1.columns, columns=df2.columns)
        
        # Calculate the correlations
        for col1 in df1.columns:
            for col2 in df2.columns:

                # Filter non-null values
                valid_idx = df1[col1].notna() & df2[col2].notna()
                x = df1.loc[valid_idx, col1]
                y = df2.loc[valid_idx, col2]

                method = "spearman"

                # Check the conditions required to do Pearson correlation (continuous and normal data)
                if (x.dtype == float) & (y.dtype == float):
                    if (sp.stats.shapiro(x)[1] > 0.05) & (sp.stats.shapiro(y)[1] > 0.05): 
                        method = "pearson"

                # Perform the corresponding correlation
                if method == "spearman": 
                    corr, pval = sp.stats.spearmanr(x, y)
                else:
                    corr, pval = sp.stats.pearsonr(x, y)

                # Store the results in the matrices
                corr_matrix.at[col1, col2] = corr
                pval_matrix.at[col1, col2] = pval
                n_matrix.at[col1, col2] = len(x)

        # Store the matrices in the annotated data (adata) object
        self.adata.uns[f"{name}_Corr"] = corr_matrix
        self.adata.uns[f"{name}_Corr_pval"] = pval_matrix
        self.adata.uns[f"{name}_Corr_N"] = n_matrix



    
    def do_correlations(self, Variables_A, Variables_B, name, return_df = False):
        """
        Uses the __correlation() method from this class to create the correlation matrix 
        from the two groups of variables and generates a report from it.

        Args:
            Variables_A (list): List of variables to correlate from the first group.
            Variables_B (list): List of variables to correlate from the second group.
            name (str): Name to prefix the stored results and generated report.
            return_df (bool, optional): Whether to return the resultant df or not.

        Returns:
            pd.DataFrame: DataFrame containing the generated correlation report.
        """

        # Define the dataframes with the selected columns
        df1 = self.find_var(Variables_A)
        if len(df1.shape)==1:
            df1=df1.to_frame()

        df2 = self.find_var(Variables_B)
        if len(df2.shape)==1:
            df2=df2.to_frame()
        
        # Create the correlation matrix
        self.__correlations(df1, df2, name=name)
        
        # Generate the report from the correlation matrix
        results_df = self.generate_corr_report(name=name)

        if return_df:
            return results_df




    def order_comparisons(self,index):
        """
        Args:
            index (pd.Index): index of the results DataFrame.

        Returns:
            Values of the index strings oredered.
        """
        parts = index.split(" vs ")
        parts.sort()

        return " vs ".join(parts)


    def generate_corr_report(self, name, select_vars=None):
        """
        Generates a report from a correlation analysis previously run.

        Parameters:
            name (str): The name of the correlation analysis.
            select_vars (list, optional): List of selected variables to include in the report.

        Returns:
            pd.DataFrame: DataFrame containing the correlation report.
        """
        
        matrix_dict = self.adata.uns

        var1, var2, corr, pval, num = [], [], [], [], []

        df_bool = matrix_dict[f"{name}_Corr"] < 0.999  # DataFrame indicating the self-correlation variables

        if select_vars is not None:           
            select_vars = set(select_vars)
            assert len(select_vars.intersection(set(df_bool.columns.tolist() + df_bool.index.tolist()))) == len(select_vars), "One or more variables not found in the correlation matrix"
            
            cols = list(select_vars.intersection(df_bool.columns))
            if len(cols) > 0:
                df_bool = df_bool[cols]
            
            rows = list(select_vars.intersection(df_bool.index))
            if len(rows) > 0:
                df_bool = df_bool.loc[rows]

        for col in df_bool.columns:
            for term in df_bool.index[df_bool[col]].tolist():
                var1.append(col)
                var2.append(term)
                corr.append(matrix_dict[f"{name}_Corr"].loc[term, col])
                pval.append(matrix_dict[f"{name}_Corr_pval"].loc[term, col])
                num.append(matrix_dict[f"{name}_Corr_N"].loc[term, col])

        # Convert the lists to a DataFrame for better visualization
        results_df = pd.DataFrame({
            'Variable_1': var1,
            'Variable_2': var2,
            'Correlation': corr,
            'P-value': pval,
            'N': num
        })

        results_df["Correlated Variables"] = results_df.Variable_1 + " vs " + results_df.Variable_2
        results_df = results_df.dropna().set_index("Correlated Variables")

        # Apply the function to the indexes
        results_df.index = results_df.index.map(self.order_comparisons)

        # Eliminate duplicated indexes
        results_df = results_df[~results_df.index.duplicated(keep='first')]

        # Perform false discovery rate control
        results_df["FDR"] = sp.stats.false_discovery_control(results_df["P-value"], method="bh")
        results_df["Significative"] = results_df["FDR"] < 0.05

        # If select_vars is specified, handle the selection logic
        if select_vars is not None:
            prefix = 1 + max([int(i.split("_")[0]) if i.split("_")[0].isdigit() else 0 for i in self.adata.uns.keys()])
            name = f"{prefix}_Selection"

            # If the report has not been generated yet
            if np.all(np.array([self.adata.uns[result].index.tolist() != results_df.index.tolist() for result in self.adata.uns.keys()])):
                self.adata.uns[f"{name}_Corr_report"] = results_df
        else:
            self.adata.uns[name] = results_df

        return results_df
    

    def anova(self,target,condition,name="ANOVA_RESULTS"):
        
        df_test = self.find_var([target,condition]).dropna()

        df_test = self.__remove_unwanted_chars(df_test)

        N = df_test.shape[0]
        
        model= sm.formula.ols(formula=f"{df_test.columns[0]} ~ {df_test.columns[1]}", data = df_test).fit()

        # Realizar ANOVA
        anova_table = sm.stats.anova_lm(model, typ=2)
        
       
        tukey = None

        # If ANOVA is significant, perform Tukey's HSD

        resT = ""
        if anova_table['PR(>F)'].iloc[0] < 0.05:

            # Create and format the post-hoc results dataframe
            tukey =  pd.DataFrame(sm.stats.multicomp.pairwise_tukeyhsd(endog=df_test[df_test.columns[0]], 
                                                                       groups=df_test[df_test.columns[1]], alpha=0.05).summary().data)
            tukey.columns = tukey.iloc[0]
            tukey = tukey[1:].reset_index(drop=True)
            tukey[tukey["reject"].astype(bool)]

            for i in tukey[tukey["reject"]].iterrows():
                resT += f"{i[1]["group1"]} vs {i[1]["group2"]}; pval:{i[1]["p-adj"]}; meandiff:{i[1]["meandiff"]}||"

        

        if f"{name}_dict" not in self.adata.uns.keys():

            self.adata.uns[f"{name}_dict"] = {}
            
        self.adata.uns[f"{name}_dict"]["Name"] = self.dict_anova["Name"]
        self.adata.uns[f"{name}_dict"][": ".join([target,condition])] = anova_table.loc[df_test.columns[-1]][["F","PR(>F)"]].tolist() + [N,resT]

        # Create a DataFrame from the comparison results dictionary
        results_df = pd.DataFrame(self.adata.uns[f"{name}_dict"]).set_index("Name").T

                      
        if results_df.shape[0]>1:

            results_df["FDR"] = sp.stats.false_discovery_control(results_df["P-value"].tolist(), method="bh")

            results_df["Significant"] = results_df["FDR"]<0.05
        
       
       # Save the data in adata.uns
        self.adata.uns[name] = results_df 


        return results_df,tukey


    

    def chi_sq(self, col: str, expected_proportions=[0.5, 0.5]):
        """
        Perform a chi-squared test with the values of a metadata column.

        Parameters:
            col (str): The name of the metadata column to perform the chi-squared test on, or list of var names to compare, 
            the first is the reference.
            expected_proportions (list, optional): The expected proportions for each category. Default is [0.5, 0.5].

        Returns:
            p_value (float): Value P of the chi-squared test.
            None: The results of the chi-squared test are stored in the `adata.uns` dictionary and printed.

        """

        # If two variables are going to be compared
        if type(col)!=str:

            name = "-".join(col)

            chi_df = self.find_var(col).dropna()
            df_contingency = pd.crosstab(chi_df[col[0]],chi_df[col[1]])
            chi2, p_value, dof, expected  = sp.stats.chi2_contingency(df_contingency)
            N = chi_df.shape[0]

        else: 
            name = col
            # Extract the number of observed counts
            data = self.find_var(col).dropna()
            observed_counts = [i for i in data.value_counts()] 
            total_count = sum(observed_counts) 
            # Calculate the expected counts
            expected_counts = [total_count * p for p in expected_proportions] 

            # Perform the chi-squared test
            chi2, p_value = sp.stats.chisquare(f_obs=observed_counts, f_exp=expected_counts)

            N = len(data)


        # Store the comparison results in a dictionary
        self.dict_chi_sq[name] =[chi2, p_value, N]

        # Create a DataFrame from the comparison results dictionary
        results_df = pd.DataFrame(self.dict_chi_sq).set_index("Name").T
        
        # Evaluate fdr
        results_df["FDR"] = sp.stats.false_discovery_control(results_df["P-value"].tolist(), method="bh")
        results_df["Significant"] = results_df["FDR"] < 0.05
        
        # Save the data in adata.uns
        self.adata.uns["CHI_SQ_TABLE"] = results_df 

        return results_df
    

    def clustering(self,variables:list=[], name = "K_means",kind="k-means",clusters=2,from_obsm=False):

        if not from_obsm:
            df = self.find_var(variables)

        else:
            df = self.adata.obsm[from_obsm]

        if kind=="k-means":
            # Predict the clusters
            self.adata.obs[name] = KMeans(n_clusters=clusters,random_state=123).fit_predict(df)
        else:
            print("Option not implemented.")

    def generar_diccionario_colores(self,valores):
        # Genera un colormap que cubre una cantidad de colores igual a la longitud de la lista
        colores = plt.cm.tab10(range(len(valores)))  # Usa colormap 'tab10', que tiene 10 colores distintos
        
        # Crear el diccionario con cada valor de la lista como clave y su color asignado como valor
        diccionario_colores = {valor: colores[i] for i, valor in enumerate(valores)}
        return diccionario_colores


        
    def __get_cov_ellipse(self,cov, nstd):
        """
        Create an ellipse based on a covariance matrix and center point.
        """
        eigvals, eigvecs = np.linalg.eigh(cov)
        order = eigvals.argsort()[::-1]
        eigvals, eigvecs = eigvals[order], eigvecs[:, order]
        vx, vy = eigvecs[:,0]
        theta = np.degrees(np.arctan2(vy, vx))
        width, height = 2 * nstd * np.sqrt(eigvals)

        return width, height, theta
    
    
    def plot_pcoa(self,obsm_key,condition_name,save=False,show=True,ax=None,exp_var = None, palette=None):
        # Charge the data to plot
        df_plot = self.adata.obsm[obsm_key][["PC1","PC2"]]

        # Create a column with the condition info
        condition = self.find_var(condition_name).dropna().astype(str)

        df_plot[condition_name] = df_plot.index.map(dict(zip(condition.index, condition)))

        df_plot = df_plot.dropna()
        
        centroids = df_plot.groupby(condition_name).mean().reset_index()

        #centroids[condition_name] =  "CENTER "+centroids[condition_name]
        centroids.set_index(condition_name,inplace=True)

        if palette is None:
            # Obtener el color de las elipses
            palette = self.generar_diccionario_colores(self.find_var(condition_name).unique())


        if ax is None:  # If no axis is provided, create a new figure
            fig, ax = plt.subplots(1, 1, figsize=(10, 6))

        sns.scatterplot(ax=ax, data=df_plot, x="PC1", y="PC2", hue=condition_name, s=40,palette=palette)


        #color_dict = {label: palette[idx] for idx, label in enumerate(df_plot[condition_name].unique())}

            # Añadir las elipses
        for label, group in df_plot.groupby(condition_name):
            cov = np.cov(group[['PC1', 'PC2']].T)
            center = group[['PC1', 'PC2']].mean()
            width, height, angle = self.__get_cov_ellipse(cov, 2)

            ellipse = matplotlib.patches.Ellipse(xy=(center["PC1"], center["PC2"]), width=width, 
                            height=height, angle=angle, edgecolor=palette[label], 
                            facecolor='none', linestyle='--')
            ax.add_patch(ellipse)

             # Asegúrate de graficar el centroide en el eje correcto
            ax.scatter(centroids["PC1"].loc[label], centroids["PC2"].loc[label], color=palette[label], s=200, marker="s")
            
  
        # Set the texts on the plot
        if exp_var is None:
            ax.set_xlabel("PCo1")
            ax.set_ylabel("PCo2")    
        else:

            ax.set_xlabel(f"PCo1 (Explained Variance {exp_var["PC1"]}%)")
            ax.set_ylabel(f"PCo2 (Explained Variance {exp_var["PC2"]}%)")   

        plt.legend( bbox_to_anchor=(1, 1),fontsize=9)
        plt.tight_layout()
        if save:
            plt.savefig(save,bbox_inches="tight")
        
        if show:
            plt.show()

        if not ax is None:
            return ax

    def plot_merged_abundances(self, condition, taxon, palette=None, ax = None,
                        level="Genus",violin=False,y_log=False,save=False,
                        show=True):
        """
        Creates a plot to compare the abundances of a specified taxon at a particular 
        taxonomic level (e.g., Genus, Family) across groups defined by a condition.

        Args:
            condition (str): The metadata column name that defines the groups to compare.
            taxon (str): The name of the taxon to plot (e.g., a genus or family).
            palette (dict): A dictionary mapping group names to colors for the plot.
            ax (matplotlib.axes._axes.Axes, optional): An existing matplotlib Axes object 
                to plot on. If None, a new figure and Axes will be created. Defaults to None.
            level (str, optional): The taxonomic level to consider (e.g., "Genus", "Family"). 
                Defaults to "Genus".
            violin (bool, optional): Whether to create a violin plot instead of a box plot. 
                Defaults to False.
            y_log (bool, optional): Whether to use a logarithmic scale for the y-axis. 
                Defaults to False.
            save (bool, optional): Whether to save the plot to a file. Defaults to False.
            show (bool, optional): Whether to display the plot. Defaults to True.

        Returns:
            matplotlib.axes._axes.Axes: The Axes object containing the plot, which allows 
                for further customization or saving.
        """

        targetvars = self.adata.var[level].dropna()[self.adata.var[level].dropna()==taxon].index
        targetvar = self.find_var(targetvars).sum(axis=1)
        df_plot = pd.DataFrame([targetvar,self.find_var(condition)],index=[taxon,condition]).T

        # Create an axe for plotting
        if ax is None:
            fig,ax = plt.subplots(1,figsize=(5,5))

        if violin:
            sns.violinplot(df_plot,x=condition,y=taxon,palette=palette,ax=ax)

        else:
            sns.boxplot(df_plot,x=condition,y=taxon,fliersize=0,palette=palette,ax = ax)

        sns.stripplot(df_plot,x=condition,y=taxon,color="black",alpha=0.7,ax=ax)

        if y_log:
            ax.set_yscale("log")

        if save:
            plt.savefig(save,bbox_inches="tight")
        
        # Show the result
        if show:
            plt.show()

        # Return the ax object
        else:
            return ax 

    def plot_differences(self, condition, vars, kind="Box", ylab="", xlab=" ", tick_label_names= [],title = "",obsm_name=False,
                         ylog=False,show_n=True, save=False, show="Show",theme=False,palette="deep",figsize=(10, 6)):
        """
        Method to create violin or boxplots comparing different states of the data.

        Args:
            condition (str): The condition to compare the variables by.
            vars (list): Name of the variables to compare (list, np.array, or pd.Series).
            kind (str, optional): "Box" for boxplot, "Violin" for violin plot, "Bar" for barplot. Defaults to "Box".
            ylab (str, optional): Label of the y-axis. Defaults to "".
            xlab (str, optional): The name of the variables, appears in the x label. Defaults to " ".
            theme (str,optional): Set a sns theme, options are ['paper', 'notebook', 'talk', 'poster']. Defaults to False.
            palette (str|list): Set the colours of the plot (hue differences). Defoults to 'deep'.
            show_n (bool, optional): If False it hides the number of samples being plotted from the tick label name. Defaults to True.
            tick_label_names (list, optional): Change the name of the tick labels shown.
            save (bool, optional): Path to save the figure. Defaults to False.
            show (str, optional): Whether to show the plot, return it for user edition or close it directly. Defaults to "Show".
            obsm_name = bool, if a str is provided uses this as the reference for plotting

        Raises:
            AttributeError: If an unsupported kind of plot is requested.
        """
        if theme:
            sns.set_theme(theme)

        if obsm_name:

            adata_df = self.adata.obsm[obsm_name]

        else:
            adata_df = None


        if type(vars)!=str:
            df_var = self.find_var(var=list(vars)+[condition],adata_df=adata_df)
        else:
            df_var = self.find_var(var=[vars,condition],adata_df=adata_df)
            
        df_var = df_var.melt(id_vars=condition, value_name="value", var_name=xlab.upper())

        # Create the plot
        fig = plt.figure(figsize=figsize)

        if kind == "Violin":
            sns.violinplot(x=xlab, y="value", data=df_var, inner='quart', hue=condition,palette=palette, split=False)

        elif kind == "Box":
            sns.boxplot(x=xlab, y="value", data=df_var, hue=condition,palette=palette, fliersize=0)
            if not ylog:
                # Enhance aesthetics
                plt.grid(True, which='both', linestyle='--', linewidth=0.5)  # Add a subtle grid

        elif kind== "Bar":
            sns.barplot(x=xlab, y="value", data=df_var,
                        err_kws={'linewidth': 1.6},capsize=0.1, edgecolor='black', palette=palette,hue=condition)
        
        else:
            print("Option not implemented yet")
            raise AttributeError
        
        ax = sns.stripplot(x=xlab, y="value", data=df_var, hue=condition, dodge=True, jitter=True, palette='dark:k', alpha=0.5, legend=False)


        if ylog:
            plt.yscale("log")

        plt.title(title,fontweight=600)
          

        # Add labels
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.xticks(rotation=25, ha='right')

        # Place the legend
        plt.legend(title=condition, loc='best')

        sns.despine()

        if show_n:
            # Extraer los tick labels actuales del eje x y generar unos nuevos con el valor de la n
            x_tick_labels = [f"{item.get_text()} (n={df_var.dropna()[" "].value_counts()[item.get_text()]}) " for item in ax.get_xticklabels()]
        else:
            x_tick_labels = [item.get_text() for item in ax.get_xticklabels()]

        # Change the tick label names
        if type(tick_label_names)==list and len(tick_label_names)>0:
            x_tick_labels = [re.sub(".*\(",f"{tick_label_names[i]} (",item) for i, item in enumerate(x_tick_labels)]
            if not show_n:
                x_tick_labels = [tick_label_names[i] for i, item in enumerate(x_tick_labels)]
                
        # Aplicar los nuevos tick labels al eje x
        ax.set_xticks([i[0] for i in enumerate(x_tick_labels)],labels = x_tick_labels )

        sns.reset_orig()
        if save:
            plt.savefig(save, bbox_inches="tight")

        if show=="Show":
            plt.show()

        elif show == "Edit":
            return fig,ax # Returns the figure so the user can edit it.
        
        else:
            plt.close()
            
        return df_var
        

    def plot_pca(self,vars:list,hue:str,title="",figsize=(8,6),save=False,show=True, ret_data = False):
        """Calculates the pca of a set of variables and plot the 2 first components, colouring by a condition.

        Args:
            vars (list): list of variables to use.
            hue (str): Defaults to ""
            title (str, optional): title of the plot. Defaults to "".
            figsize (tuple, optional): size of the plot. Defaults to (8,6).
            save (bool, optional): Path where the plot is saved or not to save it. Defaults to False.
            show (bool, optional): Show the plot or not. Defaults to True.
        Returns:
            If ret_data==True returns a dataframe with the 2 PC and the hue.
        """
        # Prepare the data
        dat = self.find_var(vars).dropna()  # Extract the data in a dataframe
        data = StandardScaler().fit_transform(dat) # Scale the data
        pca = PCA().fit(data)# Fit the PCA  
        exp_var = pca.explained_variance_ratio_
        data = pd.DataFrame(data=pca.transform(data)[:,:2],columns=["PC1","PC2"] ,index=dat.index)
        data[hue]=self.find_var(hue).loc[dat.index] # Add the condition column
        
        # Crea una figura y un objeto de ejes
        fig, ax = plt.subplots(figsize=figsize)
        sns.scatterplot(data=data,x="PC1",y="PC2",hue=hue,ax = ax)
        plt.title(title,fontweight=800)
        plt.xlabel(f"PC1 (Exp variance: {round(exp_var[0]*100,2)} %)")
        plt.ylabel(f"PC2 (Exp variance: {round(exp_var[1]*100,2)} %)")
        if save:
            plt.savefig(save, bbox_inches="tight")

        if show:
            plt.show()
        else:
            plt.close()

        if ret_data:
            return data
        


    def plot_venn(self,columns:list,make_dummys=True,figsize= (20,20),show=True,save=False):

        """Plot a venn diagram of the samples intersection in 3 variables (booleans)

        Args:
            columns (list): a list of columns names, the columns must have bool values.
        """

        df = self.find_var(columns) 

        if make_dummys:
            df = pd.concat([self.make_dummy(var = col) for col in columns],axis=1).dropna()
            sub_plts = df.columns[df.columns.str.contains(columns[0])].tolist()
            sub_plts_opt = df.columns[df.columns.str.contains(columns[1])].tolist()

            fig,axs = plt.subplots(ncols=len(sub_plts),nrows=1,figsize=figsize)

            for ax,plot in enumerate(sub_plts):
                df_plot = df[[plot]+sub_plts_opt]
                venn({i:set(df_plot[df_plot[i]].index) for i in df_plot.columns},ax=axs[ax])

        else:
            venn({i:set(df[df[i]].index) for i in df.columns})

        if save:
            plt.savefig(save,bbox_inches="tight")
        if show:
            plt.show()

    def roc_curve(self,target:str,label:str,title=False,save=False,show=True):

        """Finds the optimal threshold and 

        Args:
            target (str): the continious value.
            label (str): the categorycal value.
            save (bool, optional): Path where the plot is saved or not to save it. Defaults to False.
            show (bool, optional): Show the plot or not. Defaults to True.
        Returns:
            optimal_threshold (np.float): value of the optimal threshold. 
        """

        df_roc = self.find_var(var=[target,label]).dropna() # Import the data for the roc curve
        df_roc1= df_roc.copy() # Backup_copy

        dic = {val:i for i,val in enumerate(df_roc[label].unique())} # codyfy the labels into numeric values

        df_roc[label]= df_roc[label].map(dic) # convert into a dummy variable

        auc = roc_auc_score(df_roc[label],df_roc[target])

        if auc<0.5: # Checks that the number 1 is placed in the positives
            dic = {val:i for i,val in enumerate(df_roc1[label].sort_values(ascending=False).unique())} # codyfy the labels into numeric values
            df_roc[label]= df_roc1[label].map(dic) # convert into a dummy variable
            auc = roc_auc_score(df_roc[label],df_roc[target])

        fpr, tpr, thresholds = roc_curve(df_roc[label],df_roc[target]) # Calculate the roc curve parameters
        
        # Encontrar el punto de corte óptimo (Youden's J)
        youden_index = tpr - fpr
        optimal_idx = np.argmax(youden_index)
        optimal_threshold = thresholds[optimal_idx]


        # Graficar la curva ROC
        plt.figure()
        plt.plot(fpr, tpr, color='darkorange', lw=2, label='Curva ROC (AUC = %0.2f)' % auc)
        plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')

        plt.scatter(fpr[optimal_idx], tpr[optimal_idx], color='red', label='Optimal Threshold  (%.2f)' % optimal_threshold)

        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')

        if title:
            plt.title(title,fontweight=800)
        else:  
            plt.title(f'Curva ROC {target}',fontweight=800)

        plt.legend(loc="lower right")
        
        if save:
                plt.savefig(save, bbox_inches="tight")

        if show:
            plt.show()
        else:
            plt.close()
        
        return optimal_threshold

            

    def plot_correlation(self, var_x, var_y,hue=False, save=False, show="Show",figsize=(5, 1.5),
                         theme=False,palette="deep"):
        """
        Plot the correlation between two variables using a scatter plot.

        Args:
            var_x (str): The name of the first variable (x axis).
            var_y (str): The name of the second variable (y axis).
            hue (str, optional): The condition to split the data by.  
            palette (str or list, optional): The colors in the plot. Defaults to "deep".
            theme (str,optional): Set a sns theme, options are ['paper', 'notebook', 'talk', 'poster']. Defaults to False.            
            save (bool or str, optional): Path to save the figure. Defaults to False.
            show (str, optional): Whether to show the plot, return it ("Edit") or close it directly. Defaults to "Show".
        """
        if theme:
            sns.set_theme(theme)

        if hue:
            df_plot = self.find_var([var_x, var_y,hue]).dropna()
            ax = sns.lmplot(x=var_x, y=var_y, hue=hue, height=figsize[0],aspect=figsize[1],
                        data=df_plot,palette=palette,legend_out=False)

        else:
            df_plot = self.find_var([var_x, var_y]).dropna()
            ax = sns.lmplot(x=var_x, y=var_y,data=df_plot,
                       palette=palette,legend_out=False,height=figsize[0],aspect=figsize[1])

        sns.despine()

        plt.title(f"Correlation: {var_y} vs {var_x} (n={df_plot.shape[0]})", fontweight=800)

        sns.reset_orig()

        if save:
            plt.savefig(save, bbox_inches="tight")

        if show=="Show":
            plt.show()

        elif show == "Edit":
            return ax # Returns the figure so the user can edit it.
        
        else:
            plt.close()