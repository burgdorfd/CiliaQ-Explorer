import os
import pandas as pd
import tkinter as tk
from tkinter import *
from tkinter import filedialog, ttk
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import shapiro
from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multitest import multipletests
import scikit_posthocs as sp
import sklearn
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import itertools
import pingouin as pg
import umap
import re

class Statistics():
    def __init__(self, groups_param, data_param, selected_measurements_param):
        self.all_groups = groups_param
        self.all_data = data_param
        self.tukey_df = pd.DataFrame()
        self.dunn_df = pd.DataFrame()
        self.column_names = selected_measurements_param
        self.measurement_selection = self.column_names
        #self.measurement_selection = self.column_names.remove(["Replicate", "Group"])
        # prepare the dataframe
        for measurement in self.measurement_selection:
            if (measurement not in ["Replicate", "Group"]):
                try:
                    if(not self.all_data[measurement].empty):
                        self.all_data[measurement] = pd.to_numeric(self.all_data[measurement], errors="coerce") 
                except KeyError:
                    pass

    def get_dunn_results(self):
        return self.dunn_df
    
    def get_tukey_results(self):
        return self.tukey_df
    
    def hypothesis_testing(self, data_param, alpha, measurement_name):
        #data preprocessing
        conditions = data_param.columns.tolist()
        measurement = measurement_name
        data = data_param.apply(pd.to_numeric, errors='coerce')
        n_conditions = len(conditions)

        #2 .test for normal distributed data
        result_norm_dist = {}
        for condition in conditions:
            data_norm_dist = data[condition].dropna()

            if (len(data_norm_dist) >= 3):
                stat, p_shapiro = stats.shapiro(data_norm_dist)
                result_norm_dist[condition] = p_shapiro
            else: 
                print(f"The group {condition} contains not enough cilia to run hypothesis testing. {condition} was excluded from hypothesis testing")
                #result_norm_dist[condition] = np.nan
        
        norm_dist_summary = all(value > alpha for key, value in result_norm_dist.items())
        summary = {"measurement" : measurement,
                   "Norm. distribution": norm_dist_summary} 

        if n_conditions < 2: 
            print("Insert data from more conditions please.")
            return None

        if (n_conditions == 2):
            cond1 = data.iloc[:, 0].dropna()
            cond2 = data.iloc[:, 1].dropna()
            result_homogen = {}
            #3. test for homogeneity of variance
            if (len(cond1) >= 3 and len(cond2) >= 3):
                stat, p_levene = stats.levene(cond1, cond2)
                result_homogen.update({f"condition": p_levene})
            else:
                print(f"The group {condition} contains not enough cilia to run hypothesis testing. {condition} was excluded from hypothesis testing ---")

            homogen_dist_summary = all(value > alpha for key, value in result_homogen.items())
            summary.update({"Homogeneity of Variance": homogen_dist_summary})

            #perform hypothesis testing    
            if((summary["Norm. distribution"] == True) and (summary["Homogeneity of Variance"] == True)):
                stat, p = stats.ttest_ind(cond1, cond2)
                performed_test = "t-Test"
            elif((summary["Norm. distribution"] == False) and (summary["Homogeneity of Variance"] == False) or 
                 (summary["Norm. distribution"] == False) and (summary["Homogeneity of Variance"] == True)):
                stat, p = stats.mannwhitneyu(cond1, cond2)
                performed_test = "Mann-Witney-U Test"
            elif((summary["Norm. distribution"] == True) and (summary["Homogeneity of Variance"] == False)):
                stat, p = stats.ttest_ind(cond1, cond2, equal_var=False)
                performed_test = "welch's ttest"
            else:
                print(f"The group {condition} did not neet the quality standards for analysis. {condition} was excluded from hypothesis testing ---")

            pairwise = {
                f"{cond1.name} vs {cond2.name}": [self.get_asterisk(p), p]}
            
            summary.update({"hypothesis test": performed_test, 
                            "test-stat": stat,
                            "p-value": p,
                            "compared_groups": f"{cond1.name}-{cond2.name}",
                            "asterisks" : self.get_asterisk(p),
                            "pairwise": pairwise
                            })
            return summary
            
        elif(n_conditions > 2):
            data_col_list = [data[col].dropna() for col in conditions]

            result_homogen = {}

            if all(len(col_data) >= 3 for col_data in data_col_list):
                stat, p_levene = stats.levene(*data_col_list)
                result_homogen.update({"All conditions": p_levene})
            else:
                print(f"At least one group contains not enough cilia to run hypothesis testing. Skipping Levene's test.")
                result_homogen.update({"All conditions": np.nan})



            #3. test for Homogeneity of variance
            #if (len(cond1) >= 3 and len(cond2) >= 3):
            #    stat, p_levene = stats.levene(cond1, cond2)
            #    result_homogen.update({f"condition": p_levene})
            #else:
            #    print(f"The group {condition} contains not enough cilia to run hypothesis testing. {condition} was excluded from hypothesis testing ---")

            homogen_dist_summary = all(value > alpha for key, value in result_homogen.items())
            summary.update({"Homogeneity of Variance": homogen_dist_summary})

            if((summary["Norm. distribution"] == True) & (summary["Homogeneity of Variance"] == True)):
                stat, p = stats.f_oneway(*data_col_list)
                performed_test = "One-way ANOVA"
            elif((summary["Norm. distribution"] == False) & (summary["Homogeneity of Variance"] == True) or 
                 (summary["Norm. distribution"] == False) & (summary["Homogeneity of Variance"] == False)):
                stat, p = stats.kruskal(*data_col_list)
                performed_test = "Kruskal-Test"
            elif((summary["Norm. distribution"] == True) & (summary["Homogeneity of Variance"] == False)):
                melted_data = data.melt(var_name="group", value_name="value")
                welch_res = pg.welch_anova(dv="value", between="group", data=melted_data)
                stat, p = welch_res["F"].values[0], welch_res["p-unc"].values[0]
                performed_test = "Welch's ANOVA"

            summary.update({"hypothesis test": performed_test, 
                            "test-stat": stat,
                            "p-value": p,
                            "compared_groups": "general",
                            "asterisks" : self.get_asterisk(p)
                            })
            
            #pairwise comparisons
            pairwise_results = {}
            if((performed_test == "One-way ANOVA") & (p < alpha)):
                melted_data = data.melt(var_name = "group", value_name = "value")
                tukey_results = pairwise_tukeyhsd(melted_data["value"], melted_data["group"], alpha=alpha)
                tukey_summary = tukey_results.summary().data[1:]
                for i in range(len(tukey_results.pvalues)): 
                    for comparison in tukey_summary:
                        cond1 = comparison[0]
                        cond2 = comparison[1]            
                        p_value = float(comparison[3])           
                        asterisk = self.get_asterisk(p_value)
                        pair = f"{cond1} vs {cond2}"
                        pairwise_results[pair] = [asterisk, p_value]
                summary.update({"pairwise": pairwise_results})
                return summary

            elif ((performed_test == "Kruskal-Test") and (p < alpha)):
                melted_data = data.melt(var_name = "group", value_name = "value")              
                dunn_results = sp.posthoc_dunn(melted_data, val_col="value", group_col="group", p_adjust='bonferroni')
                conditions = dunn_results.columns.tolist()
                for i, cond1 in enumerate(conditions):
                    for cond2 in conditions[:i]:
                        p_value = float(dunn_results.loc[cond1, cond2])
                        asterisk = self.get_asterisk(p_value)
                        pair = f"{cond1} vs {cond2}"
                        pairwise_results[pair] = [asterisk, p_value]
                summary.update({"pairwise": pairwise_results})
                return summary

            elif((performed_test == "Welch's ANOVA") and (p < alpha)):
                games_howell_results = pg.pairwise_gameshowell(data = melted_data, dv = "value", between = "group")         
                for index, row in games_howell_results.iterrows():
                    cond1 = row["A"]
                    cond2 = row["B"]            
                    p_value = float(row["pval"])           
                    asterisk = self.get_asterisk(p_value)
                    pair = f"{cond1} vs {cond2}"
                    pairwise_results[pair] = [asterisk, p_value] 
                summary.update({"pairwise": pairwise_results})
                return summary
        
    def get_asterisk(self, p_value):
        if(p_value <= 0.001):
            return "***"
        elif(p_value <= 0.01):
            return "**"
        elif(p_value <= 0.05):
            return "*"
        else:
            return "ns"

    def check_norm_dist(self):
        pvalues = {}
        for column in self.data:
            p_val = shapiro(self.data[column])[1]
            pvalues[column] = p_val
        if all(value >= 0.05 for value in pvalues.values()):
            return True
        else:
            return False
        
    def calculate_des_Statistics(self):
        des_stat = pd.DataFrame(columns = ["Measurement", "Group", "Mean", "Median", "Std", "Min", "Max", "Count"])
        for measurement in self.measurement_selection:
            if(measurement not in ["Replicate", "Group"]):
                stats_temp = self.all_data.groupby("Group")[measurement].agg(["mean", "median", "std", "min", "max", "count"])
                stats_temp.reset_index(inplace=True)
                stats_temp.columns = stats_temp.columns.str.capitalize()
                stats_temp.insert(0, "Measurement", measurement)
                stats_temp = stats_temp.round(2)
                des_stat = pd.concat([des_stat, stats_temp], ignore_index = True)
        des_stat.sort_values(["Measurement", "Group"], inplace=True)
        return des_stat
    
    def perform_UMAP(self):
        data = self.all_data
        col_to_drp = ['ID', 
                          'Replicate', 
                          "x center [micron]",
                          "y center [micron]",
                          "z center [micron]",
                          "Intensity threshold A",
                          "Intensity threshold B",
                          "Intensity threshold Basal Stain",
                          "Volume [micron^3]",
                          "Surface [micron^2]",
                          #"orientation vector x [micron] (vector from first to last skeleton point)",
                          #"orientation vector y [micron] (vector from first to last skeleton point)",
                          #"orientation vector z [micron] (vector from first to last skeleton point)",
                          ]
        data = data.drop([col for col in col_to_drp if col in data.columns], axis=1)
        data.dropna(how='all', axis=1, inplace=True)
        data.dropna(axis=0, how="any", inplace=True)
        group = data['Group']
        data = data.drop('Group', axis=1)  
        data.replace("�", 0, inplace=True)
        scaler = StandardScaler()
        data_scaled = scaler.fit_transform(data)

        reducer = umap.UMAP(
                            n_neighbors=15,    # Local neighborhood size (try 10-50)
                            min_dist=0.1,     # How tightly points are packed (0.0–0.5)
                            n_components=2,   # 2D or 3D output
                            random_state=42
                            )
        embedding = reducer.fit_transform(data_scaled)
        plt.figure(figsize=(8, 6))
        sns.scatterplot(x=embedding[:, 0], y=embedding[:, 1], hue=group, palette="Set1", s=100, alpha=0.7)
        plt.title('UMAP projection (unsupervised)')
        plt.xlabel('UMAP 1')
        plt.ylabel('UMAP 2')
        plt.show()

    def perform_PCA(self):
        #remove unwanted columns
        data = self.all_data
        col_to_drp = ['ID', 
                          'Replicate', 
                          "x center [micron]",
                          "y center [micron]",
                          "z center [micron]",
                          "Intensity threshold A",
                          "Intensity threshold B",
                          "Intensity threshold Basal Stain",
                          ]
        data = data.drop([col for col in col_to_drp if col in data.columns], axis=1)
        data.dropna(how='all', axis=1, inplace=True)
        data.dropna(axis=0, how="any", inplace=True)
        group = data['Group']
        data_temp = data.drop('Group', axis=1)

        scaler = StandardScaler()
        data_scaled = scaler.fit_transform(data_temp)
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(data_scaled)
        pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'])
        pca_df.reset_index(drop=True, inplace=True)
        group.reset_index(drop=True, inplace=True)
        pca_df.insert(0, "Group", group)
        plt.figure(figsize=(8, 6))
        sns.set_theme(style="whitegrid")
        sns.scatterplot(x=pca_df['PC1'], y=pca_df['PC2'], s=100, hue = pca_df['Group'], palette = "Set1", alpha = 0.7)
        plt.gca().set(title = 'Principle Component Analysis', xlabel = 'Principal Component 1', ylabel = 'Principal Component 2')
        plt.show()
