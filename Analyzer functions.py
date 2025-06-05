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


if __name__ == "__main__":
    pass

def import_ciliaQ_data(directories, column_names):
    all_groups = []
    all_data = pd.DataFrame(columns=column_names)
    for directory in enumerate(directories):
        group_name = slice_filename(directory[1])
        dataGroup = DataGroup(group_name, directory[1], column_names)
        dataGroup.import_group(directory[1], column_names, group_name)
        all_groups.append(dataGroup)
    #for group in all_groups:
        all_data = pd.concat([all_data, dataGroup.getConcatenatedData()])
    return all_groups, all_data

#def poolReplicates(filepath, column_list):
#	group = pd.DataFrame(columns = column_list, index = None)
#	return group

def appendToMeasurement(listOfMeasurements, measurement):
    listOfMeasurements.append(measurement)

def selectFolder(directoryList, label):
    filepath = filedialog.askdirectory()
    directoryList.append(filepath)
    label.configure(text = "\n".join(directoryList), background = "white")

def slice_filename(str):
    slash_to_Slice = str.rfind("/")
    sliced_string =  str[slash_to_Slice + 1 :]
    return sliced_string
    
def trim_Group_Replicate(data_list):
    if "Group" in data_list: 
        data_list.remove("Group")
    if "Replicate" in data_list: 
        data_list.remove("Replicate") 
    return data_list       

class DataGroup:
    def __init__(self, name_param, directory_param, classes):
        self.directory = directory_param
        self.classList = classes
        self.pooledData = pd.DataFrame(columns = self.classList)
        self.singleReplicates = []
        self.name = name_param
        self.ciliaQ_versions = []

    def get_ciliaQ_versions(self):
        return self.ciliaQ_versions

    def import_group(self, dir, classList, group_name):
        replicate_ID = 1
        for file in os.scandir(dir):
            if "CQ.txt" in file.name:
                rep = Replicate(file, classList, replicate_ID, group_name)
                rep.import_metadata()
                rep.import_Data()
                rep_data = rep.get_data() 
                self.singleReplicates.append(rep) 
                rep_data.reset_index(drop=True, inplace = True)
                aligned_data = rep_data.reindex(columns=self.pooledData.columns)
                self.pooledData = pd.concat([self.pooledData, aligned_data], ignore_index=True)
                self.ciliaQ_versions.append(rep.get_ciliaQ_version())
                replicate_ID = replicate_ID + 1

    def getConcatenatedData(self):
        return self.pooledData
    
    def get_Replicate_Data(self):
        return self.singleReplicates
    
    def getIndex(self):
        return self.index_in_anaysis

    def getName(self):
        return self.name

class Replicate(DataGroup):
    def __init__(self, file_param, classes, id, group):
        self.file = file_param   
        self.classes = classes
        self.rep_ID = id
        self.group_name = group        
        self.image_name = None
        self.skl_gauss_sigma_XY = None
        self.skl_gauss_sigma_Z = None
        self.excluded_cilia = None
        self.analyzed_cilia = None
        self.calibration_XY = None
        self.calibration_Z = None
        self.data_df = pd.DataFrame()
        self.ciliaQ_version = None

    def import_metadata(self):
        meta_dic = {}
        with open(self.file.path, "r", encoding = "utf-8") as f:
            lines = [line.strip().split("\t") for line in f]
            if (len(lines)) != 0:
                self.image_name = str(lines[1][1])
                #print("image name: " + self.image_name) 
                self.skl_gauss_sigma_XY = float(lines[20][1])
                #print("self.skl_gauss_sigma_XY: " + str(self.skl_gauss_sigma_XY))
                self.skl_gauss_sigma_Z = float(lines[21][1])
                #print("self.skl_gauss_sigma_Z: " + str(self.skl_gauss_sigma_Z))
                self.excluded_cilia = int(lines[18][3])
                #print("self.excluded_cilia: " + str(self.excluded_cilia))
                self.analyzed_cilia = int(lines[18][5])
                #print("self.analyzed_cilia: " + str(self.analyzed_cilia))
                self.calibration_XY = float(lines[8][1])
                #print("self.calibration_XY: " + str(self.calibration_XY))
                self.calibration_Z = float(lines[8][3])
                #print("self.calibration_Z: " + str(self.calibration_Z))

    def import_Data(self):
        try:
            dataset = pd.read_csv(
                    self.file.path, 
                    index_col= False, 
                    sep = "\t", 
                    skiprows = 28, 
                    usecols=range(60)).dropna(axis="columns", how = "all")
            self.ciliaQ_version = dataset.iloc[-1, 1]
            dataset = dataset.loc[:, ~dataset.columns.str.contains("^Unnamed", na=False)]
            dataset.drop(dataset.tail(3).index, inplace = True) 
            dataset["Group"] = str(self.group_name)
            dataset["Replicate"] = self.rep_ID
            self.data_df = dataset
        except FileNotFoundError:
            print("File not found. Please enter a folder that contains .CQ files")

    def get_skl_gauss_XY(self):
        return self.skl_gauss_sigma_XY
    
    def get_skl_gauss_Z(self):
        return self.skl_gauss_sigma_Z
    
    def get_image_name(self):
        return self.image_name
    
    def get_ciliaQ_version(self):
        return self.ciliaQ_version
    
    def get_excluded_cilia(self):
        return self.excluded_cilia
    
    def get_analyzed_cilia(self):
        return self.analyzed_cilia
    
    def get_cal_XY(self):
        return self.calibration_XY
    
    def get_cal_Z(self):
        return self.calibration_Z
    
    def get_data(self):
        return self.data_df
    
    def get_id(self):
        return self.rep_ID
    
class Plotter:
    def __init__(self, colors):
        self.colour_palette = colors
    
    def plotBoxplot(self, data_df, measurement_name, y_axis_names_dic):
        data_df = data_df.reset_index(drop=True)
        if data_df.empty:
            print(f"Measurement {measurement_name} not found!")
        else:
            data_df["value"] = pd.to_numeric(data_df["value"], errors='coerce')
            data_df = data_df.dropna(subset = ["value"])
            plt.figure(figsize=(6, 4))
            sns.set_theme(style="whitegrid")
            ax_box = sns.boxplot(data = data_df, hue = "variable", x = "variable", y = "value", palette="magma") 
            ax_box.set_title(measurement_name, fontweight = "bold")
            ax_box.set(xlabel = "", ylabel = y_axis_names_dic[measurement_name])
            ax_box.tick_params(labelrotation = 30)
            sns.stripplot(x = "variable", y = "value", hue = "variable", data=data_df, palette='dark:black', size=3, ax = ax_box)
            plt.xlabel(None)
            plt.show()
            #plt.savefig(f"{measurement_name}", format="svg")

    def plotViolinplot(self, data_df, measurement_name, y_axis_names_dic):
        if data_df.empty:
            print(f"Measurement {measurement_name} not found!")
        else:
            plt.figure(figsize=(6, 4))
            sns.set_theme(style="whitegrid")
            ax_vio = sns.violinplot(data = data_df, x = "variable", y = "value", hue= "variable", palette="magma")
            ax_vio.set_title(measurement_name)
            ax_vio.set(xlabel = "", ylabel = y_axis_names_dic[measurement_name])
            sns.stripplot(x="variable", y = "value", data=data_df, palette='dark:black', size=3, ax=ax_vio)
            plt.xlabel(None)
            plt.show()        
            #plt.savefig(f"{measurement_name}", format="svg")

    #def superplot(self, combined, measurement, y_axis_names_dic):
    #    sns.set_theme(style="whitegrid")
    #    combined = combined.reset_index(drop = True)
    #    combined[measurement] = pd.to_numeric(combined[measurement], errors='coerce')
    #    combined = combined.dropna(subset=[measurement])
    #    if(not combined[measurement].empty):
    #        ReplicateAverages = combined.groupby(['Group', 'Replicate'], as_index=False).agg({str(measurement) : "mean"}).reset_index()
    #        sns.stripplot(x="Group", y=measurement, hue="Replicate", data=combined)
    #        ax = sns.swarmplot(x = "Group", y=measurement, hue="Replicate", size=15, edgecolor="k", linewidth=2, data=ReplicateAverages)
    #        ax.legend_.remove()
    #        ax.set_title(measurement)
    #        ax.set(xlabel = "", ylabel = y_axis_names_dic[measurement])
    #        ax.tick_params(labelrotation = 45)
    #        plt.show()
    #    else:
    #        print(measurement + " was not found")

    def superplot4(self, combined, measurement, y_axis_names_dic, stat_list):
        #stat_measurement = [dic for dic in stat_list if dic.get("measurement") == measurement]
        stat_measurement = next((dic for dic in stat_list if dic.get("measurement") == measurement), None)
        print(stat_measurement)
        sns.set_theme(style="whitegrid")
        combined = combined.reset_index(drop = True)
        combined[measurement] = pd.to_numeric(combined[measurement], errors='coerce')
        combined = combined.dropna(subset=[measurement])
        if(not combined[measurement].empty):
            ReplicateAverages = combined.groupby(['Group', 'Replicate'], as_index=False).agg({str(measurement) : "mean"}).reset_index()
            sns.stripplot(x="Group", y=measurement, hue="Replicate", data=combined)
            ax = sns.swarmplot(x = "Group", y=measurement, hue="Replicate", size=15, edgecolor="k", linewidth=2, data=ReplicateAverages)
            ax.legend_.remove()
            ax.set_title(measurement)
            ax.set(xlabel = "", ylabel = y_axis_names_dic[measurement])
            ax.tick_params(labelrotation = 45)

            if stat_measurement is not None:
                if stat_measurement["pairwise"] is not None:

                    x_ticks_group = ax.get_xticks()
                    xticklabels = [tick.get_text() for tick in ax.get_xticklabels()]
                    group_positions = {grp: pos for grp, pos in zip(xticklabels, x_ticks_group)}
                    max_y_val = combined.groupby("Group")[measurement].max() 
                    overall_range = combined[measurement].max() - combined[measurement].min()
                    offset = overall_range * 0.10
                    for comparison, values in stat_measurement["pairwise"].items():
                        try:
                            group1, group2 = comparison.split(" vs ")
                        except ValueError:
                            print(f"Comparison '{comparison}' not formatted as 'Group1 vs Group2'. Skipping this comparison.")
                            continue  
                    # Get x positions for the groups
                    x1 = group_positions.get(group1)
                    x2 = group_positions.get(group2)
                    if x1 is None or x2 is None:
                        print(f"Could not find x positions for comparison '{comparison}'.")  
                    y1 = max_y_val.get(group1, 0)
                    y2 = max_y_val.get(group2, 0)
                    y = max(y1, y2) + offset  ###Remove if necessary

                    #add padding above statistics bar
                    current_ax_lim = ax.get_ylim()  # returns (lower, upper)
                    lower, upper = current_ax_lim
                    new_ax_lim = upper + offset * 2  # offset is already a float
                    ax.set_ylim(lower, new_ax_lim)

                    ax.plot([x1, x2], [y, y], lw=1.5, c='black')
                    ax.plot([x1, x1], [y, y - offset/2], lw=1.5, c='black')
                    ax.plot([x2, x2], [y, y - offset/2], lw=1.5, c='black')
                    ax.text((x1 + x2) / 2, y + offset/10, values[0], ha='center', va='bottom', color='black', fontsize=15, weight='bold')
            plt.show()
        else:
            print(measurement + " was not found")



    def superplot5(self, combined, measurement, y_axis_names_dic, stat_list, directory_for_saving):
        stat_measurement = next((dic for dic in stat_list if dic.get("measurement") == measurement), None)
        print(stat_measurement)
        sns.set_theme(style="whitegrid")
        combined = combined.reset_index(drop = True)
        combined[measurement] = pd.to_numeric(combined[measurement], errors='coerce')
        combined = combined.dropna(subset=[measurement])

        if(not combined[measurement].empty):
            ReplicateAverages = combined.groupby(['Group', 'Replicate'], as_index=False).agg({str(measurement) : "mean"}).reset_index()

            filename = re.sub(r'[\\/*?:"<>|]', "_", measurement)
            write_path = directory_for_saving + f"/{filename}.csv"
            with open(write_path, "w") as csv: 
                csv.write(f"Plotting results for {measurement}\n")
                csv.write(f"{os.path.abspath(write_path)} \t")
                csv.write("\n")
                combined[[measurement, 'Group', 'Replicate']].to_csv(write_path, mode = "a", sep=",", index = False)
                      
            fig, ax = plt.subplots()
            sns.violinplot(x="Group", y=measurement, data=combined, ax=ax, color="lightgrey", inner="quart", alpha = 0.6, density_norm='count')
            sns.stripplot(x="Group", y=measurement, hue="Replicate", data=combined)
            sns.swarmplot(x = "Group", y=measurement, hue="Replicate", size=15, edgecolor="k", linewidth=2, data=ReplicateAverages, ax = ax)
            ax.legend_.remove()
            ax.set_title(measurement, weight = "bold")
            ax.set(xlabel = "", ylabel = y_axis_names_dic[measurement])
            ax.tick_params(labelrotation = 45)
            
            if stat_measurement is not None:
                if stat_measurement["pairwise"] is not None:
                    label_performed_test = f"Statistical test: {stat_measurement['hypothesis test']}"
                    #fig.text(0.01, 0.01, label_performed_test, ha='left', va='bottom', fontsize=10, ) # NEW
                    fig.suptitle(label_performed_test, x=0.01, y=-0.03, ha='left', va='bottom', fontsize=10)
                    x_ticks_group = ax.get_xticks()
                    xticklabels = [tick.get_text() for tick in ax.get_xticklabels()]
                    group_positions = {grp: pos for grp, pos in zip(xticklabels, x_ticks_group)}
                    max_y_val = combined.groupby("Group")[measurement].max() 
                    overall_range = combined[measurement].max() - combined[measurement].min()
                    offset = overall_range * 0.10
                    for comparison, values in stat_measurement["pairwise"].items():
                        try:
                            group1, group2 = comparison.split(" vs ")
                        except ValueError:
                            print(f"Comparison '{comparison}' not formatted as 'Group1 vs Group2'. Skipping this comparison.")
                            continue  
                    # Get x positions for the groups
                    x1 = group_positions.get(group1)
                    x2 = group_positions.get(group2)
                    if x1 is None or x2 is None:
                        print(f"Could not find x positions for comparison '{comparison}'.")  
                    y1 = max_y_val.get(group1, 0)
                    y2 = max_y_val.get(group2, 0)
                    y = max(y1, y2) + offset * 1.5

                    #add padding above statistics bar
                    current_ax_lim = ax.get_ylim()  # returns (lower, upper)
                    lower, upper = current_ax_lim
                    new_ax_lim = upper + offset * 2  # offset is already a float
                    ax.set_ylim(lower, new_ax_lim)
                    ax.plot([x1, x2], [y, y], lw=1.5, c='black')
                    ax.plot([x1, x1], [y, y - offset/2], lw=1.5, c='black')
                    ax.plot([x2, x2], [y, y - offset/2], lw=1.5, c='black')
                    ax.text((x1 + x2) / 2, y + offset/10, values[0], ha='center', va='bottom', color='black', fontsize=15)
            plt.show()
        else:
            print(measurement + " was not found")









    def superplot3(self, combined, measurement, y_axis_names_dic, statistics_df):
        sns.set_theme(style="whitegrid")
        combined = combined.reset_index(drop = True)
        combined[measurement] = pd.to_numeric(combined[measurement], errors='coerce')
        combined = combined.dropna(subset=[measurement])
        if(not combined[measurement].empty):
            ReplicateAverages = combined.groupby(['Group', 'Replicate'], as_index=False).agg({str(measurement) : "mean"}).reset_index()
            sns.stripplot(x="Group", y=measurement, hue="Replicate", data=combined)
            ax = sns.swarmplot(x = "Group", y=measurement, hue="Replicate", size=15, edgecolor="k", linewidth=2, data=ReplicateAverages)
            ax.legend_.remove()
            ax.set_title(measurement)
            ax.set(xlabel = "", ylabel = y_axis_names_dic[measurement])
            ax.tick_params(labelrotation = 45)
            plt.show()
        else:
            print(measurement + " was not found")

    #def superplot2(self, combined, measurement, y_axis_names_dic):
    #    sns.set_theme(style="whitegrid")
    #    combined = combined.reset_index(drop = True)
    #    combined[measurement] = pd.to_numeric(combined[measurement], errors='coerce')
    #    combined = combined.dropna(subset=[measurement])
    #    if(not combined[measurement].empty):
    #        ReplicateAverages = combined.groupby(['Group', 'Replicate'], as_index=False).agg({str(measurement) : "mean"}).reset_index()
    #        sns.boxplot(
    #                x="Group",
    #                y=measurement,
    #                data=combined,
    #                showfliers=True,           # Keep outliers visible
    #                whis=1.5,                  # 1.5 * IQR rule
    #                boxprops={"facecolor": "none", "edgecolor": "none"},
    #                whiskerprops={"visible": False},
    #                capprops={"visible": False},
    #                medianprops={"visible": False},
    #                flierprops={"marker": "x", "markerfacecolor": "none", "markersize": 7},
    #                )
    #        sns.stripplot(x="Group", 
    #                      y=measurement, 
    #                      hue="Replicate", 
    #                      data=combined)
    #        ax = sns.swarmplot(x = "Group", y=measurement, hue="Replicate", size=15, edgecolor="k", linewidth=2, data=ReplicateAverages)
    #        ax.legend_.remove()
    #        ax.set_title(measurement)
    #        ax.set(xlabel = "", ylabel = y_axis_names_dic[measurement])
    #        ax.tick_params(labelrotation = 45)
    #        plt.show()
    #    else:
    #        print(measurement + " was not found")

class Measurements: 
    def __init__(self):
        self.measures_dic = {
            "ID" : ["count", "morph"],
            "x center [micron]": ["x-center coordinate [µm]", "spatial"],
            "y center [micron]": ["y-center coordinate [µm]", "spatial"],
            "z center [micron]": ["z-center coordinate [µm]", "spatial"],
            "Volume [voxel]": ["Cilium Volume [voxel]", "morph"],
            "Volume [micron^3]": ["Cilium Volume [µm^3]", "morph"],
            "# surface voxels": ["Cilium surface [voxels]", "morph"],
            "Surface [micron^2]": ["Cilium surface voxels [µm^2]", "morph"],
            "Shape complexity index": ["Shape complexity [a.u.]", "morph"],
            "Sphere radius [micron]": ["Sphere radius [µm]", "morph"],
            "Maximum span [micron]": ["Maximum span [µm]", "morph"],
            "A: Colocalized volume [micron^3] (if channel in input image was background-removed)": ["Colocalized volume [µm^3]", "channel_A_B"],
            "A: Colocalized volume [% total volume] (if channel in input image was background-removed)": ["Colocalized volume [%]", "channel_A_B"],
            "B: Colocalized volume [micron^3] (if channel in input image was background-removed)": ["Colocalized volume [µm^3]", "channel_A_B"],
            "B: Colocalized volume [% total volume] (if channel in input image was background-removed)": ["Colocalized volume [%]", "channel_A_B"],
            "A: Colocalized compared to BG volume [micron^3]": ["Colocalized compared to BG volume [µm^3]", "channel_A_B"],
            "A: Colocalized compared to BG volume [% total volume]": ["Colocalized compared to BG volume [%]", "channel_A_B"],
            "B: Colocalized compared to BG volume [micron^3]": ["Colocalized compared to BG volume [µm^3]", "channel_A_B"],
            "B: Colocalized compared to BG volume [% total volume]": ["Colocalized compared to BG volume [%]", "channel_A_B"],
            "minimum intensity (in reconstruction channel)": ["Minimum intensity [a.u.]", "cilia_marker"],
            "maximum intensity (in reconstruction channel)": ["Maximum intensity [a.u.]", "cilia_marker"],
            "average intensity of the 10% of voxels with highest intensity (in reconstruction channel)": ["average intensity [a.u.]", "cilia_marker"],
            "average intensity (in reconstruction channel)": ["average intensity [a.u.]", "cilia_marker"],
            "SD of intensity (in reconstruction channel)": ["SD of intensity [a.u.]", "cilia_marker"],
            "minimum A intensity": ["Minimum intensity [a.u.]", "channel_A_B"],
            "maximum A intensity": ["Maximum intensity [a.u.]", "channel_A_B"],
            "average A intensity of the 10% of voxels with highest A intensity": ["average intensity [a.u.]", "channel_A_B"],
            "average A intensity": ["Average intensity [a.u.]", "channel_A_B"],
            "SD of A intensity": ["SD of intensity [a.u.]", "channel_A_B"],
            "minimum B intensity": ["Minimum intensity [a.u.]", "channel_A_B"],
            "maximum B intensity": ["Maximum intensity[a.u.]", "channel_A_B"],
            "average B intensity of the 10% of voxels with highest B intensity": ["[a.u.]", "channel_A_B"],
            "average B intensity": ["[a.u.]", "channel_A_B"],
            "SD of B intensity": ["[a.u.]", "channel_A_B"],
            "# of found skeletons (quality parameter)" : ["[a.u.]", "morph"],
            "# branches (quality parameter)" : ["[a.u.]", "morph"],
            "tree length [micron] (quality parameter)" : ["[a.u.]", "morph"],
            "cilia length [micron] (largest shortest path of largest skeleton)" : ["[a.u.]", "morph"],
            "orientation vector x [micron] (vector from first to last skeleton point)" : ["[a.u.]", "spatial"],
            "orientation vector y [micron] (vector from first to last skeleton point)" : ["[a.u.]", "spatial"],
            "orientation vector z [micron] (vector from first to last skeleton point)" : ["[a.u.]", "spatial"],
            "cilia bending index (arc length of cilium / eucledian distance of first and last skeleton point)" : ["[a.u.]", "morph"],
            "Intensity threshold A" : ["[a.u.]", "channel_A_B"],
            "Intensity threshold B" : ["[a.u.]", "channel_A_B"],
            "Intensity threshold Basal Stain" : ["[a.u.]", "channel_A_B"],
            "Integrated A intensity" : ["[a.u.]", "channel_A_B"],
            "Average A intensity on center line" : ["[a.u.]", "channel_A_B"],
            "Integrated B intensity" : ["[a.u.]", "channel_A_B"],
            "Average A intensity on center line" : ["[a.u.]", "channel_A_B"],
            "A: Colocalized on centerline compared to BG volume [micron]" : ["[a.u.]", "channel_A_B"],
            "A: Colocalized on centerline compared to BG volume [% total length]" : ["[a.u.]", "channel_A_B"],
            "B: Colocalized on centerline compared to BG volume [micron]" : ["[a.u.]", "channel_A_B"],
            "B: Colocalized on centerline compared to BG volume [% total length]" : ["[a.u.]", "channel_A_B"],
        }   

    def get_measurement_dic(self):
        dic_temp = {key: values[0] for key, values in self.measures_dic.items()}
        return dic_temp
    
    def get_measurement_category(self):
        dic_temp = {key: values[1] for key, values in self.measures_dic.items()}
        return dic_temp 
   
    def get_measurement_list(self):
        measurement_list = list(self.measures_dic.keys())
        measurement_list += ["Group", "Replicate"]
        return list(measurement_list)
    
class Statistics():
    def __init__(self, groups_param, data_param, selected_measurements_param):
        self.all_groups = groups_param
        self.all_data = data_param
        self.tukey_df = pd.DataFrame()
        self.dunn_df = pd.DataFrame()
        self.column_names = selected_measurements_param
        print(self.column_names)
        self.measurement_selection = self.column_names
        #self.measurement_selection = self.column_names.remove(["Replicate", "Group"])
        # prepare the dataframe
        for measurement in self.measurement_selection:
            if (measurement not in ["Replicate", "Group"]):
                if(not self.all_data[measurement].empty):
                    self.all_data[measurement] = pd.to_numeric(self.all_data[measurement], errors="coerce") 

    def get_dunn_results(self):
        return self.dunn_df
    
    def get_tukey_results(self):
        return self.tukey_df

    
    def hypothesis_testing(self, data, alpha, measurement_name):
        #data preprocessing
        conditions = data.columns.tolist()
        measurement = measurement_name
        data = data.apply(pd.to_numeric, errors='coerce')
        n_conditions = len(conditions)


        #2 .test for normal distributed data
        result_norm_dist = {}
        for condition in conditions:
            data_norm_dist = data[condition].dropna()
            #print("now comes data_norm_dist")
            #print(data_norm_dist)
            #print("Datanormdist end")
            if (len(data_norm_dist) >= 3):
                stat, p_shapiro = stats.shapiro(data_norm_dist)
                result_norm_dist
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
            #3. test for homogeniety of variance
            if (len(cond1) >= 3 and len(cond2) >= 3):
                stat, p_levene = stats.levene(cond1, cond2)
                result_homogen.update({f"condition": p_levene})
            else:
                print(f"The group {condition} contains not enough cilia to run hypothesis testing. {condition} was excluded from hypothesis testing ---")

            homogen_dist_summary = all(value > alpha for key, value in result_homogen.items())
            summary.update({"Homogeniety of Variance": homogen_dist_summary})

            #perform hypothesis testing    
            if((summary["Norm. distribution"] == True) & (summary["Homogeniety of Variance"] == True)):
                stat, p = stats.ttest_ind(cond1, cond2)
                performed_test = "t-Test"
            elif((summary["Norm. distribution"] == False) & (summary["Homogeniety of Variance"] == False) | 
                 (summary["Norm. distribution"] == False) & (summary["Homogeniety of Variance"] == True)):
                stat, p = stats.mannwhitneyu(cond1, cond2)
                performed_test = "Mann-Witney-U Test"
            elif((summary["Norm. distribution"] == True) & (summary["Homogeniety of Variance"] == False)):
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

            #3. test for homogeniety of variance
            if (len(cond1) >= 3 and len(cond2) >= 3):
                stat, p_levene = stats.levene(cond1, cond2)
                result_homogen.update({f"condition": p_levene})
            else:
                print(f"The group {condition} contains not enough cilia to run hypothesis testing. {condition} was excluded from hypothesis testing ---")

            homogen_dist_summary = all(value > alpha for key, value in result_homogen.items())
            summary.update({"Homogeniety of Variance": homogen_dist_summary})

            if((summary["Norm. distribution"] == True) & (summary["Homogeniety of Variance"] == True)):
                stat, p = stats.f_oneway(*data_col_list)
                performed_test = "One-way ANOVA"
            elif((summary["Norm. distribution"] == False) & (summary["Homogeniety of Variance"] == True)):
                stat, p = stats.kruskal(*data_col_list)
                performed_test = "Kruskal-Test"
            elif((summary["Norm. distribution"] == False) & (summary["Homogeniety of Variance"] == False)):
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
            if((performed_test == "One-way ANOVA") & (p < 0.05)):
                melted_data = data.melt(var_name = "group", value_name = "value")
                tukey_results = pairwise_tukeyhsd(melted_data["value"], melted_data["group"], alpha=alpha)
                tukey_summary = tukey_results.summary().data[1:]
                #print(tukey_summary)
                for i in range(len(tukey_results.pvalues)): #this is nothing else than number of comparisons
                    for comparison in tukey_summary:
                        cond1 = comparison[0]
                        cond2 = comparison[1]            
                        p_value = float(comparison[3])           
                        asterisk = self.get_asterisk(p_value)
                        pair = f"{cond1} vs {cond2}"
                        pairwise_results[pair] = [asterisk, p_value]
                        #print(f" This is the reuslt for: {cond1} vs {cond2}")
                        #print(f"performed test: {performed_test} and Tukey's")
                        #print(f"p-value: {p_value}")
                        #print(f"asterik: {asterisk}")
                 
                    #pair = f"{tukey_results.groupsunique[tukey_results.pairindices[i, 0]]} vs {tukey_results.groupsunique[tukey_results.pairindices[i, 1]]}"
                    #pairwise_results.update({pair: [self.get_asterisk(tukey_results.pvalues[i]), tukey_results.pvalues[i]]})

            elif (performed_test == "Kruskal-Wallis" and p < alpha):
                comparisons = list(itertools.combinations(conditions, 2)) 
                p_values = []
                for cond1, cond2 in comparisons:
                    stats_temp, p = stats.mannwhitneyu(data[cond1].dropna(), data[cond2].dropna(), alternative = "two-sided")
                    p_values.append(p)
            
                #bonferroni correction
                corrected_p_values = multipletests(p_values, alpha=alpha, method='bonferroni')[1]
                for (cond1, cond2), p_corrected in zip(comparisons, corrected_p_values):
                    pair = f"{cond1} vs {cond2}"
                    pairwise_results.update({pair: [self.get_asterisk(p_corrected), p_corrected]})
            summary.update({"pairwise": pairwise_results})
   
            return summary

        
    def get_asterisk(self, p_value):
        if(p_value < 0.001):
            return "***"
        elif(p_value < 0.01):
            return "**"
        elif(p_value < 0.05):
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
        des_stat = pd.DataFrame(columns = ["Measurement", "Group", "Mean", "Median", "Min", "Max", "Count"])
        for measurement in self.measurement_selection:
            if(measurement not in ["Replicate", "Group"]):
                stats_temp = self.all_data.groupby("Group")[measurement].agg(["mean", "median", "std", "min", "max", "count"])
                stats_temp.reset_index(inplace=True)
                stats_temp.columns = stats_temp.columns.str.capitalize()
                stats_temp.insert(0, "Measurement", measurement)
                #for group, values in stats_temp.to_dict(orient="index").items():
                #to_concat = pd.DataFrame(columns= [""Measurement", "Group", "Mean", "Median", "Min", "Max", "Count""])
                #to_add_temp = {"Measurement": measurement, "Group": group,   **values}
                des_stat = pd.concat([des_stat, stats_temp], ignore_index = True)
        return des_stat
    
    def perform_UMAP(self):
        data = self.all_data
        data = data.drop(['ID', 
                          'Replicate', 
                          "x center [micron]",
                          "y center [micron]",
                          "z center [micron]",
                          "Intensity threshold A",
                          "Intensity threshold B",
                          "Intensity threshold Basal Stain",
                          ], axis=1)
        data.dropna(how='all', axis=1, inplace=True)
        data.dropna(axis=0, how="any", inplace=True)
        group = data['Group']
        data = data.drop('Group', axis=1)

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
        data = data.drop(['ID', 
                          'Replicate', 
                          "x center [micron]",
                          "y center [micron]",
                          "z center [micron]",
                          "Intensity threshold A",
                          "Intensity threshold B",
                          "Intensity threshold Basal Stain",
                          ], axis=1)
        data.dropna(how='all', axis=1, inplace=True)
        data.dropna(axis=0, how="any", inplace=True)
        group = data['Group']
        data = data.drop('Group', axis=1)

        #print(data.isna().sum())
        scaler = StandardScaler()
        data_scaled = scaler.fit_transform(data)
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(data_scaled)
        pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'])
        pca_df.reset_index(drop=True, inplace=True)
        group.reset_index(drop=True, inplace=True)
        pca_df.insert(0, "Group", group)
        plt.figure(figsize=(8, 6))
        sns.set_theme(style="whitegrid")
        sns.scatterplot(x=pca_df['PC1'], y=pca_df['PC2'], s=100, hue = pca_df['Group'], palette = "Set1", alpha = 0.7)
        plt.gca().set(title = 'PCA', xlabel = 'Principal Component 1', ylabel = 'Principal Component 2')
        plt.show()


class QC:   
    thresholds = {
        "max cilia length": 15,
        "max branches": 3,
        "max volume": 25,
        "max span": 1
    }
    def __init__(self, all_groups):
        self.groups = all_groups

    def perform_qc(self):
        errors = []
        ciliaQ_versions_all_groups = []
        for group in self.groups:
            ciliaQ_versions_all_groups.append(group.get_ciliaQ_versions())
            df = group.getConcatenatedData()
            #check for cilia length outliers
            cilia_length_outliers = self.detect_outlier(df, "cilia length [micron] (largest shortest path of largest skeleton)")
            if(len(cilia_length_outliers) != 0):  
                errors.append("The cilia in group " + group.getName() + " contains " + str(len(cilia_length_outliers[0])) + " outliers in cilia length. The affected cilia are " + str(cilia_length_outliers[0]) + ". Please verify the segmentation to prevent accidentally fused cilia")
            
            #check for high-volume cilia
            cilia_volume_outliers = self.detect_outlier(df, "Volume [micron^3]")
            if(len(cilia_volume_outliers) != 0):  
                errors.append("The cilia in group " + group.getName() + " contains " + str(len(cilia_volume_outliers[0])) + " outliers in cilia length. The affected cilia are " + str(cilia_volume_outliers[0]) + ". Please verify the segmentation to prevent accidentally fused cilia")
            
            
            #high_vol_cilia = self.compare_to_threshold(df, "Volume [micron^3]", self.thresholds["max volume"], ">")
            #if(len(high_vol_cilia) != 0):
            #    errors.append("The group " + group.getName() + " contains " + len(high_vol_cilia) + " cilia with exceptionally high volume. The affected cilia are " + str(high_vol_cilia) + ". Please verify the segmentation to prevent accidentally fused cilia.")
            
            #high_length_cilia = self.compare_to_threshold(df, "cilia length [micron] (largest shortest path of largest skeleton)", self.thresholds["max cilia length"], ">")
            #if(len(high_length_cilia) != 0):
            #    errors.append("The group " + group.getName() + " contains " + len(high_length_cilia) + " cilia with exceptionally high length. The affected cilia are " + str(high_length_cilia) + ". Please verify the segmentation to prevent accidentally fused cilia.")
            
            #I would include this
            high_branched_cilia = self.compare_to_threshold(df, "# branches (quality parameter)", self.thresholds["max branches"], ">")
            if(len(high_branched_cilia) != 0):
                errors.append("The group " + group.getName() + " contains " + str(len(high_branched_cilia)) + " exceptionally branched cilia. The affected cilia are " + str(high_branched_cilia) + ". Please verify the segmentation to prevent accidentally fused cilia.")

            cilia_branch_outliers = self.detect_outlier(df, "# branches (quality parameter)")
            if(len(cilia_branch_outliers) != 0):  
                errors.append("The cilia in group " + group.getName() + " contains " + str(len(cilia_branch_outliers[0])) + " outliers in cilia length. The affected cilia are " + str(cilia_branch_outliers[0]) + ". Please verify the segmentation to prevent accidentally fused cilia")

            cilia_maxspan_outliers = self.detect_outlier(df, "Maximum span [micron]")
            if(len(cilia_maxspan_outliers) != 0):  
                errors.append("The cilia in group " + group.getName() + " contains " + str(len(cilia_maxspan_outliers[0])) + " outliers in cilia length. The affected cilia are " + str(cilia_maxspan_outliers[0]) + ". Please verify the segmentation to prevent accidentally fused cilia")

            #high_maxspan_cilia = self.compare_to_threshold(df, "Maximum span [micron]", self.thresholds["max span"], ">") # also take outliers
            #if(len(high_maxspan_cilia) != 0):
            #    errors.append("The group " + group.getName() + " contains " + len(high_maxspan_cilia) + " exceptionally branched cilia. The affected cilia are " + str(high_maxspan_cilia) + ". Please verify the segmentation to prevent accidentally fused cilia.")

            averageInt_outliers_cilia = self.detect_outlier(df, "average A intensity of the 10% of voxels with highest A intensity")
            if(len(averageInt_outliers_cilia[0]) != 0):
                errors.append("The group " + group.getName() + " contains " + str(len(averageInt_outliers_cilia[0])) + " exceptionally branched cilia. The affected cilia are " + str(averageInt_outliers_cilia[0]) + ". Please verify the segmentation to prevent accidentally fused cilia.")

            underexposed_cilia_px = self.compare_to_threshold(df, "minimum intensity (in reconstruction channel)", 0, "=")
            if(len(underexposed_cilia_px) != 0):
                errors.append("The group " + group.getName() + " contains " + str(len(underexposed_cilia_px)) + " cilia that have underexposed pixels as part of their mask. The affected cilia are cilia " + str(underexposed_cilia_px) + ". This might indicate that the ciliary masks are too large. Please check the segmentation.")

            avInt_variance_outliers_cilia = self.detect_outlier(df, "average intensity (in reconstruction channel)")
            if(len(avInt_variance_outliers_cilia[0]) != 0):
                errors.append("The group " + group.getName() + " contains " + str(len(avInt_variance_outliers_cilia[0])) + " exceptionally branched cilia. The affected cilia are " + str(avInt_variance_outliers_cilia[0]) + ". Please verify the segmentation to prevent accidentally fused cilia.")

        #if (len(set(ciliaQ_versions_all_groups)) > 1): 
        #	errors.append("You are about to pool data from different CiliaQ-versions. Since this is not advisable, please recheck your ciliaQ analysis.")
        self.raiseError(errors)

    def raiseError(self, text):
        root = Tk()
        root.lift()
        root.geometry
        root.title = ("Raised Errors:")
        frameForMessage = ttk.Frame(root, padding = 10)
        frameForMessage.pack(fill="x", expand=True)
        for error in text:
            errorLabel = Label(frameForMessage, text = error, bg="white", anchor="w").pack()
        ButtonSelectFolder = Button(root, text = "Ignore", command= root.destroy, font=("Helvetica", 14, "bold"), bd=2, relief="raised").pack()
        root.attributes("-topmost", True)
        root.mainloop()

    #returns a list of cilia-IDs of all outliers, the first quantile q1 and third q3  
    def detect_outlier(self, data, column):
        q1, q3 = np.percentile(data[column], [25, 75])
        iqr = q3 - q1
        #filtered_data = data.query(f"`{column}` < @q1 - 1.5 * @iqr | `{column}` > @q3 + 1.5 * @iqr")
        filtered_data = data[(data[column] < q1 - 1.5 * iqr) | (data[column] > q3 + 1.5 * iqr)]
        return [filtered_data["ID"].tolist(), q1, q3]
    
    # returns a list of cilia-IDs that fulfill the filter criteria. 
    # Example call: compare_to_threshold(df, 'length', 10, ">"). Returns all cilia IDs with more than 10uM length
    def compare_to_threshold(self, data, column, threshold, larger_smaller):
        #filtered_data = data.query(f"'{column}' {larger_smaller} {threshold}")
        filtered_data = data[(data[column] < threshold) if larger_smaller == "<" else (data[column] > threshold)]
        return filtered_data["ID"].tolist()

class setup():
    def __init__(self, all_measurements_param, measurement_category_param):
        self.all_measurements = all_measurements_param
        self.measurement_category = measurement_category_param
        self.measurement_selection = []
        self.stat_selection = "automated statistics (see publication)"
        self.perform_stat = False
        self.metafileDir = None

    def getDirectories(self):
        
        def selectFolder(directoryList, label):
            filepath = filedialog.askdirectory()
            directoryList.append(filepath)
            label.configure(text = "\n".join(directoryList), background = "white")

        directories = []
        root = Tk()
        root.geometry("600x600")
        root.lift()
        root.title = ("Please select folders")
        root.columnconfigure((0,1), weight = 1)
        root.rowconfigure((0,2), weight = 1)
        root.rowconfigure(1, weight = 10)
        frameForDirectories = Frame(root, bg= "white", relief="sunken", bd=2)
        frameForDirectories.grid(row = 1, column = 0, sticky = "news", columnspan=2)
        displayDirectoriesText = Label(root, text = "Selected directories")
        displayDirectoriesText.grid(row = 0, column = 0,columnspan=2, sticky = "nw")
        displayDirectories = Label(frameForDirectories)
        displayDirectories.grid(pady = 10, sticky = "n")
        ButtonSelectFolder = Button(root, text = "Select Folder", command= lambda: selectFolder(directories, displayDirectories), font=("Helvetica", 14, "bold"), bd=2, relief="raised").grid(row = 2, column = 0, sticky = "se")
        ButtonSelectFolder = Button(root, text = "Apply", command= root.destroy, font=("Helvetica", 14, "bold"), bd=2, relief="raised").grid(row = 2, column = 1, sticky = "se")
        root.attributes("-topmost", True)
        root.mainloop()
        return directories

    def setup_window(self):
        root = tk.Tk()
        root.resizable(True, True)
        root.title("Analyzer setup:")

        def select_stat(event): 
            self.stat_selection = drop_down.get()

        def selectMetafile():
            filepath = filedialog.askdirectory()
            self.metafileDir = filepath

        tk.Label(root, text =  "Include measurements into analysis:", font=("Helvetica", 12, "bold")).pack()
        Button(root, text =  "Select measurements", font=("Helvetica", 12, "bold"), bd=2, relief="raised", command = lambda: self.set_measurement_window()).pack()

        tk.Label(root, text =  "Upload FOR5547-Metafile", font=("Helvetica", 12, "bold")).pack(pady= 20)
        Button(root, text = "Select Folder", command= lambda : selectMetafile(), font=("Helvetica", 14, "bold"), bd=2, relief="raised").pack()

        tk.Label(root, text =  "Method for hypothesis testing", font=("Helvetica", 12, "bold")).pack(pady = 10)
        options = ["No statistics", "automated statistics (see publication)"]
        drop_down = ttk.Combobox(root, values = options, state = "readonly")
        drop_down.pack()
        drop_down.current(1)
        drop_down.bind("<<ComboboxSelected>>", select_stat)
        Button(root, text = "Apply", command= root.destroy, font=("Helvetica", 14, "bold"), bd=2, relief="raised").pack( padx = 10)

        root.lift()
        root.attributes("-topmost", True)
        root.mainloop()
        return self.measurement_selection, self.stat_selection, self.perform_stat, self.metafileDir
    
    def set_measurement_window(self):
        #self.all_measurements = self.all_measurements[0:5] #only for testing purposes
        window = tk.Toplevel()
        window.resizable(True, True)
        window.title("Please select Measurements")

        checkbox_state = {}

        def updateSelection(measurement, state):
            if state.get() == 1 and measurement not in self.measurement_selection:
                self.measurement_selection.append(measurement)
            if state.get() == 0 and measurement in self.measurement_selection:
                self.measurement_selection.remove(measurement)
            
        def choose_all():
            if all(checkbox_val.get() == 1 for checkbox_val in checkbox_state.values()):
                for measurement, result in checkbox_state.items():
                    result.set(0)
                self.measurement_selection.clear()
            else:
                for measurement, result in checkbox_state.items():
                    result.set(1)
                self.measurement_selection = self.all_measurements[:]

        # sort dictionary:
        categories_sorted = sorted(set(self.measurement_category.values()))
        print(categories_sorted)

        category_columns = {"channel_A_B": 1, "morph": 0, "spatial":0, "cilia_marker": 0}
        row_count = {0: 0, 1: 0} 

        for category in categories_sorted:
            #headline: 
            column = category_columns.get(category, 0)
            tk.Label(window, text = f"{category}:", font=("Helvetica", 12, "bold")).grid(row=row_count[column], column = column, sticky="w")

            row_count[column] = row_count[column]+1
            for measurement, cat in self.measurement_category.items():
                if(cat == category):
                    result = tk.IntVar()
                    checkbox_state[measurement] = result
                    Checkbutton(window, text= measurement,  variable = result, onvalue= 1, offvalue=0, command = lambda m=measurement, r = result: updateSelection(m, r)).grid(row=row_count[column], column = column, sticky = "w", padx = 5, pady = 0)
                    row_count[column] = row_count[column]+1
        Button(window, text = "Apply", command= window.destroy, font=("Helvetica", 14, "bold"), bd=2, relief="raised").grid(column = 1, row=row_count[0]+2, columnspan = 1, padx = 10, sticky = "s")
        Button(window, text =  "Select all", font=("Helvetica", 14, "bold"), bd=2, relief="raised", command = lambda: choose_all()).grid(row=row_count[1], column = 0, columnspan = 1, padx = 10, sticky = "s")
            
        window.lift()
        window.attributes("-topmost", True)
        window.mainloop()

    
    
    
    












