from tkinter import *
from tkinter import filedialog, ttk
import numpy as np
import csv

class QC:   
    thresholds = {
        "max cilia length": 15,
        "max branches": 3,
        "max volume": 25,
        "max span": 1
    }
    
    def __init__(self, all_groups, write_path_param, time_param):
        self.groups = all_groups
        self.write_path = write_path_param
        self.time = time_param
        self.errors = []

    def perform_qc(self):
        
        ciliaQ_versions_all_groups = []
        for group in self.groups:
            ciliaQ_versions_all_groups.append(group.get_ciliaQ_versions())
            df = group.getConcatenatedData()

            cilia_length_outliers = self.detect_outlier(df, "cilia length [micron] (largest shortest path of largest skeleton)")
            if(len(cilia_length_outliers[0]) != 0):  
                self.errors.append("The cilia in group " + group.getName() + " contains " + str(len(cilia_length_outliers[0])) + " outliers in cilia length. The affected cilia are " + str(cilia_length_outliers[0]))
            
            cilia_volume_outliers = self.detect_outlier(df, "Volume [micron^3]")
            if(len(cilia_volume_outliers[0]) != 0):  
                self.errors.append("The cilia in group " + group.getName() + " contains " + str(len(cilia_volume_outliers[0])) + " outliers in cilia length. The affected cilia are " + str(cilia_volume_outliers[0]))

            high_branched_cilia = self.compare_to_threshold(df, "# branches (quality parameter)", self.thresholds["max branches"], ">")
            if(len(high_branched_cilia) != 0):
                self.errors.append("The group " + group.getName() + " contains " + str(len(high_branched_cilia)) + " exceptionally branched cilia. The affected cilia are " + str(high_branched_cilia))

            cilia_branch_outliers = self.detect_outlier(df, "# branches (quality parameter)")
            if(len(cilia_branch_outliers[0]) != 0):  
                self.errors.append("The cilia in group " + group.getName() + " contains " + str(len(cilia_branch_outliers[0])) + " outliers in cilia length. The affected cilia are " + str(cilia_branch_outliers[0]))

            cilia_maxspan_outliers = self.detect_outlier(df, "Maximum span [micron]")
            if(len(cilia_maxspan_outliers[0]) != 0):  
                self.errors.append("The cilia in group " + group.getName() + " contains " + str(len(cilia_maxspan_outliers[0])) + " outliers in cilia length. The affected cilia are " + str(cilia_maxspan_outliers[0]))

            averageInt_outliers_cilia = self.detect_outlier(df, "average A intensity of the 10% of voxels with highest A intensity")
            if(len(averageInt_outliers_cilia[0]) != 0):
                self.errors.append("The group " + group.getName() + " contains " + str(len(averageInt_outliers_cilia[0])) + " exceptionally branched cilia. The affected cilia are " + str(averageInt_outliers_cilia[0]))

            underexposed_cilia_px = self.compare_to_threshold(df, "minimum intensity (in reconstruction channel)", 0, "=")
            if(len(underexposed_cilia_px) != 0):
                self.errors.append("The group " + group.getName() + " contains " + str(len(underexposed_cilia_px)) + " cilia that have underexposed pixels as part of their mask. The affected cilia are cilia " + str(underexposed_cilia_px))

            avInt_variance_outliers_cilia = self.detect_outlier(df, "average intensity (in reconstruction channel)")
            if(len(avInt_variance_outliers_cilia[0]) != 0):
                self.errors.append("The group " + group.getName() + " contains " + str(len(avInt_variance_outliers_cilia[0])) + " exceptionally branched cilia. The affected cilia are " + str(avInt_variance_outliers_cilia[0]) + ". Please verify the segmentation to prevent accidentally fused cilia.\n")
        
        if len(self.errors) > 0:
            self.errors.append("Please verify the segmentation to prevent accidentally fused cilia.\n")
        flattened_ciliaQ_versions = [item for sublist in ciliaQ_versions_all_groups for item in sublist]
        if (len(set(flattened_ciliaQ_versions)) > 1): 
            self.errors.append("You are about to pool data from different CiliaQ-versions. Since this is not advisable, please recheck your ciliaQ analysis. \n")
        self.raiseError(self.errors)
        self.write_qc_summary()

    def raiseError(self, text):
        root = Tk()
        root.lift()
        root.title("Raised Errors:")
        frameForMessage = ttk.Frame(root, padding = 10)
        frameForMessage.pack(fill="x", expand=True)
        for error in text:
            Label(frameForMessage, text = error, bg="white", anchor="w").pack()
        ButtonSelectFolder = Button(root, text = "Ignore", command= root.destroy, font=("Helvetica", 14, "bold"), bd=2, relief="raised", wraplength = 1000).pack()
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
    
    def write_qc_summary(self):
        self.write_path = self.write_path + "/qc_summary.txt" 
        with open(self.write_path, "w") as file: 
            file.write(f"This is a summary of all quality control errors that were arised on {self.time}:\n")
            file.write("\n")
            writer = csv.writer(file)
            if len(self.errors) == 0:
                csv.write(f"The analysis did not raise any quality control issues.")
            else:
                for error in self.errors:
                    writer.writerow([error])

