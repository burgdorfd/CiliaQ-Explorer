import os
import pandas as pd
import functions
from datagroup import DataGroup
from plotter import Plotter
from qc import QC
from setup import Setup
from stat_analysis import Statistics
from measurement import Measurements
from datetime import datetime

if __name__ == "__main__":
    pass

class Explorer:
    def __init__(self):
        self.all_measurements = []
        self.directories = []
        self.measurement_selection = []
        self.stat_selection = []
        self.perform_stat = ""
        self.metafileDir = None
        self.all_groups = []
        self.all_data = pd.DataFrame()
        self.statistics = None
        self.stat_list = []
        self.m = None
        self.setup = None
        self.time = datetime.now().strftime("%Y-%m-%d %H:%M:%S") 
        self.common_columns = None
        self.des_stat = None
        self.save_SVG = False
        self.save_PNG = False
        self.save_directory = None
        self.qc = None

    def setup_explorer(self):
        self.m = Measurements()
        self.all_measurements = self.m.get_measurement_list() 
        self.setup = Setup(self.all_measurements, self.m.get_measurement_category())
        self.setup.welcome_window()
        self.directories = self.setup.getDirectories()
        self.save_directory = f"{self.directories[0]}/Analysis"
        os.makedirs(self.save_directory, exist_ok=True)

        self.save_directory_plots = f"{self.save_directory}/Plots"
        os.makedirs(self.save_directory_plots, exist_ok=True)

        self.save_directory_plot_summary = f"{self.save_directory}/Data"
        os.makedirs(self.save_directory_plot_summary, exist_ok=True)

        self.measurement_selection, self.stat_selection, self.perform_stat, self.metafileDir, self.save_SVG, self.save_PNG = self.setup.setup_window()        
    
    def import_data(self):
        self.all_groups, self.all_data, self.common_columns = functions.import_ciliaQ_data(self.directories, self.all_measurements)
        neglected_measurements = [measurement for measurement in self.measurement_selection if measurement not in self.common_columns] 
        for neglected_measurement in neglected_measurements:
            print(f"Parameter {neglected_measurement} was not found in all replicates that were feeded to the analysis. Statistics and Plotting is not possible.")    
        self.measurement_selection = [measurement for measurement in self.measurement_selection if measurement in self.common_columns]
    
    def perform_statistics(self):
        self.statistics = Statistics(self.all_groups, self.all_data, self.measurement_selection) 
        self.des_stat = self.statistics.calculate_des_Statistics()
        measurement_selection_temp = functions.trim_Group_Replicate(self.measurement_selection)
        for measurement in measurement_selection_temp:
            isolate_measurement = pd.DataFrame()
            try:
                for group in self.all_groups:
                    data = group.getConcatenatedData()
                    if measurement in data:
                        temp_df = data[[measurement]].copy()
                        temp_df.columns = [group.getName()]
                        isolate_measurement = pd.concat([isolate_measurement, temp_df], axis=1)
                    else: 
                        print(f"error, measurement {measurement} not found in group {group.getName()}" )
                if(self.stat_selection == "automated statistics (see publication)"):
                    stat = self.statistics.hypothesis_testing(isolate_measurement, 0.05, measurement)
                    self.stat_list.append(stat)
                
            except UnboundLocalError:
                print(f"error, measurement {measurement} not found")
            except KeyError:
                print(f"error, measurement {measurement} not found")

    def perform_PCA(self):
        self.statistics.perform_PCA()
    
    def perform_UMAP(self):
        self.statistics.perform_UMAP()

    def plot_superplot(self):
        if "Replicate" in self.measurement_selection:
            self.measurement_selection.remove("Replicate")
        if "Group" in self.measurement_selection:
            self.measurement_selection.remove("Group")
        for measurement in self.measurement_selection:
            isolate_measurement = pd.DataFrame()
            plotter = Plotter(colors = [])
            try:
                for group in self.all_groups:
                    groupData = group.getConcatenatedData()
                    isolate_measurement = pd.concat([isolate_measurement, groupData[measurement]], axis=1)
                    isolate_measurement.rename(columns={measurement: group.getName()}, inplace = True)
                isolate_measurement = pd.melt(isolate_measurement)
                plotter.superplot5(self.all_data, measurement, self.m.get_measurement_dic(), self.stat_list, self.save_directory_plots, self.des_stat, self.save_SVG, self.save_PNG)
            except UnboundLocalError:
                print(f"error, measurement {measurement} not found")
            except KeyError:
                print(f"error, measurement {measurement} not found")

    def get_metadata(self):
        if self.metafileDir is not None:
            with open(self.metafileDir, "r", encoding="utf-8") as f:
                lines = [line.strip().split(";") for line in f]
                if (len(lines)) != 0:
                    self.research_group = str(lines[9][2])
                    print("this is research group:" + str(self.research_group))
                    self.researcher = str(lines[10][2])
                    self.exp_start = str(lines[11][2])
                    self.exp_title = str(lines[12][2])
                    self.used_SOP = str(lines[24][2])

                    self.model_organism = str(lines[16][2]) 
                    self.model_line_name = str(lines[17][2])
                    self.model_genotype = str(lines[19][2])
                    
                    self.microscopy_type = str(lines[37][2])
                    self.microscopy_name = str(lines[38][2])
                    self.magnification = str(lines[39][2])
                    self.num_aperture = str(lines[40][2])

                    self.group_names = [str(lines[30][i]) for i in range(2, len(lines[30])) if lines[30][i] is not None]
                    self.number_of_exp_groups = str(lines[26][2])
                    self.group_treatments = [str(lines[31][i]) for i in range(2, len(lines[31])) if lines[31][i] is not None]
                    self.group_duration = [str(lines[32][i]) for i in range(2, len(lines[32])) if lines[32][i] is not None]
                    self.group_concentrations = [str(lines[33][i]) for i in range(2, len(lines[33])) if lines[33][i] is not None]

                    self.immunostaining = True if (str(lines[42][2]) == "yes" or str(lines[42][2]) == "Yes" or str(lines[42][2]) == "YES") else False 
                    self.transgene = True if (str(lines[43][2]) == "yes" or str(lines[43][2]) == "Yes" or str(lines[43][2]) == "YES") else False
                    self.fixation = str(lines[40][2]) if (str(lines[40][2]) != "No" or str(lines[40][2]) is not None) else "no fixation used"
                    self.exp_replicate_names = [str(lines[i][2]) for i in range(64, len(lines)) if lines[i][2] is not None]

    def write_for5547_summary(self):
        self.get_metadata()
        write_path = self.save_directory_plot_summary + "/FOR5547_summary.txt"

        with open(write_path, "w") as csv: 
            csv.write("This is a summary of the metadata provided by a MFOR5547-metafile. This file serves serves to ensure reproducability across experiments and analysis runs.\n")
            csv.write("\n")
            if self.metafileDir is not None:
                csv.write("\n")
                
                csv.write(f"Experiment Name:\t{self.exp_title}\n")
                csv.write(f"Research group:\t{self.research_group}\n")
                csv.write(f"Researcher:\t{self.researcher}\n")
                csv.write(f"Experiment started on:\t{self.exp_start}\n")
                csv.write(f"used SOP:\t{self.used_SOP}\n\n")

                csv.write(f"Microscopy Type:\t{self.microscopy_type}\n")
                csv.write(f"Microscopy Name:\t{self.microscopy_name}\n")
                csv.write(f"Magnification:\t{self.magnification}\n")
                csv.write(f"Numerical Aperture:\t{self.num_aperture}\n\n")
                
                csv.write(f"Immunostaining:\t{self.immunostaining}\n")
                csv.write(f"Transgene:\t{self.transgene}\n")
                csv.write(f"Fixation:\t{self.fixation}\n\n")

                csv.write(f"Model Organism:\t{self.model_organism}\n")
                csv.write(f"Line Name:\t{self.model_line_name}\n")
                csv.write(f"genotype (if applicable):\t{self.model_genotype}\n\n")
                
                csv.write(f"Number of exp. conditions:\t{self.number_of_exp_groups}\n")
                csv.write("Treatments:\t" + "\t".join(self.group_treatments) + "\n")
                csv.write("Drug concentration:\t" + "\t".join(self.group_concentrations) + "\n")
                csv.write("Treatment duration:\t" + "\t".join(self.group_duration) + "\n\n")

    def perform_qc(self):
        self.qc = QC(self.all_groups, self.save_directory_plot_summary, self.time, self.setup.get_root())
        self.qc.perform_qc()

    def write_summary(self):
        self.get_metadata()
        write_path = self.save_directory_plot_summary + "/summary.txt"

        with open(write_path, "w") as csv:   
            csv.write(f"This is a summary of all replicates that were pooled on {self.time}:\n")
            csv.write("\n")

            if self.metafileDir is not None:
                csv.write("\n")
                csv.write(f"Experiment Name:\t{self.exp_title}\n")
                csv.write(f"Research group:\t{self.research_group}\n")
                csv.write(f"Researcher:\t{self.researcher}\n")
                csv.write(f"Experiment started on:\t{self.exp_start}\n")
                csv.write(f"used SOP:\t{self.used_SOP}\n\n")

                csv.write(f"Microscopy Type:\t{self.microscopy_type}\n")
                csv.write(f"Microscopy Name:\t{self.microscopy_name}\n")
                csv.write(f"Magnification:\t{self.magnification}\n")
                csv.write(f"Numerical Aperture:\t{self.num_aperture}\n\n")
                
                csv.write(f"Immunostaining:\t{self.immunostaining}\n")
                csv.write(f"Transgene:\t{self.transgene}\n")
                csv.write(f"Fixation:\t{self.fixation}\n\n")

                csv.write(f"Model Organism:\t{self.model_organism}\n")
                csv.write(f"Line Name:\t{self.model_line_name}\n")
                csv.write(f"genotype (if applicable):\t{self.model_genotype}\n\n")
                
                csv.write(f"Number of exp. conditions:\t{self.number_of_exp_groups}\n")
                csv.write("Treatments:\t" + "\t".join(self.group_treatments) + "\n")
                csv.write("Drug concentration:\t" + "\t".join(self.group_concentrations) + "\n")
                csv.write("Treatment duration:\t" + "\t".join(self.group_duration) + "\n\n")
 

        for group in self.all_groups:
            for rep in group.get_Replicate_Data():
                with open(write_path, "a") as csv:

                    csv.write("\n")
                    csv.write(f"group:\t{group.getName()}\n")
                    csv.write(f"Replicate:\t{rep.get_id()}\n")
                    csv.write(f"image name\t{rep.get_image_name()}\n")
                    csv.write(f"CiliaQ version\t{rep.get_ciliaQ_version()}\n")
                    csv.write("\n")

                    rep.get_data().to_csv(write_path, mode = "a", sep="\t", index = False)
                    print(rep.get_data())
                    
                    with open(write_path, "a") as csv:
                        csv.write("\n")
    
    def write_config(self):
        write_path = self.save_directory_plot_summary + "/config.txt"
        with open(write_path, "w") as csv: 
            csv.write("This is a config-file for the CiliaQ Analyzer:\n")
            csv.write("\n")
            for directory in self.directories:
                csv.write(f"{os.path.abspath(directory)} \t")
            for measurement in self.measurement_selection:
                csv.write(f"{measurement} \t")

    #def plot_single_Reps(self):
        #Plot single reps
        #colors = ["blue","orange","cyan","red","green","yellow",]
        #measurement_selection_temp = functions.trim_Group_Replicate(measurement_selection)
        #print(measurement_selection_temp)
        #for measurement in measurement_selection_temp:
        #    isolate_measurement = pd.DataFrame()
        #    plotter = functions.Plotter(colors)
        #    try:
        #        for group in all_groups:
        #            replicates = group.get_Replicate_Data()
        #            i = 1
        #            for rep in replicates:
        #                if measurement in rep.get_data():
        #                    temp_df = rep.get_data()[[measurement]].copy()
        #                    temp_df.columns = ["value"]
        #                    temp_df["variable"] = f"{group.getName()} - Rep. {i+1}"
        #                    isolate_measurement = pd.concat([isolate_measurement, temp_df], axis=0)
        #                    i = i + 1
        #                else: 
        #                    print(f"error, measurement {measurement} not found in group {group.getName()}, rep. {i}" )
        #                    i = i + 1
        #        plotter.plotBoxplot(isolate_measurement, measurement, m.get_measurement_dic())
        #        #plotter.plotViolinplot(isolate_measurement, measurement, m.get_measurement_dic())
        #    except UnboundLocalError:
        #        print(f"error, measurement {measurement} not found")

