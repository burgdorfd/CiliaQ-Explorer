import os
import pandas as pd

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
        with open(self.file.path, "r", encoding = "ISO-8859-1", errors = "replace") as f:
            lines = [line.strip().split("\t") for line in f]
            if (len(lines)) != 0:
                self.image_name = str(lines[1][1])
                self.skl_gauss_sigma_XY = float(lines[20][1])
                self.skl_gauss_sigma_Z = float(lines[21][1])
                self.excluded_cilia = int(lines[18][3])
                self.analyzed_cilia = int(lines[18][5])
                self.calibration_XY = float(lines[8][1])
                self.calibration_Z = float(lines[8][3])


    def import_Data(self):
        try:
            dataset = pd.read_csv(
                    self.file.path, 
                    index_col= False, 
                    sep = "\t", 
                    skiprows = 28, 
                    usecols=range(60),
                    encoding = "ISO-8859-1",
                    encoding_errors = "replace"
                    ).dropna(axis="columns", how = "all")
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
    
    def set_data(self, data_param):
        self.data_df = data_param