import pandas as pd
from tkinter import *
from tkinter import filedialog, ttk
from datagroup import DataGroup

def import_ciliaQ_data(directories, column_names):
    all_groups = []
    
    for directory in enumerate(directories):
        group_name = slice_filename(directory[1])
        dataGroup = DataGroup(group_name, directory[1], column_names)
        dataGroup.import_group(directory[1], column_names, group_name)
        all_groups.append(dataGroup)

    all_replicates = []
    for group in all_groups:
        for rep in group.get_Replicate_Data():
            all_replicates.append(rep)
       
    common_columns = set(all_replicates[0].get_data().columns)
    for rep in all_replicates:
        common_columns = common_columns.intersection(rep.get_data().columns)
    common_columns = list(common_columns) 
    all_data = pd.DataFrame(columns=common_columns)
    for group in all_groups:
        for rep in group.get_Replicate_Data():
            rep.set_data(rep.get_data()[common_columns])

        group.pooledData = pd.DataFrame(columns=common_columns)
        for rep in group.get_Replicate_Data():
            trimmed_data = rep.get_data()
            group.pooledData = pd.concat([group.pooledData, trimmed_data], ignore_index=True)
        
        all_data = pd.concat([all_data, group.getConcatenatedData()])
    
    return all_groups, all_data, common_columns

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
