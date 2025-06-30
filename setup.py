import os
import pandas as pd
import tkinter as tk
from tkinter import *
from tkinter import filedialog, ttk, messagebox
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
import re
from PIL import Image, ImageTk
import webbrowser

class Setup():

    def __init__(self, all_measurements_param, measurement_category_param):
        self.all_measurements = all_measurements_param
        self.measurement_category = measurement_category_param
        self.measurement_selection = []
        self.stat_selection = "automated statistics (see publication)"
        self.perform_stat = False
        self.metafileDir = None
        self.save_SVG = False
        self.save_PNG = False
        self.root = tk.Tk()
        self.root.geometry("0x0")
        self.root.withdraw()
        self.root.lift()
        self.root.update()
        self.root.attributes('-topmost', True)

    def close_window(self, window):
        window.grab_release()
        window.destroy()
        self.root.deiconify()
        self.root.withdraw()

    def get_root(self):
        return self.root
    
    def getDirectories(self):
        directories = []
        self.root.withdraw()
        def selectFolder(directoryList, var):
            filepath = filedialog.askdirectory(parent = window)
            if filepath:
                directoryList.append(filepath)
                var.set("\n".join(directoryList))
                #label.configure(text = "\n".join(directoryList), background = "white")
                #label.update_idletasks()

        window = tk.Toplevel(self.root)
        window.geometry("600x400")
        window.lift()
        window.title = ("Please select folders")

        window.columnconfigure((0,1), weight = 1)
        window.rowconfigure((0,2), weight = 1)
        window.rowconfigure(1, weight = 10)

        frame = Frame(window, bg= "white", relief="sunken", bd=2)
        frame.grid(row = 1, column = 0, sticky = "news", columnspan=2)
        Label(window, text = "Selected directories").grid(row = 0, column = 0,columnspan=2, sticky = "nw")
        display_var = StringVar()
        display = Label(frame, textvariable=display_var)
        display.grid(pady = 10, sticky = "n")

        ButtonSelectFolder = Button(window, text = "Select Folder", command= lambda: selectFolder(directories, display_var), font=("Helvetica", 14, "bold"), bd=2, relief="raised").grid(row = 2, column = 0, sticky = "se")
        ButtonSelectFolder = Button(window, text = "Apply", command = lambda: self.close_window(window), font=("Helvetica", 14, "bold"), bd=2, relief="raised").grid(row = 2, column = 1, sticky = "se")
        
        #window.lift()
        window.protocol("WM_DELETE_WINDOW", self.close_window) # intercept close button
        window.transient(self.root)   # dialog window is related to main
        window.wait_visibility()
        window.grab_set()
        window.wait_window()
        #self.root.deiconify()
        return directories
    
    def welcome_window(self):
        window = tk.Toplevel(self.root)
        window.lift()
        window.resizable(False, False)
        window.title("Welcome to the CiliaQ Analyzer")
        window.attributes('-topmost', True)

        def open_link(event):
            webbrowser.open_new("10.1140/epje/s10189-021-00031-y")

        tk.Label(window, text =  "CiliaQ Analyser, Version 1.1.0, 2025", font=("Helvetica", 16, "bold"), anchor='w').pack()
        tk.Label(window, text =  "Welcome to the CiliaQ Analyser. This notebook is an addition to the CiliaQ Plugin for Fiji ImageJ (Hansen JN et al, 2021) \n \n", font=("Helvetica", 12), wraplength = 600, justify = "left", anchor='w').pack()
        tk.Label(window, text = "This notebook provides a step-by-step pipeline that can process CiliaQ-derived analysis results from different replicates and experimental conditions. This notebook allows pooling CiliaQ-derived CQ files from different experimental replicates and conditions, performing quality control and statistical analysis of the data and plotting all ciliary parameters in a superplot-format. \n \n", wraplength = 600, justify = "left", font=("Helvetica", 12), anchor='w').pack()
        tk.Label(window, text =  "In the following, you will be guided through the workflow. For more detailed instructions, read the jupyter notebook or refer to our protocol paper (Burgdorf et al., 2025)",  wraplength = 600, justify = "left", font=("Helvetica", 12), anchor='w').pack()
        Button(window, text = "Let's start!", command = lambda: self.close_window(window), font=("Helvetica", 14, "bold"), bd=2, relief="raised").pack()

        window.lift()
        window.grab_set()
        window.wait_window()

        # Re-enable root window after closing the setup window
        self.root.deiconify()

    def setup_window(self):
        
        window = tk.Toplevel(self.root)
        window.resizable(True, True)
        window.title("Analyzer setup:")
        window.lift()
        
        def select_stat(event): 
            self.stat_selection = drop_down.get()

        def selectMetafile():
            file_path = filedialog.askopenfilename(
            title="Select a metasheet (optional, but recommended)",
            filetypes=[("CSV Files", "*.csv"), ("All Files", "*.*")]
            )
            self.metafileDir = file_path

        def save_Image_Preference_SVG(state):
            if state.get() == 1:
                self.save_SVG = True
            if state.get() == 0:
                self.save_SVG = False
                
        def save_Image_Preference_PNG(state):
            if state.get() == 1:
                self.save_PNG = True
            if state.get() == 0:
                self.save_PNG = False


        tk.Label(window, text =  "Include measurements into analysis:", font=("Helvetica", 12)).pack()
        Button(window, text =  "Select measurements", font=("Helvetica", 12, "bold"), bd=2, relief="raised", command = lambda: self.set_measurement_window(window)).pack()

        tk.Label(window, text =  "Upload Metafile", font=("Helvetica", 12)).pack(pady= 20)
        Button(window, text = "Select Metafile", command= lambda : selectMetafile(), font=("Helvetica", 14, "bold"), bd=2, relief="raised").pack()

        tk.Label(window, text =  "Method for hypothesis testing", font=("Helvetica", 12)).pack(pady = 10)
        options = ["No statistics", "automated statistics (see publication)"]
        drop_down = ttk.Combobox(window, values = options, state = "readonly")
        drop_down.pack()
        drop_down.current(1)
        drop_down.bind("<<ComboboxSelected>>", select_stat)
        

        result_svg = tk.IntVar()
        Checkbutton(window, text= "Save plots as .SVG",  variable = result_svg, onvalue= 1, offvalue=0, command = lambda r = result_svg: save_Image_Preference_SVG(r)).pack(padx = 10, pady = 10)

        result_png = tk.IntVar()
        Checkbutton(window, text= "Save plots as .PNG",  variable = result_png, onvalue= 1, offvalue=0, command = lambda r = result_png: save_Image_Preference_PNG(r)).pack(padx = 10)

        Button(window, text = "Apply", command = lambda: self.close_window(window), font=("Helvetica", 14, "bold"), bd=2, relief="raised").pack(padx = 10, pady = 10)
        #window.lift()        
        window.grab_set()
        window.wait_window()
        # Re-enable root window after closing the setup window
        self.root.deiconify()
        return self.measurement_selection, self.stat_selection, self.perform_stat, self.metafileDir, self.save_SVG, self.save_PNG
    
    def set_measurement_window(self, parent):
        #self.all_measurements = self.all_measurements[0:5] #only for testing purposes
        window = tk.Toplevel(parent)
        window.resizable(True, True)
        window.title("Please select Measurements")
        self.root.withdraw()
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
        category_columns = {"Channel A/B:": 1, "Morphology": 0, "Spatial Parameters":0, "Cilia Markers": 0}
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
        bottom_row = max(row_count[0], row_count[1]) + 2
        
        Button(window, text = "Apply", command = lambda: self.close_window(window), font=("Helvetica", 14, "bold"), bd=2, relief="raised").grid(row=bottom_row, column = 1, columnspan = 1, padx = 10, sticky = "s")
        Button(window, text =  "Select all", font=("Helvetica", 14, "bold"), bd=2, relief="raised", command = choose_all).grid(row=bottom_row, column = 0, columnspan = 1, padx = 10, sticky = "s")
        parent.update_idletasks()
        window.grab_set()
        #window.lift()
        window.focus_set()
        self.root.deiconify()
        #window.wait_window()
