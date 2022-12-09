#!/usr/bin/env python
# coding: utf-8

import tkinter as tk
from tkinter import ttk
from tkinter import filedialog as fd
from tkinter.messagebox import showerror, showwarning, showinfo
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg,
    NavigationToolbar2Tk
)

import os 
from os.path import join, basename
import time 
from copy import deepcopy
from collections import defaultdict

import numpy as np
import pandas as pd

from panel_classes import *
from panel_functions import *
from unmix_functions import *

from IPython.display import clear_output

import json

from tqdm import tqdm

import datetime

CONFIG_PARAM = {}

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Panel Maker")

        # get the screen dimension
        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()
        # set window dimension 
        window_width = int(screen_width * 0.8)
        window_height = int(screen_height * 0.7)
        # find the center point
        center_x = int(screen_width/2 - window_width / 2)
        center_y = int(screen_height/2 - window_height / 2)
        # set the position of the window to the center of the screen
        self.geometry(f'{window_width}x{window_height}+{center_x}+{center_y}')

        CONFIG_PARAM["load_button_width"] = int(window_width*0.025)
        CONFIG_PARAM["spinbox_width"] = int(window_width*0.01)

class Model():
    def __init__(self,
                 file_input_paths,
                 file_inputs,
                 params,
                ):
        self.samples = None
        self.panels = None
        self.file_input_paths = file_input_paths
        self.file_inputs = file_inputs 
        self.params = params 
        self.panel_scores = []
        self.unmixed_samples = []
        self.panel_recalls = []
        

class View(ttk.Frame):
    def __init__(self, parent):
        super().__init__(parent)
        
        self.input_frame = Input_Frame(self)

        self.control_frame = Control_Frame(self)
        self.param_frame = Params_Frame(self.control_frame)
        self.run_frame = Run_Frame(self.control_frame)
        self.plot_frame = Plot_Frame(self)
        self.option_frame = Options_Frame(self.plot_frame)
        self.download_frame = Download_Frame(self.plot_frame)
        self.canvas_frame = Canvas_Frame(self.plot_frame)

        self.input_frame.grid(row = 0, column = 0, rowspan=2, sticky="news")
        self.control_frame.grid(row = 0, column =1)
        self.param_frame.grid(row = 0, column = 0)
        self.run_frame.grid(row = 0, column = 1)
        self.plot_frame.grid(row = 1, column =1)
        self.canvas_frame.grid(row = 0, column = 1, rowspan=2)
        self.option_frame.grid(row = 0, column = 0)
        self.download_frame.grid(row = 1, column = 0)
        
        sep = ttk.Separator(self.param_frame,orient='vertical')
        sep.grid(column = 6, row = 0, rowspan = 2, sticky ="n", padx= 10)
    
    def set_controller(self, controller):
        self.controller = controller
        self.input_frame.controller = controller
        self.param_frame.controller = controller
        self.run_frame.controller = controller 
        self.option_frame.controller = controller
        self.plot_frame.controller = controller
        self.download_frame.controller = controller


class Input_Frame(ttk.Frame):
    def __init__(self, parent):
        super().__init__(parent)
        
        self.controller = None 
        
        self.load_config_button = ttk.Button(self,
                                             text = "Load Config File",
                                             command = self.load_config_file,
                                             width = CONFIG_PARAM["load_button_width"])
        self.load_config_button.pack()
        self.config_file_label = ttk.Label(self)
        self.config_file_label.config(text = "Config File:")
        self.config_file_label.pack(anchor="w")
        #Abundace
        self.abundance_button = ttk.Button(self, 
                                          text="Load Marker Abundance File",
                                          command = lambda: self.button_load_file("Abundance"),
                                          width = CONFIG_PARAM["load_button_width"])
        
        self.abundance_button.pack()
        
        self.abundance_file_label = ttk.Label(self)
        self.abundance_file_label.config(text = "Abundance File:")
        self.abundance_file_label.pack(anchor="w")
        
        #Brightness
        self.brightness_button = ttk.Button(self, 
                                           text="Load Fluor Brightness File",
                                           command = lambda: self.button_load_file("Brightness"),
                                           width = CONFIG_PARAM["load_button_width"])
        self.brightness_button.pack()
        
        self.brightness_file_label = ttk.Label(self)
        self.brightness_file_label.config(text = "Brightness File:")
        self.brightness_file_label.pack(anchor="w")
        
        #Target Markers 
        self.target_button = ttk.Button(self, 
                                          text="Load Marker Target File",
                                          command = lambda: self.button_load_file("Target"),
                                          width = CONFIG_PARAM["load_button_width"])
        self.target_button.pack()
        
        self.target_file_label = ttk.Label(self)
        self.target_file_label.config(text = "Target File:")
        self.target_file_label.pack(anchor="w")
        #Spectra 
        self.spectra_button = ttk.Button(self, 
                                          text="Load Spectra File",
                                          command = lambda: self.button_load_file("Spectra"),
                                          width = CONFIG_PARAM["load_button_width"])
        self.spectra_button.pack()
        
        self.spectra_file_label = ttk.Label(self)
        self.spectra_file_label.config(text = "Spectra File:")
        self.spectra_file_label.pack(anchor="w")
        #Database 
        self.database_button = ttk.Button(self, 
                                          text="Load Database File",
                                          command = lambda: self.button_load_file("Database"),
                                          width = CONFIG_PARAM["load_button_width"])
        self.database_button.pack()
        
        self.database_file_label = ttk.Label(self)
        self.database_file_label.config(text = "Database File:")
        self.database_file_label.pack(anchor="w")
        
    def set_controller(self, controller):
        self.controller = controller
        
    def button_load_file(self, file_type):
        file_name = fd.askopenfilename()
        base_file_name = basename(file_name)        
        if file_type == "Abundance":
            self.load_input_file(file_type="Abundance",
                                path=file_name)
            
        elif file_type == "Brightness":
            self.load_input_file(file_type="Brightness",
                                path=file_name)
        
        elif file_type == "Target":
            self.load_input_file(file_type="Target",
                                path=file_name)
        
        elif file_type == "Spectra":
            self.load_input_file(file_type="Spectra",
                                path=file_name)
            
        elif file_type == "Database":
            self.load_input_file(file_type="Database",
                                path=file_name)
    
    def load_input_file(self, file_type, path):
        if self.controller:
            self.controller.load_input(file_type, path)
        
    def load_config_file(self):
        if self.controller:
            self.controller.load_config_file()


class Params_Frame(ttk.Frame):
    def __init__(self, parent):
        super().__init__(parent)    
        
        self.controller = None


        # Search Params
        self.intensity_weight = tk.StringVar()
        self.corr_weight = tk.StringVar()
        self.random_prob = tk.StringVar()
        self.risk_prob = tk.StringVar()
        self.frontier_size = tk.StringVar()
        self.branch_num = tk.StringVar()
         
        intensity_param_label = ttk.Label(self, text="Intensity Weight")
        corr_param_label = ttk.Label(self, text="Correlation Weight")
        random_param_label = ttk.Label(self, text="Random Prob")
        risk_param_label = ttk.Label(self, text="Risk Prob")
        frontier_param_label = ttk.Label(self, text="Frontier Size")
        branch_param_label = ttk.Label(self, text="Branch Num")
        
        intensity_param_label.grid(row=0, column=0)
        corr_param_label.grid(row=1, column=0)
        random_param_label.grid(row=0, column=2)
        risk_param_label.grid(row=1,column=2)
        frontier_param_label.grid(row=0, column=4)
        branch_param_label.grid(row=1,column=4)

        intensity_spin_box = ttk.Spinbox(
            self,
            width=CONFIG_PARAM["spinbox_width"],
            from_=0,
            to=1,
            textvariable=self.intensity_weight,
            wrap=True,
            increment=0.1,
            command = lambda: self.update_params("intensity_weight",\
                                                 float(self.intensity_weight.get())))
        corr_spin_box = ttk.Spinbox(
            self,
            width=CONFIG_PARAM["spinbox_width"],
            from_=0,
            to=1,
            textvariable=self.corr_weight,
            wrap=True,
            increment=0.1,
            command = lambda: self.update_params("corr_weight",\
                                                 float(self.corr_weight.get())))
        random_spin_box = ttk.Spinbox(
            self,
            width=CONFIG_PARAM["spinbox_width"],
            from_=0,
            to=1,
            textvariable=self.random_prob,
            wrap=True,
            increment=0.1,
            command = lambda: self.update_params("random_prob",\
                                                 float(self.random_prob.get())))
        risk_spin_box = ttk.Spinbox(
            self,
            width=CONFIG_PARAM["spinbox_width"],
            from_=0,
            to=1,
            textvariable=self.risk_prob,
            wrap=True,
            increment=0.1,
            command = lambda: self.update_params("risk_prob",\
                                                 float(self.risk_prob.get())))
        frontier_spin_box = ttk.Spinbox(
            self,
            width=CONFIG_PARAM["spinbox_width"],
            from_=0,
            to=10,
            textvariable=self.frontier_size,
            wrap=True,
            increment=1,
            command = lambda: self.update_params("frontier_size",\
                                                 int(self.frontier_size.get())))
        branch_spin_box = ttk.Spinbox(
            self,
            width=CONFIG_PARAM["spinbox_width"],
            from_=0,
            to=10,
            textvariable=self.branch_num,
            wrap=True,
            increment=1,
            command = lambda: self.update_params("branch_num",\
                                                 int(self.branch_num.get())))
        
        intensity_spin_box.grid(column=1, row=0)
        corr_spin_box.grid(column=1, row=1)
        random_spin_box.grid(column=3, row=0)
        risk_spin_box.grid(column=3, row=1)
        frontier_spin_box.grid(column=5, row=0)
        branch_spin_box.grid(column=5, row=1)
        
        sep = ttk.Separator(self,orient='horizontal')
        sep.grid(column = 0, row = 2, columnspan = 6, sticky ="ew", pady= 5)
        
         # Unmix Params
        self.sample_size = tk.StringVar()
        self.intensity_factor = tk.StringVar()
        self.auto_intensity = tk.StringVar()
        self.noise = tk.StringVar()
        self.positive_fraction = tk.StringVar()
         
        sample_size_label = ttk.Label(self, text="Sample Size")
        intensity_factor_label = ttk.Label(self, text="Intensity Factor")
        auto_intensity_label = ttk.Label(self, text="Autofluor Intensity")
        noise_label = ttk.Label(self, text="Noise")
        positive_fraction_label = ttk.Label(self, text="Positive Fraction")
        
        sample_size_label.grid(row=3, column=0)
        intensity_factor_label.grid(row=3, column=2)
        auto_intensity_label.grid(row=3, column=4)
        noise_label.grid(row=4,column=0)
        positive_fraction_label.grid(row=4, column=2)

        sample_size_spin_box = ttk.Spinbox(
            self,
            width=CONFIG_PARAM["spinbox_width"],
            from_=0,
            to=100000,
            textvariable=self.sample_size,
            wrap=True,
            increment=100,
            command = lambda: self.update_params("sample_size",\
                                                 float(self.sample_size.get())))
        intensity_factor_spin_box = ttk.Spinbox(
            self,
            width=CONFIG_PARAM["spinbox_width"],
            from_=0,
            to=100000,
            textvariable=self.intensity_factor,
            wrap=True,
            increment=100,
            command = lambda: self.update_params("intensity_factor",\
                                                 float(self.intensity_factor.get())))
        auto_intensity_spin_box = ttk.Spinbox(
            self,
            width=CONFIG_PARAM["spinbox_width"],
            from_=0,
            to=100000,
            textvariable=self.auto_intensity,
            wrap=True,
            increment=100,
            command = lambda: self.update_params("auto_intensity",\
                                                 float(self.auto_intensity.get())))
        noise_spin_box = ttk.Spinbox(
            self,
            width=CONFIG_PARAM["spinbox_width"],
            from_=0,
            to=1,
            textvariable=self.noise,
            wrap=True,
            increment=0.1,
            command = lambda: self.update_params("noise",\
                                                 float(self.noise.get())))
        positive_fraction_spin_box = ttk.Spinbox(
            self,
            width=CONFIG_PARAM["spinbox_width"],
            from_=0,
            to=1,
            textvariable=self.positive_fraction,
            wrap=True,
            increment=.1,
            command = lambda: self.update_params("positive_fraction",\
                                                 int(self.positive_fraction.get())))
        
        sample_size_spin_box.grid(column=1, row=3)
        intensity_factor_spin_box.grid(column=3, row=3)
        auto_intensity_spin_box.grid(column=5, row=3)
        noise_spin_box.grid(column=1, row=4)
        positive_fraction_spin_box.grid(column=3, row=4)
        
    def set_controller(self, controller):
            self.controller = controller 
        
    def update_params(self, param, value):
        if self.controller:
            self.controller.update_params(param, value)
        

class Run_Frame(ttk.Frame):
    def __init__(self, parent):
        super().__init__(parent) 
        
        self.controller = None 
        result_num_label = ttk.Label(self, text="#Panels to Search: ")
        unmix_num_label = ttk.Label(self, text="#Panels to Unmix: ")
        self.result_num = tk.StringVar(value=10)
        self.unmix_num = tk.StringVar(value=10)
        result_spin_box = ttk.Spinbox(
            self,
            width=CONFIG_PARAM["spinbox_width"],
            from_=1,
            to=100000,
            textvariable=self.result_num,
            wrap=True,
            command = lambda: self.update_params("result_num",int(self.result_num.get())))
        unmix_spin_box = ttk.Spinbox(
            self,
            width=CONFIG_PARAM["spinbox_width"],
            from_=1,
            to=100000,
            textvariable=self.unmix_num,
            wrap=True,
            command = lambda: self.update_params("unmix_num",int(self.unmix_num.get())))
            
        run_beam_button = ttk.Button(self, 
                                    text = "Search",
                                    command = lambda: self.search(int(self.result_num.get())))
        run_unmix_button = ttk.Button(self, 
                                     text = "Unmix",
                                     command = lambda: self.unmix(int(self.unmix_num.get())))
        
        result_num_label.grid(row = 0, column = 0)
        unmix_num_label.grid(row = 2, column = 0)
        result_spin_box.grid(row = 1, column = 0)
        unmix_spin_box.grid(row = 3, column = 0)
        run_beam_button.grid(row = 0, column = 1, rowspan=2, sticky="ewns")
        run_unmix_button.grid(row = 2, column = 1, rowspan=2, sticky="ewns")        
        
    def set_controller(self, controller):
        self.controller = controller 
        
    def search(self, result_num):
        print("Search")
        if self.controller:
            self.controller.search(result_num)
            
    def unmix(self, unmix_num):
        if self.controller:
            self.controller.unmix(unmix_num)
    
    def update_params(self, param, value):
        if self.controller:
            self.controller.update_params(param, value)
        
class Control_Frame(ttk.Frame):
    def __init__(self, parent):
        super().__init__(parent)
        self.columnconfigure(0, weight = 3)
        self.columnconfigure(1, weight = 2)        
                
class Options_Frame(ttk.Frame):
    def __init__(self, parent):
        super().__init__(parent)
        
        self.controller = None 
        
        intensity_eval_weight = tk.StringVar(value=0.)
        corr_eval_weight = tk.StringVar(value=0.)
        plot_approx_button = ttk.Button(self,                                         text="Plot Approx Score",                                        command = lambda: self.plot_approx(float(intensity_eval_weight.get()),                                                                           float(corr_eval_weight.get())))
        plot_approx_button.grid(row = 0, column = 0)
        
        intensity_eval_label = ttk.Label(self, text="Intensity Eval")
        correlation_eval_label = ttk.Label(self, text="Correlation Eval")
        intensity_eval_label.grid(row=1, column=0)
        correlation_eval_label.grid(row=3, column=0)
        

        intensity_eval_spin_box = ttk.Spinbox(
            self,
            width=CONFIG_PARAM["spinbox_width"],
            from_=0,
            to=1,
            textvariable=intensity_eval_weight,
            wrap=True,
            increment=0.1, 
            command = lambda: self.update_params("intensity_eval_weight",\
                                                 float(intensity_eval_weight.get())))
        corr_eval_spin_box = ttk.Spinbox(
            self,
            width=CONFIG_PARAM["spinbox_width"],
            from_=0,
            to=1,
            textvariable=corr_eval_weight,
            wrap=True,
            increment=0.1,
            command = lambda: self.update_params("corr_eval_weight",\
                                                 float(corr_eval_weight.get())))
        
        intensity_eval_spin_box.grid(row=2, column=0)
        corr_eval_spin_box.grid(row=4, column=0)

    def plot_approx(self, intensity_eval_weight, corr_eval_weight):
        if self.controller:
            self.controller.plot_approx(intensity_eval_weight, corr_eval_weight)
    
    def update_params(self, param, value):
        if self.controller:
            self.controller.update_params(param, value)
        
class Plot_Frame(ttk.Frame):
    def __init__(self, parent):
        super().__init__(parent)  
        
        self.columnconfigure(0,weight=1)
        self.columnconfigure(1,weight=2)
        

        
class Download_Frame(ttk.Frame):
    def __init__(self, parent):
        super().__init__(parent)  
        
        self.controller = None
        self.download_panel_button = ttk.Button(self, text = "Download Panels", command=self.download_panel)
        
        self.output_name_label = ttk.Label(self, text="Output File Name")
        self.output_name = tk.StringVar()
        self.output_name_entry = ttk.Entry(self, textvariable=self.output_name)

        
        
        self.download_panel_button.grid(column=0, row=2)
        self.output_name_entry.grid(column=0, row=1)
        self.output_name_label.grid(column=0, row=0)
        
    def download_panel(self):
        
        output_name = self.output_name.get()
        if not output_name:
            output_name = f"panel_{str(datetime.datetime.now())}".replace(" ", "-").replace(".", "-").replace(":", "-")
        self.controller.download_panel(output_name)


class Canvas_Frame(ttk.Frame):
    def __init__(self, parent):
        super().__init__(parent)  
        
        figure = Figure(figsize=(6, 4), dpi=100)
        figure_canvas = FigureCanvasTkAgg(figure, self)
        NavigationToolbar2Tk(figure_canvas, self)
        axes = figure.add_subplot()
        figure_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        axes.plot([1,2,3])      
    
    def plot(self,plot_type, data, **kwargs):
        # try: 
        for child in self.winfo_children():
            child.pack_forget()
        if plot_type == "Line":
            figure = Figure(figsize=(6, 4), dpi=100)
            figure_canvas = FigureCanvasTkAgg(figure, self)
            NavigationToolbar2Tk(figure_canvas, self)
            axes = figure.add_subplot()
            figure_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
            if "sort" in kwargs.keys():
                df = pd.DataFrame(data, columns = ["Value"]).reset_index()
                df = df.sort_values("Value", ascending = True).reset_index(drop=True)
                for i in range(len(df)):
                    axes.text(y = df.Value[i], x = i,                               s = "#"+str(df["index"][i]), size = 10)
                axes.plot(df.Value)      
        # except Exception as error: 
        #     print("error")


class Controller:
    def __init__(self, model, view):
        self.model = model 
        self.view = view 
    
    def load_config_file(self):
        try: 
            config_path = fd.askopenfilename()
            config = {}
            with open(config_path, 'r') as fp:
                config = json.load(fp)
            assert "Panel_Config" == list(config.keys())[0], "Not a Panel Config File"
            for pair in list(config.items())[1:6]:
                self.load_input(pair[0], pair[1])
            for pair in list(config.items())[6:]:
                self.update_params(pair[0], pair[1])

            self.view.param_frame.intensity_weight.set(config["intensity_weight"])
            self.view.param_frame.corr_weight.set(config["corr_weight"])
            self.view.param_frame.random_prob.set(config["random_prob"])
            self.view.param_frame.risk_prob.set(config["risk_prob"])
            self.view.param_frame.frontier_size.set(config["frontier_size"])
            self.view.param_frame.branch_num.set(config["branch_num"])
            
            self.view.param_frame.sample_size.set(config["sample_size"])
            self.view.param_frame.intensity_factor.set(config["intensity_factor"])
            self.view.param_frame.auto_intensity.set(config["auto_intensity"])
            self.view.param_frame.noise.set(config["noise"])
            self.view.param_frame.positive_fraction.set(config["positive_fraction"])
            
            self.view.run_frame.result_num.set(config["result_num"])
            self.view.run_frame.unmix_num.set(config["unmix_num"])
            
            file_name = basename(config_path)
            self.view.input_frame.config_file_label["text"] = "Config file: " + file_name

        except Exception as error:
            showerror(                      title = "Error",
                      message = error)
            return
            
        
    def load_input(self, file_type, path):
        base_file_name = basename(path)
        extension = base_file_name.split(".")[-1]
        if extension == "csv":
            df = pd.read_csv(path)
        elif extension == "xlsx":
            df = pd.read_excel(path)
        else:
            showerror(            title='Error',            message='FIle extensions must be .csv or .xlsx')
            return
        
        if file_type == "Abundance":
            if self.check_abundance_format(df):
                model.file_inputs[file_type] = df
                model.file_input_paths[file_type] = path 
                self.view.input_frame.abundance_file_label.config(                    text = "Abundance file: " + base_file_name)
        elif file_type == "Brightness":
            if self.check_brightness_format(df):
                model.file_inputs[file_type] = df
                model.file_input_paths[file_type] = path 
                self.view.input_frame.brightness_file_label.config(                    text = "Brightness file: " + base_file_name)
        elif file_type == "Target":
            if self.check_target_format(df):
                model.file_inputs[file_type] = df
                model.file_input_paths[file_type] = path 
                self.view.input_frame.target_file_label.config(                    text = "Target file: " + base_file_name)
        elif file_type == "Spectra":
            if self.check_spectra_format(df):
                model.file_inputs[file_type] = df
                model.file_input_paths[file_type] = path 
                self.view.input_frame.spectra_file_label.config(                    text = "Spectra file: " + base_file_name)
        elif file_type == "Database":
            if self.check_database_format(df):
                model.file_inputs[file_type] = df
                model.file_input_paths[file_type] = path 
                self.view.input_frame.database_file_label.config(                    text = "Database file: " + base_file_name)
    
    
    def check_abundance_format(self, df):
        try: 
            assert len(df.columns) == 2, "Abundance file must contains 2 columns"
            assert df.columns[0] == "Marker", "Abundance file: 1st column must be named 'Marker'"
            assert df.columns[1] == "Abundance", "Abundance file: 2nd column must be named 'Abundance'"
        
        except AssertionError as error:
            showerror(            title='Error',            message=error)
            return False 
        
        return True
            
    def check_brightness_format(self, df):
        try: 
            assert len(df.columns) == 2, "Brigtness file must contains 2 columns"
            assert df.columns[0] == "Fluor", "Brightness file: 1st column must be named 'Fluor'"
            assert df.columns[1] == "Brightness", "Brightness file: 2nd column must be named 'Brightness'"
        
        except AssertionError as error:
            showerror(            title='Error',            message=error)
            return False 
        return True
    def check_target_format(self, df):
        try: 
            assert len(df.columns) == 1, "Target file must contains 1 column"
            assert df.columns[0] == "Marker", "Target file: column must be named 'Marker'"
        
        except AssertionError as error:
            showerror(            title='Error',            message=error)
            return False
        return True
    
    def check_spectra_format(self, df):
        try: 
            pass
        except AssertionError as error:
            showerror(            title='Error',            message=error)
            return False
        return True
    
    def check_database_format(self, df):
        try: 
            pass
        except AssertionError as error:
            showerror(            title='Error',            message=error)
            return False
        return True
    
    def update_params(self, param, value):
        self.model.params[param] = value
        
    def search(self, panel_num):
        
        try:
            assert len(model.file_inputs.values()) == 5, "Missing input files"
            assert len(model.params.keys()) >= 6, "Missing parameters"

            params = model.params 
            abundance_df, brightness_df, target_df, spectra_df, database_df = model.file_inputs.values()

            database_df, marker_name_list, marker_abundance_list,             fluor_name_list, fluor_brightness_list, fluor_channel_values_list =             panel_data_preprocessing(abundance_df, brightness_df, target_df, spectra_df, database_df)

            start = time.time()
            self.afhulym = Fluor(name = "AfHuLym", brightness = 1,                            channel_values= spectra_df["AfHuLym"].values)
            real_marker_strand = create_marker_strand(marker_name_list, marker_abundance_list)
            real_fluor_strand = create_fluor_strand(fluor_name_list,                                                    fluor_brightness_list, fluor_channel_values_list)
            real_dna = create_dna(real_fluor_strand, real_marker_strand,                                   database_df["Fluor"].values, database_df["Marker"].values,                                 marker_abundance_list, fluor_brightness_list)
            results = constrained_beam_search(init_panel= real_dna, result_num= panel_num,                                               frontier_size=params["frontier_size"], branch_num=params["branch_num"],                                              random_prob=params["random_prob"], risk_prob=params["risk_prob"],                                               auto_fluor=self.afhulym, alpha = params["intensity_weight"], beta = params["corr_weight"])



            end = time.time()
            print('{:.4f} s'.format(end-start))
            panels_flattened = [i for b in map(lambda x:[x]                                                if not isinstance(x, list) else x,                                               results) for i in b]
            model.panels = panels_flattened

            model.panel_scores = self.calculate_approx_scores(model.panels, model.params["intensity_eval_weight"],                                                              model.params["corr_eval_weight"])
            
            model.panel_file = self.generate_panel_file()
            
        except Exception as error:
            showerror(                      title="Error",                      message=error)
            return

    def calculate_approx_scores(self, panels, intensity_eval_weight, corr_eval_weight):
        panel_scores = []
        for panel in panels:
            score = calculate_panel_approx_score(markers=panel.marker_history, fluors=panel.fluor_history,                                          auto_fluor=self.afhulym, intensity_eval_weight=intensity_eval_weight                                                 , corr_eval_weight=corr_eval_weight)
            panel_scores.append(score)
        return panel_scores
        
    def plot_approx(self, intensity_eval_weight, corr_eval_weight):
        try:
            assert isinstance(self.model.panel_scores, list), "No Panels Found"
        except AssertionError as error:
            showerror(            title='Error',            message=error)
            return
        
        model.panel_scores = self.calculate_approx_scores(model.panels, intensity_eval_weight, corr_eval_weight)
        self.view.canvas_frame.plot(plot_type = "Line",                                     data = self.model.panel_scores, 
                                    sort = True,
                                    index = np.arange(len(self.model.panel_scores)))
        
    def unmix(self, unmix_num):
        
        recall_results = []
        df_all_results = []
        all_samples =[]
        method = "Ridge"
        params = model.params 
        model.unmixed_samples = []
        model.panel_recalls = []
        abundance_df, brightness_df, target_df, spectra_df, database_df = model.file_inputs.values()
        panel_and_score = tuple(zip(model.panels, model.panel_scores))
        panel_and_score_sorted = sorted(panel_and_score, key = lambda x: x[1])
        panel_sorted = []
        
        for panel, score in panel_and_score_sorted:
            panel_sorted.append(panel)
        for i, panel in enumerate(panel_sorted[:unmix_num]):
            sample = Sample(ID = i, cell_num=params["sample_size"], dna = panel, channel_num = len(model.file_inputs["Spectra"]))
            # Add each flour to 20% cells in sample:
            for fluor in panel.fluor_history:
                # Add 20%
                chosen_cells = np.random.choice(sample.cells,                                                 int(params["positive_fraction"]*params["sample_size"]), replace=False)
                for cell in chosen_cells:
                    intensity_factor = np.random.normal(0,1)
                    intensity_factor = np.exp(intensity_factor/2) * params["intensity_factor"]
                    cell.add_fluor(fluor = fluor, intensity_factor=intensity_factor)
                    sample.components[fluor.name].append(cell.id)
            # Add Autofluor
            for cell in sample.cells:
                auto_fluor = np.random.normal(0,1)
                auto_fluor = np.exp(auto_fluor/2) * params["auto_intensity"]
                cell.add_fluor(fluor = self.afhulym, intensity_factor = params["auto_intensity"])

            spectra_panel = spectra_df[[fluor.name for fluor in panel.fluor_history]].copy()
            df_result = pd.DataFrame(columns= spectra_panel.columns.to_list() + ["id"]) 
            for i in tqdm(range(len(sample.cells))):
                measurement = sample.cells[i].emit_light(noise=params["noise"])
                unmixed_read = unmix(measurement, spectra_panel, method , alpha = .1, max_iter = 5000)
                unmixed_read["id"] = sample.cells[i].id
                df_result = df_result.append(unmixed_read, ignore_index=True)
            all_samples.append(sample)
            df_all_results.append(df_result)
    
            #Calculate Recall
            recall_dict = {}
            for fluor in df_result.columns[:-1]:
                positive_population = set(df_result.sort_values(fluor, ascending = False)["id"]                                          [:int(params["sample_size"]*params["positive_fraction"])].index.to_list())
                match_size = len(set(sample.components[fluor]).intersection(set(positive_population)))
                recall = match_size/(int(params["sample_size"]*params["positive_fraction"]))
                recall_dict[fluor] = recall
            recall_results.append(recall_dict)
        
            print(recall_results, len(recall_results))
            model.panel_recalls=recall_results
            model.unmixed_samples = df_all_results
            model.samples = all_samples
    
    def generate_panel_file(self):
        panels = self.model.panels
        df = pd.DataFrame()
        count=0
        for panel in panels:
            count+=1
            fluors = []
            markers = []
            bond_history = panel.bond_history 
            for bond in bond_history:
                fluors.append(bond.fluor_end.name)
                markers.append(bond.marker_end.name)
            df.insert(loc = len(df.columns), column = "Marker"+str(count) , value = markers)  
            df.insert(loc = len(df.columns), column = "Fluor"+str(count), value = fluors)
        return df 
    
    def download_panel(self, output_name):
        df = self.generate_panel_file() 
        assert len(df) > 0, "No panels generated"
        save_dir = fd.askdirectory()
        df.to_csv(join(save_dir, f"{output_name}.csv"), index=None)


if __name__ == "__main__":
    app = App()
    view = View(app)
    model = Model({},{},{})
    controller = Controller(model, view)
    view.set_controller(controller)
    view.pack()
    app.mainloop()

