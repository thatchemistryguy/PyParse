"""

Copyright 2023 GlaxoSmithKline Research & Development Limited

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Authors: Joe Mason, Francesco Rianjongdee, Harry Wilders, David Fallon
"""


import os
import sys
from os import path

import time
import logging
import argparse
import re
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
import pandas as pd
import math
from statistics import mean
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw
import zipfile
from jinja2 import Environment, PackageLoader
import numpy as np
import matplotlib.cm as cm
import glob

import getWatersData
import cpAssignment
import genOutput
import plotting

def getUserReadableWell(wellno, plate_col_no):
    """
    Converts the well as a number into a user-friendly string,
    e.g. well 11 becomes "B5" for a 4*6 well plate. 

    :param wellno: An integer representing a specific well on the plate
    
    :return: A string representing a specific well on the plate
    """
    
    rowVal = math.floor((wellno-1) / plate_col_no)
    colVal = (wellno) % plate_col_no
    if colVal == 0:
        colVal = plate_col_no
    
    label = f'{chr(ord("@")+(rowVal)+1)}{colVal}'
    return label

def generateMol(smiles, name, save_dir):
    """
    Generates a 2D rendering of the given structure and saves 
    it as a .png

    :param smiles: a string (SMILES) of the compound
    :param name: a string corresponding to name of that compound
    :param save_dir: a string corresponding to the output directory
    
    :return: A .png rendering of the given compound
    """

    mol = Chem.MolFromSmiles(smiles.strip())
    _discard = AllChem.Compute2DCoords(mol)
    Draw.MolToFile(mol, f'{save_dir}structures/{name}.png', size=(200, 150))

def buildHTML(save_dir, cpTable, analysis_name, times = {}):
    """
    Build a HTML output file using jinja2 and a html_template
    that is stored in the directory "templates". 
    
    :param save_dir: A string designating the output directory
    :param cpTable: Pandas datatable containing all information on 
                        the compounds used for analysis
    :param all_compounds: a list of all compound names
    :param impurities: a list of all impurity names
    :param analysis_name: User provided name for the analysis
    :param times: Optional parameter of a list of floats related to processing time
                    for each step of the analysis. 
    :return: HTML file saved to save_dir 
    """


    
    env = Environment(
        loader = PackageLoader("PyParse", "templates")
    )

    template = env.get_template("html_template.html")

    cptablerows = cpTable.to_dict(orient="records")
    

    with open(f'{save_dir}/html_output.html', "w") as fo:
        fo.write(template.render(
            cpnames = list(cpTable.loc[cpTable["type"] != "Impurity", "name"]),
            imp_no = len(cpTable.loc[cpTable["type"] == "Impurity"].index),
            cptablerows = cptablerows, 
            save_dir = save_dir,
            path = path,
            times = times,
            round = round,
            pt = options.plot_type, 
            analysis_name = analysis_name,
            options = vars(options)
            )
        )
    fo.close() 

def main():

    """ 
    Provide an .rpt file and a csv containing the compounds in the 
    plate to analyse the wells and return an output containing plots of compounds,
    key wells, and multiple output tables. 
    """
    
    usage = "Automated LCMS Analysis and Data Extraction"

    #Start timer
    time_start = time.perf_counter()

    #Setup argparse and load default values
    parser = argparse.ArgumentParser(description = usage,
                                     formatter_class=argparse.RawTextHelpFormatter,)
    parser.set_defaults(
                        validate = "True",
                        verbose = "False",

                        analysis_name = "PyParse Analysis",

                        mass_abs_tol = 0.5,
                        time_abs_tol = 0.025,
                        uv_abs_tol = 10,

                        filter_by_rt = [],


                        detector = "UV",

                        min_peak_area = 0.5,
                        min_massconf_threshold = 10,
                        min_uv_threshold = 20,

                        min_no_of_wells = 5,

                        uv_match_threshold = 0.5,
                        uv_cluster_threshold = 0.5,
                        massconf_threshold = 0.5,
                        cluster_size_threshold = 0.8,

                        #Default values are set to 0, to easily monitor
                        #if the user has manually specified them.    
                        plate_col_no = 0,
                        plate_row_no = 0,

                        points_per_trace = 500,
                        
                        mass_or_area = "mass_conf",
                        plot_type = "Parea", 
                        calc_higherions = "True",
                        calc_boc = "True",
                        
                        gen_csv = "True",
                        gen_zip = "True",
                        gen_atable = "True", 

                        output = "output",

                        instrument = "Waters",

                        find_freq_imp = "True",
                        )
    
    parser.add_argument("input_data", help = "Input file/folder for analysis")
    
    parser.add_argument("input_csv", help="Input .csv file containing compounds in plate")
    
    parser.add_argument("-o", "--output", action="store", type=str, dest = "output",
                        help = "Location to store output.")
    
    parser.add_argument("-V", "--validate", action = "store", type=str, dest = "validate",
                       help = "True/False: Run hit validation processes. \n")

    parser.add_argument("-v","--verbose", action = "store", type = str, dest = "verbose",
                     help = "Set to 'True' to log verbosely.\n")

    parser.add_argument("-mat","--mass_abs_tol",action = "store", type = float, dest = "mass_abs_tol",
                     help = "Variable to determine how close a m/z should be to match a given mass for a structure.\n")

    parser.add_argument("-tat","--time_abs_tol",action = "store", type = float, dest = "time_abs_tol",
                     help = "Variable to determine how close a peak time should be to match the rest of a cluster of peaks.\n")

    parser.add_argument("-uat","--uv_abs_tol",action = "store", type = int, dest = "uv_abs_tol",
                     help = "Variable to determine how close the UV max should be to be counted in the same cluster.\n")

    parser.add_argument("-mpa","--min_peak_area",action = "store", type = float, dest = "min_peak_area",
                     help = "Variable to determine the minimum acceptable peak area given as a percentage of all peak areas.\n")

    parser.add_argument("-mmt","--min_massconf_threshold", action = "store",type = int,dest = "min_massconf_threshold",
                     help = "Variable to determine the lowest massConf at which a hit will be recorded.\n")

    parser.add_argument("-mut","--min_uv_threshold", action = "store", type = int, dest = "min_uv_threshold",
                     help = "Minimum threshold before a UV signal can actually be considered to be a maxima.\n")

    parser.add_argument("-umt","--uv_match_threshold", action = "store", type=float, dest = "uv_match_threshold",
                     help = "Variable to determine what fraction of UV maxima a hit needs to have to be deemed not orange.\n")

    parser.add_argument("-uct","--uv_cluster_threshold", action = "store",type = float, dest = "uv_cluster_threshold",
                     help = "Variable to determine how big a UV cluster must be compared to the largest UV cluster expressed as a fraction.\n")

    parser.add_argument("-mt","--massconf_threshold", action = "store", type=float, dest = "massconf_threshold",
                     help = "Variable for threshold at which a massConf is deemed too far below the mean for that cluster to be worthy of inclusion in cluster['green'].\n")

    parser.add_argument("-cst","--cluster_size_threshold", action = "store", type = float, dest = "cluster_size_threshold",
                     help = "Variable, defining the point at which a cluster is deemed too small to be valid, by comparison to the length of the largest cluster, as a percentage.\n")
    
    parser.add_argument("-r","--rows",action = "store", type=int, dest = "plate_row_no",
                     help = "Number of rows in the plate.\n")

    parser.add_argument("-c","--columns",action = "store", type = int, dest = "plate_col_no",
                     help = "Number of columns in the plate.\n")
    
    parser.add_argument("-now", "--min_no_of_wells", action="store", type=int, dest = "min_no_of_wells", 
                        help = "Minimum number of wells that have a hit before validation process is run.\n")
    
    parser.add_argument("-ppt", "--points_per_trace", action="store", type=int, dest = "points_per_trace", 
                        help = "Define data points per chromatogram as an integer.\n")
    
    parser.add_argument("-moa", "--mass_or_area", action="store", type=str, dest = "mass_or_area", 
                        help = "'mass_conf' or 'area': Choose how to determine which peak is chosen in the absence of validation.\n")
    
    parser.add_argument("-pt", "--plot_type", action="store", type=str, dest = "plot_type",
                        help = "Choose what to plot in the heatmap and piechart "
                               "from the following options: SMarea, Parea, STDarea, "
                               "corrSMarea, corrParea, corrSTDarea, P/SM+P, P/STD, "
                               "corrP/STD, corrP/SM+P.\n")
    
    parser.add_argument("-chi", "--calc_higherions", action="store", type=str, dest = "calc_higherions",
                        help = "Look for [M+2H]2+ and [M+3H]3+ to find hits and calculate the mass confidence"
                        " of a peak, True/False")
                        
    parser.add_argument("-cboc", "--calc_boc", action="store", type=str, dest = "calc_boc",
                        help = "Look for Boc degradation ions, True/False")
                        
    parser.add_argument("-g", "--generate_csv", action="store", type=str, dest = "gen_csv", 
                        help = "Choose to generate and save a CSV, True/False.\n")
    
    parser.add_argument("-z", "--generate_zip", action="store", type=str, dest = "gen_zip", 
                        help = "Choose to generate and save a zip file, True/False.\n")
    
    parser.add_argument("-a", "--generate_atable", action="store", type=str, dest = "gen_atable", 
                        help = "Choose to generate and save an analytical table, True/False.\n")
    
    parser.add_argument("-f", "--filter_by_rt", nargs = "+", action="store", type=str, dest = "filter_by_rt",
                        help = "Provide retention time ranges to ignore, in the form 'mintime-maxtime'"
                        " where the two times are separated by a hyphen. Multiple ranges (separated by a space) are excepted. ")

    parser.add_argument("-n", "--name", action="store", type=str, dest = "analysis_name", 
                        help = "Choose a name for the analysis.\n")

    parser.add_argument("-d", "--detector", action="store", type=str, dest = "detector",
                        help = "Choose which detector to use, UV or ELSD. Waters data only.")
    
    parser.add_argument("-ffi", "--find_freq_imp", action="store", type=str, dest = "find_freq_imp",
                        help = "Choose whether to search for and report frequently observed impurities. ")

    parser.add_argument("-i", "--instrument", action="store", type=str, dest = "instrument",
                        help = "Select the type of machine the data originated from, Waters, Agilent or Shimadzu")

    #Set options to global and parse arguments        
    times = {}
    global options
    options = parser.parse_args()

    #Ensure all required input options have a valid address
    root_names = [options.input_data, options.input_csv, options.output]
    
    for name in root_names:
        if name.startswith("/"):
            pass
        else:
            if name.startswith("./"):
                name.replace("./", os.getcwd())
            else:
                name = os.getcwd() + "/" + name    
    
    #Create the save directory if one isn't already present.
    save_dir = root_names[2]
    if not save_dir.endswith("/"):
        save_dir = save_dir + "/"
    
    error_msg = ""
    try:
        os.mkdir(save_dir)
    except OSError as error:
        error_msg = f'The directory "{save_dir}" already exists.'
	
    #Start logging, and log error_msg if one was created whilst making the output
    #directory. 
    if options.verbose == "True":
        logging.basicConfig(filename = f'{save_dir}/logfile.txt', level=logging.DEBUG)    
    else:
        logging.basicConfig(filename = f'{save_dir}/logfile.txt', level=logging.INFO)
    if error_msg != "":
        logging.debug(error_msg)
    
    #Boiler plate logs 
    logging.info(f'PyParse was run from the following directory: {os.getcwd()}.')
    logging.info(f'The output files were saved to {save_dir}.')

    #Generate a generic HTML output file so that PyParse displays
    #the log file even in cases of uncaught errors. 
    with open(f'{save_dir}html_output.html', 'w') as f:
        f.write(f'Unspecified error occurred. See logfile for details.')
        f.close()

    #make sub-directories to store all graphs and pictures of structures
    try: 
        os.mkdir(f'{save_dir}/graphs')
    except OSError as error:
        logging.debug("Graphs directory already exists.")
    try: 
        os.mkdir(f'{save_dir}/structures')
    except OSError as error:
        logging.debug("Structures directory already exists.")

        #Check to ensure the root names of the data file/folder and csv file are correct, i.e. they exist
    for name in [root_names[0], root_names[1]]:
        if not os.access(name,os.R_OK) or not os.path.exists(name):
            logging.error(f'Input file/folder {name} does not exist. Please use an appropriate input file/folder.')
            with open(f'{save_dir}html_output.html', 'w') as f:
                f.write(f'Input file/folder {name} does not exist. Please use an appropriate input file/folder.')
                f.close()
            sys.exit(2)

    #Check to make sure the platemap is in the .csv format. 
    if root_names[1][-4:].lower() != ".csv":
        logging.error(f'Plate map is not in the .csv format. Please use a .csv file format.')
        with open(f'{save_dir}html_output.html', 'w') as f:
            f.write(f'Plate map is not in the .csv format. Please use a .csv file format.')
            f.close()
        sys.exit(2)

    if options.instrument == "Waters":
        if root_names[0][-4:].lower() != ".rpt":
            logging.error(f'LCMS data is not in the .rpt format. Please use a .rpt file format.')
            with open(f'{save_dir}html_output.html', 'w') as f:
                f.write(f'LCMS data is not in the .rpt format. Please use a .rpt file format.')
                f.close()
            sys.exit(2)

    if options.instrument == "Shimadzu":
        if not root_names[0].endswith("/"):
            root_names[0] = root_names[0] + "/"

        if len(glob.glob(root_names[0] + '*.daml')) == 0:
            logging.error('No .daml files found. Please specify alternative directory with .daml files present.')
            with open(f'{save_dir}html_output.html', 'w') as f:
                f.write(f'No .daml files found. Please specify a directory with .daml files present.')
                f.close()
            sys.exit(2)


    times["Initialise Script"] = time.perf_counter() - time_start 
    
    pre_rpt = time.perf_counter()
    #Depending on the instrument specified, fetch the LCMS data and convert it to the 
    #correct format. 
    if options.instrument == "Waters":
        rawData = getWatersData.rawWatersData(root_names[0])

        #Determine how the file records the well (e.g. "A,1", 1, "1,1", etc)
        #Note that this must be done before processing the data
        rawData.getWellFormat() 
        
        #Process the data types in turn
        rawData.processDAD()
        rawData.processMS()
        rawData.processUV()
        rawData.processTrace()
        rawData.processELSD()
        times["Import RPT File"] = time.perf_counter() - pre_rpt 
        logging.info(f'{rawData.rawDADTable["well"].nunique()} wells were found in the rpt.')

    #elif options.instrument == "Shimadzu":
    #    getShimadzuData.init(options)
    #    [dataTable, chroma, sample_IDs, total_area_abs] = getShimadzuData.getData(root_names[0])
    #    logging.info(f'{len(dataTable.items())} .daml data sources were found in the folder specified.')
    
    else:
        logging.error("This instrument is not currently supported by PyParse")
        sys.exit(2)

    #If the user did not provide the dimensions of the plate, set these values now. 
    if options.plate_col_no == 0 or options.plate_row_no == 0:
        options.plate_col_no = rawData.col_no
        options.plate_row_no = rawData.row_no

    #Select the rawData tables to use for the assignment
    #note that regardless of choice of detector, if data 
    #regarding the wavelengths each peak absorbs at is found, 
    #that data will be use to refine the selection when hits come to be validated. 

    msData = rawData.rawMSTable
    uvData = rawData.rawUVTable
    if options.detector == "UV":  
        lcData = rawData.rawDADTable
        traceData = rawData.rawTraceTable.loc[rawData.rawTraceTable["detector_type"] == "DAD"]
    else:
        lcData = rawData.rawELSDTable
        traceData = rawData.rawTraceTable.loc[rawData.rawTraceTable["detector_type"] == "ELSD"]

    #Import the structure data from the comma-separated values platemap file provided
    #by the user to initiate the compoundDF. 
    pre_import = time.perf_counter()            
    cpTable = cpAssignment.Assignment(root_names[1], options.plate_col_no, options.plate_row_no)
    cpTable.generateCPTable()
    cpTable.generateEMs("True") #Generate exact masses for each compound

    logging.info(f'{len(cpTable.cpTable.index)} compounds were imported.')

    times["Compound Import"] = time.perf_counter() - pre_import

    #Check that an internal standard was provided if the plot_type expects one
    internalSTD = list(cpTable.cpTable.loc[cpTable.cpTable["type"] == "InternalSTD", "name"])
    if options.plot_type in ["P/STD", "corrP/STD"] and len(internalSTD) == 0:
        logging.error("No internal standard was provided. Halting program.")
        with open(f'{save_dir}html_output.html', 'w') as f:
            f.write("No internal standard was provided for an analysis type which required one.")
            f.close()
        sys.exit(2)

    #Check to make sure the user didn't mistakenly assign the same compound to be both starting material
    #and product

    if cpTable.cpTable.index.has_duplicates:
        logging.error("The same compound was given two different compound types (e.g. both Product and Reactant). " 
                      "Check your platemap to ensure each compound appears in only one column.")
        with open(f'{save_dir}html_output.html', 'w') as f:
            f.write("The same compound was given two different compound types (e.g. both Product and Reactant). " 
                      "Check your platemap to ensure each compound appears in only one column.")
            f.close()
        sys.exit(2)

    #Sort the compounds in the table by their type. 
    all_products = list(cpTable.cpTable.loc[cpTable.cpTable["type"] == "Product", "name"])
    all_reactants = list(cpTable.cpTable.loc[cpTable.cpTable["type"] == "Reactant", "name"])
    all_stds = list(cpTable.cpTable.loc[cpTable.cpTable["type"] == "InternalSTD", "name"])
    all_byprods = list(cpTable.cpTable.loc[cpTable.cpTable["type"] == "Byproduct", "name"])
    all_cmps = all_products + all_reactants + all_stds + all_byprods

    cpTable.cpTable["name"] = pd.Categorical(cpTable.cpTable["name"], all_cmps)
    cpTable.cpTable.sort_values("name", inplace=True)

    #For each compound, find all hits using the dataTable, validate them, and append 
    #them to a list ready for insertion into the compoundDF pandas dataframe
    pre_hit_and_val = time.perf_counter()

    #Find hits for each compound
    cpTable.findHits(msData, lcData)

    #Validate the hits for each compound
    cpTable.validateHits(lcData, msData, uvData, options.time_abs_tol, options.massconf_threshold, 
                         options.uv_abs_tol, options.uv_cluster_threshold, options.uv_match_threshold, 
                         options.cluster_size_threshold, options.min_no_of_wells, options.validate, 
                         options.mass_or_area)
    
    times["Hit Finding and Validation"] = time.perf_counter() - pre_hit_and_val
    
    #Remove any duplicated assignments. If a peak has been assigned
    #to more than one compound, de-duplication occurs such that 
    #priority is given first to the internalSTD, then to limiting reactant,
    #then to product. 

    cpTable.removeDupAssigns(lcData, options.mass_or_area)
    logging.info(f'Duplicate assignments were removed.')

    #Find impurities
    if options.find_freq_imp == "True":
        pre_imp = time.perf_counter()
        cpTable.findImpurities(lcData, msData, save_dir, options.time_abs_tol, 
                                options.min_no_of_wells, options.mass_abs_tol)
        times["Finding Impurities"] = time.perf_counter() - pre_imp

    #Generate the output dataframe, used to facilitate plotting of graphs
    #Provide this with a list of the sample IDs for each well and 
    #the dimensions of the plate. 
    pre_output_table = time.perf_counter()
    output = genOutput.Output(rawData.sample_IDs, options.plate_col_no, options.plate_row_no)
    output.generateOutputTable(cpTable.cpTable, lcData)
    times["Generate Output Table"] = time.perf_counter() - pre_output_table

    #setBestWell relys on the output dataframe to run; setBestMS, setBestTime and setBestPurity 
    # are, by definition, the strongest m/z, retention time of the peak selected or peak area in 
    #that BestWell. Thus, these four  functions can only be run once the output dataframe has been run
    cpTable.setBestWell(output.df, lcData, options.plot_type)
    cpTable.setBestMS(lcData)
    cpTable.setBestTime(lcData)
    cpTable.setBestPurity(lcData)

    #Find any cases where two compounds, expected in the same well, 
    #were found to have similar retention times and therefore could have 
    #co-eluted, but in such a way that the mass ion wasn't observed for one of those 
    #compounds because the intensity of the other compound swamped it. 
    cpTable.findPotentialConflicts()

    #Plot the heatmaps
    pre_heatmap = time.perf_counter()
    plotting.plotHeatmaps(output.df, save_dir, options.plate_row_no, options.plate_col_no)
    times["Generate Heatmap"] = time.perf_counter() - pre_heatmap
    logging.info(f'A heatmap was generated using {options.plot_type} '
                f'as the index.')

    #Return a series of piecharts to the user, as long as it's not larger than
    #a 96 well plate. For larger plates, the piecharts become too small to be
    #useful. 
    if options.plate_row_no * options.plate_col_no < 97:
        #Get a list of byproduct names so that the plotPieCharts fn knows which columns 
        #to look up in the output table. 
        byproducts = cpTable.cpTable.loc[cpTable.cpTable["type"] == "Byproduct", "name"].values
        pre_pie = time.perf_counter()
        #Alter the piechart output depending on whether an internalSTD was specified
        #in the platemap or not. 
        if internalSTD != "":
            plotting.plotPieCharts(options.plot_type, output.df, save_dir, byproducts, 
                           options.plate_row_no, options.plate_col_no)
            logging.info(f'A set of pie-charts was generated using corrP/STD '
                        f'as the index.')
        else:
            plotting.plotPieCharts(options.plot_type, output.df, save_dir, byproducts, 
                           options.plate_row_no, options.plate_col_no)
            logging.info(f'A set of pie-charts was generated using Parea '
                        f'as the index.')
        times["Generate Piechart"] = time.perf_counter() - pre_pie


    #Plot a hit validation graph for each compound in the cpTable
    pre_val_graphs = time.perf_counter()
    cpTable.cpTable.apply(plotting.plotHitValidationGraph, args = (lcData, save_dir, options.plate_col_no,
                                                                   options.plate_row_no,), axis = 1)
    times["Generate Hit Validation Graphs"] = time.perf_counter() - pre_val_graphs

    #Plot an annotated chromatogram for each compound in the cpTable
    pre_chroma = time.perf_counter()
    cpTable.cpTable.apply(plotting.plotChroma, args=(cpTable.cpTable, lcData, msData, 
                                                     traceData, save_dir, options.plate_col_no,
                                                     options.filter_by_rt,), axis = 1)
    times["Generate Annotated Chromatograms"] = time.perf_counter() - pre_chroma

    #Plot a histogram and donut chart of Parea as a measure of plate success. 
    pre_donut = time.perf_counter()
    plotting.plotHistogram(output.df, save_dir)
    plotting.plotDonut(output.df, save_dir)
    times["Generate Histogram and Donut"] = time.perf_counter() - pre_donut
    logging.info(f'A histogram and donut chart were generated.')

    #Generate a set of PNG files to depict each compound
    cpTable.cpTable.apply(lambda row: generateMol(row["canonSMILES"], row["name"], save_dir), axis = 1)

    #Generate a location map for each compound
    cpTable.cpTable.apply(plotting.genLocationHeatmaps, args=(save_dir, options.plate_col_no, options.plate_row_no,), axis = 1)

    #Determine if there are any overlapping peaks for the best well of each compound. 
    cpTable.findOverlap(lcData)

    #Generate an csv of the output table.
    pre_csvs = time.perf_counter()
    if options.gen_csv == "True":
        csv = output.df.to_csv(f'{save_dir}outputTable.csv', index = False)
        newslice = cpTable.cpTable.loc[:, ["name", "canonSMILES", "mass1", "mass2", "mass3",
                                "time", "mass-", "mass+", "best_well", "best_purity", 
                                "overlaps", "conflicts"]]

        csv2 = newslice.to_csv(f'{save_dir}compoundtable.csv', index = False)
        logging.info('The CSV outputs were generated.')
    
    #Generate a ZIP file containing all analysis files and graphs

    if options.gen_zip == "True":

        with zipfile.ZipFile(f'{save_dir}output.zip', "w", allowZip64 = True) as myzip:
            for path, directories, files in os.walk(save_dir):
                for file in files:
                    if ".zip" not in file:
                        filename = os.path.join(path, file)
                        dst = f'{path.split("/")[-1]}/{file}'
                        myzip.write(filename, arcname = dst)

    times["Save .CSV Outputs"] = time.perf_counter() - pre_csvs

    #Log the time taken to complete the run. 
    total_time = time.perf_counter() - time_start
    times["Total time"] = total_time

    #build the HTML output file
    buildHTML(save_dir, cpTable.cpTable, options.analysis_name, times = times)

    
    
    #Generates and saves an analytical table containing all information concerning inputs and
    #outputs of the plate. 
    #if options.gen_atable == "True":
    #    genAnalyticalTable(platemapDF, compoundDF, save_dir, sample_IDs, products, SMs, by_products, internalSTD)

    
    logging.info(f'The analysis was completed in {total_time} seconds.')
    print(f'The analysis was completed in {round(total_time, 2)} seconds.')



if __name__ == "__main__":
    try:
        print("Running analysis....")
        main()
    except Exception as e:
        print(e)
        logging.exception("A fatal exception occurred. Contact administrator.")