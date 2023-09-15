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
from jinja2 import Environment, FileSystemLoader




def importStructures(filename):
    """
    From a given CSV, with defined header names, and one row per well, 
    deconvolute the information to give a dataframe of one compound 
    per row, indexed by canonical SMILES. The well must be specified 
    using a capital letter from A-Z to describe the row, and a positive integer 
    to describve the column. Column numbers should be written as 1, 2, etc 
    as opposed to 01, 02, etc.  

    :param filename: file name and directory as a string
    :return: List comprising [Pandas dataframe, name of internal STD, list of names of the starting materials, 
                list of names of the products, list of names of the byproducts]
    """
    
    #read csv file into dataframe
    #replace empty cells with an empty string
    #convert all column names to lower case and remove whitespace
    inputCSV = pd.read_csv(filename)
    inputCSV.fillna("", inplace=True)
    inputCSV.columns = inputCSV.columns.str.strip().str.lower()

    compounds = {}
    internalSTD = ""
    SMs = []
    products = []
    by_products = []
    global product_count, SM_count
    product_count = 0
    SM_count = 0

    #Fn to convert a well name like B5 to machine format (11)
    def convertWellToNum(well):
        row = well[0]
        column = well[1:]
        #if the format of the row/column conforms to expectations
        if isinstance(int(column[0]), int):
            return int((ord(row) - 65) * options.plate_col_no + int(column))
        #Plates with more than 26 rows are unsupported at present. 
        else:
            logging.info("The well specified implies an unsupported plate.")
            with open(f'{save_dir}html_output.html', 'w') as f:
                f.write("The well specified implies an unsupported plate.")
                f.close()
            sys.exit(2)

    #fn to import process each compound into the right format
    def addCompound(column, row):
        if column in row:
            if row[column] != "":
                g_smiles = row[column]
                
                #generate a canonical smiles for use as an index
                mol = Chem.MolFromSmiles(g_smiles.strip())
                index = Chem.MolToSmiles(mol)
                
                #Convert well ID to an integer
                well = convertWellToNum(row["well"])
                
                if index in compounds:
                    compounds[index]["locations"].append(well)
                else:
                    compounds[index] = {"locations": [well]}
                    compounds[index]["g_smiles"] = g_smiles
                    compounds[index]["wellname"] = row["well"]
                    
                    #store the common name given by the user
                    if column == "desired product smiles":
                        name = ""
                        if "desired product name" in row and row["desired product name"] != "":
                            name = row["desired product name"]
                        else:
                            global product_count
                            product_count = product_count + 1
                            name = f'Product{product_count}'
                        compounds[index]["name"] = name

                        if "desired product rt" in row and row["desired product rt"] != "":
                            compounds[index]["rt"] = float(row["desired product rt"])
                        else:
                            compounds[index]["rt"] = 0
                        
                        compounds[index]["type"] = "Product"
                        return name

                    elif column == "internalstd smiles":
                        name = ""
                        if "internalstd name" in row and row["internalstd name"] != "":
                            name = row["internalstd name"]
                        else:
                            name = 'InternalSTD'

                        compounds[index]["name"] = name

                        if "internalstd rt" in row and row["internalstd rt"] != "":
                            compounds[index]["rt"] = float(row["internalstd rt"])
                        else:
                            compounds[index]["rt"] = 0
                        
                        compounds[index]["type"] = "InternalSTD"
                        return name
                        
                    elif column == "limiting reactant smiles":
                        name = ""
                        if "limiting reactant name" in row and row["limiting reactant name"] != "":
                            name = row["limiting reactant name"]
                        else:
                            global SM_count
                            SM_count = SM_count + 1
                            name = f'Reactant{SM_count}'

                        compounds[index]["name"] = name

                        if "limiting reactant rt" in row and row["limiting reactant rt"] != "":
                            compounds[index]["rt"] = float(row["limiting reactant rt"])
                        else:
                            compounds[index]["rt"] = 0
                        
                        compounds[index]["type"] = "Limiting Reagent"
                        return name

                    #use regex to find the byproduct columns (as multiple
                    #byproducts are permitted and each requires a separate column
                    #in the platemap.) 
                    elif re.search("^byproduct.*smiles$", column):
                        name = ""
                        byprod_num = column.split("byproduct")[1].split(" ")[0].strip()
                        if f'byproduct{byprod_num} name' in row and row[f'byproduct{byprod_num} name'] != "":
                            name = row[f'byproduct{byprod_num} name']
                        else:
                            name = f'Byproduct{byprod_num}'

                        compounds[index]["name"] = name

                        if f'byproduct{byprod_num} rt' in row and row[f'byproduct{byprod_num} rt'] != "":
                            compounds[index]["rt"] = float(row[f'byproduct{byprod_num} rt'])
                        else:
                            compounds[index]["rt"] = 0

                        compounds[index]["type"] = "Byproduct"
                        return name
        
        #Returns zero if the column isn't found in that row,
        #or if the compound has already been added from a previous well. 
        return 0

                
    #For each row (well), look to specific columns for the product, limiting reactant, etc
    #Store the wells that were in the input csv for future use
  
    columns = inputCSV.columns         
    for index, row in inputCSV.iterrows():
        #Store the returned canonical SMILES for future use
        new_prod = addCompound("desired product smiles", row)
        if new_prod != 0:
            products.append(new_prod)

        byproduct_cols = [i for i in columns if re.search("^byproduct.*smiles$", i)]
        for j in byproduct_cols:
            new_byprod = addCompound(j, row)
            if new_byprod != 0:
                by_products.append(new_byprod)
        
        new_reactant = addCompound("limiting reactant smiles", row)
        if new_reactant != 0:
            SMs.append(new_reactant)

        new_standard = addCompound("internalstd smiles", row)
        if new_standard != 0:
            internalSTD = new_standard
 
    #Convert dictionary to Pandas dataframe
    compoundDF = pd.DataFrame.from_dict(compounds, orient="index")
    
    #Location to store the isotopic mass of the compounds
    mass1 = []
    mass2 = []
    mass3 = []
    if options.calc_boc == "True":

        #Use RDkit to find the isotopic mass of each compound
        #as well as likely other masses, e.g. boc deprotection and
        #halide isotopes. 
        for index, row in compoundDF.iterrows():
            mol = Chem.MolFromSmiles(index)
            mw1 = round(Descriptors.ExactMolWt(mol), 2)
            #Add the parent mass to the list
            mass1.append(mw1)
            #Set default values for mw2 and mw3
            mw2 = 0
            mw3 = 0
            try:
                #Use rdkit to perform boc degradation
                rxn1 = AllChem.ReactionFromSmarts("[NX3,n:1][C:2](=[O:3])[O:4][C]([CH3])([CH3])[CH3]>>[*:1][C:2](=[O:3])[O:4]")
                new_mol1 = rxn1.RunReactants((mol, ))[0][0]
                #Sanitise the molecule to make sure that a sensible molecule was produced. 
                Chem.SanitizeMol(new_mol1)
                mw2 = round(Descriptors.ExactMolWt(new_mol1), 2)
                
                #Use rdkit again to get to fully deprotected molecule
                rxn2 = AllChem.ReactionFromSmarts("[NX3,n:1][C](=[O])[O][C]([CH3])([CH3])[CH3]>>[*:1][H]")
                new_mol2 = rxn2.RunReactants((mol, ))[0][0]
                Chem.SanitizeMol(new_mol2)
                mw3 = round(Descriptors.ExactMolWt(new_mol2), 2)
                
                logging.info(f'A Boc group was found and disconnected for {row["name"]}.')
                
            except Exception as e:
                
                #If the boc degradation failed, either because a boc group was
                #not found on this molecule or for other reasons, check for 
                #halide isotope patterns. 
                if "Cl" in index or "Br" in index:
                    mw2 = round(mw1 + 2, 2)
        
            #Append new masses to list 
            mass2.append(mw2)
            mass3.append(mw3)
        
    else:
        for index, row in compoundDF.iterrows():
            mol = Chem.MolFromSmiles(index)
            mw1 = round(Descriptors.ExactMolWt(mol), 2)
            mass1.append(mw1)
            
            if "Cl" in index or "Br" in index:
                mass2.append(round(mw1 + 2, 2))
                mass3.append(0)
            else:
                mass2.append(0)
                mass3.append(0)

    #Append the mass data to the dataframe
    compoundDF["mass1"] = mass1
    compoundDF["mass2"] = mass2
    compoundDF["mass3"] = mass3

    return [compoundDF, internalSTD, SMs, products, by_products]
    

def getUserReadableWell(wellno):
    """
    Converts the well as a number into a user-friendly string,
    e.g. well 11 becomes "B5" for a 4*6 well plate

    :param wellno: An integer representing a specific well on the plate
    
    :return: A string representing a specific well on the plate
    """
    
    rowVal = math.floor((wellno-1) / options.plate_col_no)
    colVal = (wellno) % options.plate_col_no
    if colVal == 0:
        colVal = options.plate_col_no
    
    label = f'{chr(ord("@")+(rowVal)+1)}{colVal}'
    return label
    
def getChromatogram(spectrum):
    """
    Takes the region of text in rpt file corresponding to the chromatogram,
    and extracts the data into two lists, corresponding to x-values and y-values. 
    
    :param spectrum: Section of rpt file as string
    
    :return: A list of 2 lists, [x-values, y-values]
    """            
    
    x_val = []
    y_val = []
    
    lineData = spectrum.split("\n")
    length = len(lineData)
    n_val = math.ceil(length / options.points_per_trace)

    for index, value in enumerate(lineData[1:]):
        if value == "}": #stop for loop if end of data section is reached
            break
        elif index % n_val == 0:
            
            if value != "{":
                data = value.split("\t")
                x_val.append(float(data[0]))
                y_val.append(float(data[1]))
    
    return [x_val, y_val]
        
        
        
def getMSData(spectrum):
    """
    Takes the specific region of the rpt file pertaining to 
    m/z data for a specific peak in a specific well, and 
    returns a 2-D list containing all m/z peaks and their 
    normalised intensity.
    
    :param spectrum: Section of rpt file as string
    
    :return: 2-D list in following format: 
        [m/z value, normalised intensity of that value]
    """
    
    masses = []
    total = 0
    
    lineData = spectrum.split(";Mass\t% BPI")[1].split("\n")
    for line in lineData[1:]:
        if line == "}": #stop the for loop if end of MSData section is reached
            break
    
        massData = line.split("\t")
        if len(massData) == 2:
            floatData = [float(i) for i in massData] #convert all data to float
            masses.append(floatData)
            
            total = total + floatData[1]
    refined_masses = []
    for i in masses:
        if math.floor((i[1]/total)*100) > 0:
            refined_masses.append([i[0], i[1]])
    
    return refined_masses
   
def getUVData(spectrum):
    """
    Takes the specific region of the rpt file pertaining to 
    UV absorbance spectrum data for a specific peak in a specific
    well, and returns a list containing all the maxima of that spectrum.
    The height of the maxima must be greater than the min_uv_threshold
    specified in options.
    
    :param spectrum: Section of rpt file as string
    
    :return: UV maxima as a list
    """

    UVmaxima = []
    
    lineData = spectrum.split(";Mass\t% BPI")[1].split("\n")
    UVx = []
    UVy = []
    for line in lineData:
        
        if line == "}":#stop for loop if end of UVData section is reached
            break
        
        UVdata = line.split("\t")
        if len(UVdata) == 2:
            UVx.append(float(UVdata[0]))
            UVy.append(abs(float(UVdata[1])))
    if UVy[0] > UVy[1] and UVy[0] > options.min_uv_threshold:
        UVmaxima.append(UVx[0])
    for i in range(1, len(UVy)-1):
        if UVy[i] > UVy[i-1] and UVy[i] > UVy[i+1] and UVy[i] > options.min_uv_threshold:
            UVmaxima.append(UVx[i])
    if UVy[-1] > UVy[-2] and UVy[-1] > options.min_uv_threshold:
        UVmaxima.append(UVx[-1])

        
    return UVmaxima            

def convertRPT2Dict(filename):
    """
    Converts the input .rpt file into a dictionary where each index corresponds to a well.
    Each dictionary value contains a list of peaks, where each peak is represented as a dictionary
    with the peak area, m/z and UV (if available) data present. 
    
    :param filename: Address of the input file
    
    :return: List comprising: [a dictionary of all peaks in all wells, a list containing each chromatogram,
                a list of the sample IDs for each well] 
    """
    
    wellData = []
    masterTable = {}
    chroma = {}
    sample_IDs = {}

    with open(filename, errors = "ignore") as f:
        fullText = f.read()
        text_split = fullText.split("[SAMPLE]") #Split the file into individual wells
        wellData = text_split[1:] 
        

    def findDataEntry(peaks, peakID):
        """
        The program should match up equivalent peaks from different 
        detectors. Detectors may have subtly different retention times, but 
        the rpt file tags equivalent peaks with the same peakID. Use this value 
        to match up peaks from different detectors. 

        :param peaks: a dictionary of peaks found so far
        :param peakID: a string; the peakID of a specific peak for a specific detector
        :return: a float corresponding to an index for the peaks dictionary, or false if none was found. 
        """
        
        matches = [index for index, i in peaks.items() if i["peakID"] == peakID]
        if len(matches) == 1:
            return matches[0]
        else:
            return False
    
    for i in range(len(wellData)):
        well = wellData[i]
        wellno = -1
        functions = well.split("[FUNCTION]") #Split data by function (Prelude, MS+, MS- or UV)
        
        #Get the well no from the specific point in the prelude section
        #There are numerous ways for an .rpt file to describe the well
        #Each method is tried in turn to account for these differences
        try:
            #If the well is simply an integer between 1 and infinity
            #Single line function to trim full string to just the well number used
            wellno = int(functions[0].split("Well")[1].split("\n")[0].split(":")[1].strip()) 
        except:
            column = functions[0].split("Well")[1].split("\n")[0].split(":")[1].split(",")[0].strip()
            row = functions[0].split("Well")[1].split("\n")[0].split(":")[1].split(",")[1].strip()
            try:
                #If the well is a combination of two integers designating column and row
                wellno = (int(row)-1) * options.plate_col_no + int(column) 
            except:
                #If the well is a combination of an integer for the column and a letter 
                #for the row
                
                try:
                    row_to_int = ord(row) - 64
                    wellno = (row_to_int-1) * options.plate_col_no + int(column)
                except:
                    #If that fails, try a combination of letter for the column and 
                    #a letter for the row
                    try:
                        row_to_int = ord(row) - 64
                        col_to_int = ord(column) - 64
                        wellno = (row_to_int-1) * options.plate_col_no + col_to_int
                    except:
                        #If that fails, try a combination of letter for the column and 
                        #an integer for the row. 
                        try:
                            col_to_int = ord(column) - 64
                            wellno = (int(row)-1) * options.plate_col_no + col_to_int
                        except:

                            logging.info("Unable to get well number. Terminating program.")
                            sys.exit(2)

                
        #get the sample id, and load it into a list
        #This list will be added as a column to the outputTable
        well_ID = functions[0].split("SampleID")[1].split("\n")[0].strip()  
        sample_IDs[wellno] = well_ID
         
        peaks = {}
        for function in functions[1:]: #Ignore the prelude data
            lines = function.split("\n")
            spectra = function.split("[SPECTRUM]")[1:] #split by spectrum (i.e. each peak)
            
            #get peakarea for this peak
            if lines[4] == "Type\tDAD ":
                chromatograms = function.split("[CHROMATOGRAM]")[1:]
                for chromatogram in chromatograms:
                    c_lines = chromatogram.split("\n")
                    if c_lines[3] == "Description\tDAD: TIC":
                        chroma[wellno] = getChromatogram(chromatogram.split("[TRACE]")[1]) #get chromatogram for this well
                        wellpeaks = chromatogram.split("[PEAK]")[1:]
                        for peak in wellpeaks:
                            
                            retTime = float(peak.split("Time")[1].split("\n")[0].strip())
                            peakWidth = peak.split("Peak\t")[1].split("\n")[0].split("\t")
                            peakArea = float(peak.split("Area %Total")[1].split("\n")[0].strip())
                            peakAreaAbs = float(peak.split("AreaAbs")[1].split("\n")[0].strip())
                            peakID = int(peak.split("Peak ID")[1].split("\n")[0].strip())

                            peak_index = False
                            if len(list(peaks)) > 0:
                                peak_index = findDataEntry(peaks, peakID)
                            
                            if peak_index == False:
                                peaks[retTime] = {
                                        "MS+": [], 
                                        "MS-": [], 
                                        "area": 0,
                                        "areaAbs": 0, 
                                        "UV": [],
                                        "time": retTime,
                                        "pStart": 0,
                                        "pEnd": 0, 
                                        "peakID": peakID
                                        }
                            elif peak_index != retTime:
                                #If data has already been entered for this peak by a different
                                #function (i.e. MS+, MS-) using a slightly different retention time,
                                #copy this data over to a new entry
                                #which uses the retention time observed in the chromatogram 
                                #i.e. we don't want some retention times from the MS+, and others
                                #from the chromatogram - we want a standard output
                                peaks[retTime] = {
                                        "MS+": peaks[peak_index]["MS+"], 
                                        "MS-": peaks[peak_index]["MS-"], 
                                        "area": peaks[peak_index]["area"], 
                                        "areaAbs": peaks[peak_index]["areaAbs"],
                                        "UV": peaks[peak_index]["UV"],
                                        "time": retTime,
                                        "pStart": peaks[peak_index]["pStart"],
                                        "pEnd": peaks[peak_index]["pEnd"], 
                                        "peakID": peaks[peak_index]["peakID"],

                                        }
                                #Delete the old entry now the data has been copied over
                                #to use the correct retention time. 
                                del peaks[peak_index]
                                
                            peaks[retTime]["area"] = peakArea
                            peaks[retTime]["areaAbs"] = peakAreaAbs
                            peaks[retTime]["pStart"] = float(peakWidth[0])
                            peaks[retTime]["pEnd"] = float(peakWidth[1])
                    
            for spectrum in spectra:
                peakID = int(spectrum.split("Peak ID")[1].split("\n")[0].strip())
                retTime = float(spectrum.split("Time")[1].split("\n")[0].strip())
                peak_index = False
                #Find matching peaks in other functions by reference to the peakID
                if len(list(peaks)) > 0:
                    peak_index = findDataEntry(peaks, peakID)
                if peak_index == False:
                    peaks[retTime] = {
                            "MS+": [], 
                            "MS-": [], 
                            "area": 0, 
                            "areaAbs": 0,
                            "UV": [],
                            "time": retTime,
                            "pStart": 0,
                            "pEnd": 0,
                            "peakID": peakID
                            }
                else:
                    retTime = peak_index
                    
                #For each function, parse all the associated data and store
                #as a float. 
                #For UV, take only the lambaMax
                if lines[3] == "IonMode\tES+ ":
                    data = getMSData(spectrum)
                    peaks[retTime]["MS+"] = data
                      
                elif lines[3] == "IonMode\tES- ":
                    data = getMSData(spectrum)
                    peaks[retTime]["MS-"] = data
                
                elif lines[4] == "Type\tDAD ":
                    data = getUVData(spectrum)
                    peaks[retTime]["UV"] = data
            
        logging.debug(f'{len(peaks)} peaks found in well {wellno}.')             
    
        masterTable[wellno] = peaks.values()

    return [masterTable, chroma, sample_IDs]
    
    
def findHits(compound, dataTable):
    """
    For a given compound, look at each peak in each well 
    to find a suitable match based on m/z data.
    
    :param compound: A Series corresponding to a specific compound
                     from the compoundDF pandas dataframe
    :param dataTable: A dictionary of all peaks in all wells for the plate
    
    :return: a list of hits, where each item in the list is a dictionary
    """

    hits = []
    for index, well in dataTable.items():
        
        if index in compound["locations"]:
            
            for peak in dataTable[index]:
                
                if peak["area"] > options.min_peak_area:
                    
                    mass_plus = []
                    mass_minus = []
                    new_hit = {
                            "time": peak["time"],
                            "area": peak["area"],
                            "areaAbs": peak["areaAbs"],
                            "UV": peak["UV"], 
                            "well": index,
                            "pStart": peak["pStart"],
                            "pEnd": peak["pEnd"],
                            "mass_conf": 0,
                            "mass+": 0,
                            "mass-": 0
                            }
                    total_ms_plus = sum([i[1] for i in peak["MS+"]])
                   
                    for [mass, intensity] in peak["MS+"]:
                        intensity_pc = intensity/total_ms_plus * 100
                        for targetmass in [compound["mass1"], compound["mass2"], compound["mass3"]]:
                            intensity_pc = intensity/total_ms_plus * 100
                            if math.isclose(targetmass+1.01, mass, abs_tol=options.mass_abs_tol):
                                mass_plus.append([mass, intensity_pc])
                            elif math.isclose((targetmass+2.02)/2, mass, abs_tol=options.mass_abs_tol) and options.calc_higherions == "True":
                                mass_plus.append([mass, intensity_pc])
                            elif math.isclose((targetmass+3.03)/3, mass, abs_tol=options.mass_abs_tol) and options.calc_higherions == "True":
                                mass_plus.append([mass, intensity_pc])

                                
                    total_ms_minus = sum([i[1] for i in peak["MS-"]])        
                    for [mass, intensity] in peak["MS-"]:
                        for targetmass in [compound["mass1"], compound["mass2"], compound["mass3"]]:
                            intensity_pc = intensity/total_ms_minus * 100
                            if math.isclose(targetmass-1.01, mass, abs_tol=options.mass_abs_tol):
                                mass_minus.append([mass, intensity_pc])
                    
                    #If the total mass confidence for this peak/compound pair 
                    #is greater than min_massconf_threshold, store this peak
                    #The peak is also stored if the retention time of the hit is close
                    #to that specified in the platemap (if one was specified) - this 
                    #allows compounds that don't ionise to still be analysed
                    
                    total_mc = sum([i[1] for i in mass_plus]) + sum([i[1] for i in mass_minus])

                    if total_mc > options.min_massconf_threshold or math.isclose(new_hit["time"], compound["rt"], abs_tol=options.time_abs_tol):
                        new_hit["mass_conf"] = total_mc
                        if len(mass_plus) > 0:
                            new_hit["mass+"] = max(mass_plus, key = lambda x: x[1])[0]
                        else: 
                            new_hit["mass+"] = "-"
                        if len(mass_minus) > 0:
                            new_hit["mass-"] = max(mass_minus, key = lambda x: x[1])[0]
                        else:
                            new_hit["mass-"] = "-"

                        hits.append(new_hit)
            
    return hits

def refineClusterByTime(cluster, comments, expected_rt):
    """
    Takes in input cluster of all the hit peaks, 
    and refines them by finding a mid-value for the retention
    time based on which hit has the greatest number of nearest neighbours. 
    Sorts the best hits into "green", uncertain ones into "orange" 
    and those where another peak closer to the mid-value was found
    in the same well into "discarded". 
    
    :param cluster: list of dictionaries, where each dictionary is a hit
    :param comments: A list of comments for that structure so far.
    
    :return: List comprising [a dictionary for the refined cluster, list of comments]
    """

    refined_cluster = {
            "green":[],
            "orange": [],
            "discarded": [],
            }
    
    mid_values = []
    mid_value = 0

    if expected_rt != 0:
        mid_value = expected_rt
    else: 
        for i in cluster:
            mid_value = i["time"]
            mid_values.append([mid_value, len([j["time"] for j in cluster if math.isclose(j["time"], mid_value, abs_tol = options.time_abs_tol/4)])])
        mid_value = max(mid_values, key = lambda x: x[1])[0]
    
    #sort the peaks by the well they occupy
    peaks_by_wells = {}
    for peak in cluster:
        if peak["well"] not in peaks_by_wells:
            peaks_by_wells[peak["well"]] = []
        peaks_by_wells[peak["well"]].append(peak)
    
    #For each well, select the peak that's closest to the mid-value in cases
    #where there was more than one hit in that cluster in one well
    for index, well in peaks_by_wells.items():
        if len(well) > 1:
            min_diff = min([abs(x["time"]-mid_value) for x in well])
            
            for peak in well:
                
                if abs(peak["time"]-mid_value) == min_diff:
                    refined_cluster["green"].append(peak)
                else:
                    refined_cluster["discarded"].append(peak)
                    comments.append(f'Peak at {peak["time"]} '
                                f'in well {getUserReadableWell(peak["well"])} was discarded '
                                'as there was an alternative peak '
                                'in the same well which was closer to the '
                                'mid-point of the cluster.')
        else:
            refined_cluster["green"].append(well[0])

    #Refine these hits further by finding those which are within time_abs_tol/2
    #of the mid-value. Any others are marked as tentative and the user is alerted.         
    ref2_cluster = {
            "green":[],
            "orange": refined_cluster["orange"],
            "discarded": refined_cluster["discarded"],
            }        
            
    for peak in refined_cluster["green"]:
        if math.isclose(peak["time"], mid_value, abs_tol = options.time_abs_tol / 2):
            ref2_cluster["green"].append(peak)
        else:
            ref2_cluster["orange"].append(peak)
            comments.append(f'Peak at {peak["time"]} in '
                           f'well {getUserReadableWell(peak["well"])} was '
                           'marked as tentative as it was found to be too '
                           'far from the mid-value of the cluster.')
    
    
    return [ref2_cluster, comments]
                        
    
    
def refineClusterByUV(cluster, UVdatafound, comments):
    """
    Takes in input cluster of all the hit peaks, 
    and refines the cluster by ensuring all peaks have a similar set
    of UV maxima. Those which do are left in "green", those which
    don't are moved to the "orange" category
    
    :param cluster: a dict, with list of dicts for each header
    :param UVdatafound: boolean for whether the rpt data contains UV data
    :param comments: A list of comments for that structure so far
    
    :return: List comprising [a dictionary for the refined cluster, list of comments]
    """
    
    refined_cluster = {
                        "green": [],
                        "orange": cluster["orange"],
                        "discarded": cluster["discarded"],
                        }
    
    if UVdatafound == True and len(cluster["green"]) > 0:
        UVclusters = []
        
        for peak in cluster["green"]:
            for UV in peak["UV"]:
                if len(UVclusters) == 0:
                    UVclusters.append([UV])
                    
                else:
                    clusterFound = False
                    for UVcluster in UVclusters:
                        if math.isclose(UVcluster[-1], UV, abs_tol=options.uv_abs_tol):
                            UVcluster.append(UV)
                            clusterFound = True
                            break
                    if not clusterFound:
                        UVclusters.append([UV])
           
        meanUV = []
        #Check to ensure UVclusters isn't an empty set
        if len(UVclusters) > 0:
            lengthOfMostCommon = max([len(UVcluster) for UVcluster in UVclusters])

            for UVcluster in UVclusters:
                if len(UVcluster) >= lengthOfMostCommon * options.uv_cluster_threshold:
                    meanUV.append(mean(UVcluster))
            
            for peak in cluster["green"]:
                
                intersection = []
                for meanvalue in meanUV:
                    for UV in peak["UV"]:
                        if math.isclose(meanvalue, UV, abs_tol = options.uv_abs_tol):
                            intersection.append(meanvalue)
                            
                if len(intersection) >= len(meanUV) * options.uv_match_threshold:
                    refined_cluster["green"].append(peak)
                else:
                    refined_cluster["orange"].append(peak)
                    comments.append(f'Peak at {peak["time"]} in well {getUserReadableWell(peak["well"])} '
                    f'was marked tentative due to mismatch in UV maxima with the rest of the cluster.')
        #If UVclusters is an empty set, pass through the original hits without further
        #validation. 
        else:
            comments.append('UV validation not permitted where cluster contains peaks with no UV data. '
                'UV validation was not performed.')
            refined_cluster["green"] = cluster["green"]
    else:
        comments.append('UV data was not found for the plate, so UV validation was not performed.')
        refined_cluster["green"] = cluster["green"]        
            
    return [refined_cluster, comments]
    
def refineClusterByMassConf(cluster, comments):
    """
    Takes in input cluster of all the hit peaks, 
    and refines them by ensuring all peaks have a similar mass confidence
    to the cluster's mean. Those which do are left in "green"; 
    those which don't are moved to the "orange" category.
    
    :param cluster: a dict, with list of dicts for each header
    :param comments: A list of comments for the compound so far
    
    :return: List comprising [a dictionary for the refined cluster, list of comments]
    """
    
    refined_cluster = {
                        "green": [],
                        "orange": cluster["orange"],
                        "discarded": cluster["discarded"],
                        }
    #find what the mean mass confidence is of all peaks
    #currently under the "green" category
    total = 0
    mean = 0
    if len(cluster["green"]) > 0:
        for peak in cluster["green"]:
            total += peak["mass_conf"]
        mean = total / len(cluster["green"])
    else:
        for peak in cluster["orange"]:
            total += peak["mass_conf"]
        mean = total / len(cluster["orange"])
    
    
    
    for peak in cluster["green"]:
        if peak["mass_conf"] < mean * options.massconf_threshold:
            refined_cluster["orange"].append(peak)
            comments.append(f'Peak at {peak["time"]} in well {getUserReadableWell(peak["well"])} '
                  f'was marked tentative due to poor mass confidence in comparison to the rest of the cluster.')
        else:
            refined_cluster["green"].append(peak)
    
    return [refined_cluster, comments]

def selectClusterByMassConf(clusters):
    """
    If more than one cluster was found for the compound, 
    this function is called to try to select a single cluster based on 
    which cluster has the highest mean massConf. If more than one cluster
    has a close-to-highest-mean massconf, take them all. 
    
    :param clusters: a list of dictionaries, with list of dictionaries for each header
    :return refined_clusters: a list of dictionaries, with list of dictionaries for each header
    
    :return discarded_clusters: a list of dictionaries, with list of dictionaries for each header
    """
    
    refined_clusters = []
    discarded_clusters = []
    
    means = []
    #find the mean massconf for each cluster
    for cluster in clusters:
        total = 0
        for peak in cluster["green"]:
            total += peak["mass_conf"]
        if total != 0:
            mean = total / len(cluster["green"])
            means.append(mean)
        else:
            means.append(0)
    #find the maximum mean value to compare all clusters against
    max_mean = max(means)
    
    for i in range(len(means)):
        if max_mean * options.massconf_threshold < means[i] or (max_mean == 0 and means[i] == 0):
            refined_clusters.append(clusters[i])
        else:
            discarded_clusters.append(clusters[i])
            
    return [refined_clusters, discarded_clusters]
            
def selectClusterBySize(clusters):
    """
    If more than one cluster was found for the compound, 
    this function is called to try to select a single cluster based on 
    which cluster is the largest. If more than one cluster
    has a close-to-largest size, take them all. 
    
    :param clusters: a list of dictionaries, with list of dictionaries for each header
    :return refined_clusters: a list of dictionaries, with list of dictionaries for each header
    
    :return discarded_clusters: a list of dictionaries, with list of dictionaries for each header
    """
    
    refined_clusters = []
    discarded_clusters = []
    
    lengths = [len(cluster["green"]) + len(cluster["orange"]) for cluster in clusters]
    max_length = max(lengths)
    
    for cluster in clusters:
        if len(cluster["green"]) + len(cluster["orange"]) > max_length * options.cluster_size_threshold:
            refined_clusters.append(cluster)
        else:
            discarded_clusters.append(cluster)
            
    #If there is still more than one cluster, repeat the process using only those peaks
    #that haven't been marked as suspicious
    if len(refined_clusters) > 1:
        
        lengths = [len(cluster["green"]) for cluster in clusters]
        max_length = max(lengths)
        
        #As long as at least one cluster had a green hit, filter the clusters by size
        if max_length != 0:
            refined_clusters = []
            
            for cluster in clusters:
                if len(cluster["green"]) > max_length * options.cluster_size_threshold:
                    refined_clusters.append(cluster)
                else:
                    discarded_clusters.append(cluster)
            
    return [refined_clusters, discarded_clusters]    

def validateHits(cpname, peakList, expected_rt):
    """
    Looks at all hits for a given compound, and refines
    the list based on retention time, massConf and UV data. 
    The goal is to end up with only high confidence hits, and no more than
    one hit per well. 
    
    :param cpname: string of the compound name
    :param peakList: A list of all hits for (peaks assigned to) that compound
    :param expected_rt: An integer for the expected retention time of the compound
    
    :return: A dictionary, where each header contains a list of dictionaries (peaks), 
            along with other useful variables to aid plot generation. 
    """
    #Initialise variables which end up in the final output
    output = {
            "green": [],
            "discarded": [],
            "discarded_by_cluster": []
        }
    comment_text = []
    mean_time = 0
    
    
    #Store the mean retention times of each 
    #cluster to enable annotation in the output plots. 
    cluster_bands = []
    
    #check to see if UV data is available to validate against
    UVdatafound = True if len(peakList[0]["UV"]) > 0 else False

    #find how many wells have at least one hit:
    hitwells = set([i["well"] for i in peakList])
    
    comment_text.append(f'{len(peakList)} hits were found for {cpname}.')
    
    clusters = []
        
    #lambda: inline function definition
    #sort the peaks in order of ascending retention time
    
    peakList.sort(key = lambda x: x["time"])
    
    #cluster the peaks by retention time into groups to better
    #classify which hits are likely to be genuine and which
    #are false positives. 
    for i in range(len(peakList)):
        
        if len(clusters) == 0:
            clusters.append([peakList[i]])
        else:
            clusterFound = False
            for cluster in clusters:
                mean_rt = mean([i["time"] for i in cluster])
                
                if math.isclose(mean_rt, peakList[i]["time"], abs_tol=options.time_abs_tol):
                    cluster.append(peakList[i])
                    clusterFound = True
                    break
            if not clusterFound:
                clusters.append([peakList[i]])
                
    comment_text.append(f'{len(clusters)} clusters were found for {cpname}.')
    
    #find the average retention time for each cluster, and store in cluster_bands
    #use this to label the hit validation graph and to refine by expected_rt.  
    
    for i in range(len(clusters)):
        avg_time = sum([j["time"] for j in clusters[i]]) / len(clusters[i])
        cluster_bands.append(round(avg_time,2))
        
        for peak in clusters[i]:
            peak["cluster"] = i
    
    #If the user has specified a retention time, we should select only the cluster 
    #that is closest to that retention time, and within options.time_abs_tol
    if expected_rt != 0:
        suitable_clusters = [index for index, i in enumerate(cluster_bands) if math.isclose(i, expected_rt,
                                                                abs_tol=options.time_abs_tol)]
        
        #If there is more than one cluster close to the specified retention time
        #take only the cluster which is closest 
        if len(suitable_clusters) > 1:
            diffs = [cluster_bands[index]-expected_rt for index in suitable_clusters]
            index_min = min(range(len(diffs)), key=diffs.__getitem__)
            clusters = [clusters[suitable_clusters[index_min]]]
            #Update the cluster bands to only include the correct label
            cluster_bands = [cluster_bands[suitable_clusters[index_min]]]

            comment_text.append("<strong>Multiple clusters were found close the specified"
                                " retention time.</strong>")
            comment_text.append(f'<strong>Cluster {index_min} was selected as it was closest'
                                ' to the specified retention time.</strong>')
        elif len(suitable_clusters) == 1:
            clusters = [clusters[suitable_clusters[0]]]
            cluster_bands = [cluster_bands[suitable_clusters[0]]]
            comment_text.append("<strong>A single cluster was found close the specified"
                                " retention time and this was selected for analysis.</strong>")
        else:
            comment_text.append("<strong>No cluster was found near to the specified "
                                "retention time. Proceeding with analysis using all "
                                f'{len(clusters)} clusters.</strong>')
    
    
    #Note: if options.validate is set to true, validation may still run 
    #using the one cluster identified from the given retention time
    #(assuming one was)      
    if len(hitwells) > options.min_no_of_wells and options.validate == "True":
           
        #iterate through the clusters and refine each cluster
        #by time, then by UV, then by massconf.
        refined_clusters = []
        for i in range(len(clusters)):
            
            #Refine using the retention time of the hits
            ref_cluster = refineClusterByTime(clusters[i], comment_text, expected_rt)
            comment_text = ref_cluster[1]   
            comment_text.append(f'{len(ref_cluster[0]["discarded"])} hits were '
                            f'discarded from cluster {i} for {cpname} as '
                            'a result of refinement by retention time.')
            
            #Refine by the UV maxima of the hits
            ref2_cluster = refineClusterByUV(ref_cluster[0], UVdatafound, comment_text)
            comment_text = ref2_cluster[1]
            initial_length = len(ref_cluster[0]["orange"])
            comment_text.append(f'{len(ref2_cluster[0]["orange"]) - initial_length}'
                                f' hits were marked as tentative from cluster {i} '
                                f'for {cpname} as a result of refinement by UV maxima.') 
            
            #Refine using the mass confidence of the hits
            ref3_cluster = refineClusterByMassConf(ref2_cluster[0], comment_text)
            comment_text = ref3_cluster[1]
            initial_length = len(ref2_cluster[0]["orange"])
            comment_text.append(f'{len(ref3_cluster[0]["orange"]) - initial_length}'
                                f' hits were marked as tentative from cluster {i} '
                                f'for {cpname} as a result of refinement by mass confidence.')
            
            #Put the final refined cluster into the refined_clusters list
            refined_clusters.append(ref3_cluster[0])
        

        #If there is more than one cluster, compare their mean mass confidences
        #Remove any where that mean is less than the maximum mean * massConf_threshold
        if len(refined_clusters) > 1:
            after_ref = selectClusterByMassConf(refined_clusters)
            refined_clusters = after_ref[0]
            output["discarded_by_cluster"] += after_ref[1]
            comment_text.append(f'{len(after_ref[1])} clusters were discarded after refinement'
                                 ' by mass confidence. ')
            
        #If there is still more than one cluster, compare their size.
        #Take only those which are larger than (max_size * cluster_size_threshold)
        if len(refined_clusters) > 1:
            after_ref = selectClusterBySize(refined_clusters)
            refined_clusters = after_ref[0]
            output["discarded_by_cluster"] += after_ref[1]
            comment_text.append(f'{len(after_ref[1])} clusters were discarded after refinement'
                                 ' by size comparison. ')
        
        #If there is still more than one cluster, the peak used in the heatmap
        #will be selected based only on peakarea. Inform the user that 
        #manual selection is required. 
        if len(refined_clusters) > 1:
            comment_text.append(f'<strong>High Priority: More than one cluster was found for {cpname}. ' 
                                          'User MUST perform their own analysis.</strong>')

        
        #Convert the format of the data so that it is indexed
        #by the well it's in        
        cluster_by_well = {}
        for cluster in refined_clusters:
            for peak in cluster["green"]:
                if peak["well"] not in cluster_by_well:
                    cluster_by_well[peak["well"]] = {
                            "green": [],
                            "orange": [],
                            "discarded": [],
                            }
            
                cluster_by_well[peak["well"]]["green"].append(peak)
                
            for peak in cluster["orange"]:
                if peak["well"] not in cluster_by_well:
                    cluster_by_well[peak["well"]] = {
                            "green": [],
                            "orange": [],
                            "discarded": [],
                            }
            
                cluster_by_well[peak["well"]]["orange"].append(peak)
            
            for peak in cluster["discarded"]:
                if peak["well"] not in cluster_by_well:
                    cluster_by_well[peak["well"]] = {
                            "green": [],
                            "orange": [],
                            "discarded": [],
                            }
            
                cluster_by_well[peak["well"]]["discarded"].append(peak)
        
        #For each well, select which peaks are green and 
        #which should be discarded. 
        #This is done using the results of the above prioritisation
        #Peaks under "green" are chosen first. If none are available for that
        #well, a orange peak is used instead. It is important to note that
        #even a orange peak must conform to tight retentionTime criteria.
        for key, well in cluster_by_well.items():
            output["discarded"] += well["discarded"]
            if len(well["green"]) > 0:
                if len(well["green"]) == 1:
                    output["green"].append(well["green"][0])
                else:
                    id_max = options.mass_or_area
                    
                    max_val = max(peak[id_max] for peak in well["green"])
                    peak_added = False
                    
                    #Sort the peaks by their peak area, so that if two
                    #peaks have the same mass_conf when options.mass_or_area is
                    #equal to "mass_conf", the largest peak is selected in preference. 

                    for peak in sorted(well["green"], key = lambda x: x["area"], reverse=True):
                        if peak_added == False and peak[id_max] == max_val:
                            output["green"].append(peak)
                            comment_text.append(f'Largest {id_max} selected in preference to others available for well '
                                                f'{getUserReadableWell(peak["well"])}.')
                            peak_added = True
                        else:
                            output["discarded"].append(peak)
                            comment_text.append(f'The peak at {peak["time"]} in well {getUserReadableWell(peak["well"])} '
                                                      f'was discarded as it had a smaller {id_max} than an '
                                                      f'equally likely alternative.')
                
                if len(well["orange"]) > 0:
                    output["discarded"] += well["orange"]
                    for peak in well["orange"]:
                        comment_text.append(f'The peak at {peak["time"]} in well '
                                    f'{getUserReadableWell(peak["well"])} for {cpname} was discarded '
                                    f'because a better match was found.')
            elif len(well["orange"]) == 1:
                peak = well["orange"][0]
                output["green"].append(well["orange"][0])
                comment_text.append(f'<strong>The tentative peak at {peak["time"]} in well {getUserReadableWell(peak["well"])} was '
                                          f'used as there was no better option. User should check this well.</strong>')
            elif len(well["orange"]) > 1:
                id_max = options.mass_or_area
                max_val = max(peak[id_max] for peak in well["orange"])
                peak_added = False

                #Sort the peaks by their peak area, so that if two
                #peaks have the same mass_conf when options.mass_or_area is
                #equal to "mass_conf", the largest peak is selected in preference. 

                for peak in sorted(well["orange"], key = lambda x: x["area"], reverse=True):
                    if peak_added == False and peak[id_max] == max_val:
                        output["green"].append(peak)
                        comment_text.append(f'Largest {id_max} selected in preference to other tentative hits available for well '
                                            f'{getUserReadableWell(peak["well"])}.')
                        peak_added = True
                    else:
                        output["discarded"].append(peak)
                        comment_text.append(f'<strong>The tentative peak at {peak["time"]} in well {getUserReadableWell(peak["well"])} '
                                                      f'was discarded as it had a smaller {id_max} than an '
                                                      f'equally likely alternative.</strong>')
        
        mean_time = sum([i["time"] for i in output["green"]]) / len(output["green"])
                
    elif len(peakList) > 0:
        if not options.validate == "True":
            comment_text.append("Validation was not performed as requested by the user.")
        else:
            comment_text.append(f'Validation was not performed for {cpname} as '
                                'there were insufficient hits.')
        #organise the peaks by the well they occupy
        cluster_by_well = {}
        
        for cluster in clusters:
            for peak in cluster:
                if peak["well"] not in cluster_by_well:
                    cluster_by_well[peak["well"]] = []
                cluster_by_well[peak["well"]].append(peak)
     
        logging.info(f'Peaks were chosen for {cpname} by {options.mass_or_area}.')
        for index, well in cluster_by_well.items():
            if len(well) > 1:
                id_max = options.mass_or_area
                max_val = max(well, key = lambda x: x[id_max])[id_max]

                #Sort the peaks by their peak area, so that if two
                #peaks have the same mass_conf when options.mass_or_area is
                #equal to "mass_conf", the largest peak is selected in preference. 
                peak_added = False
                for peak in sorted(well, key = lambda x: x["area"], reverse=True):
                    
                    if peak_added == False and peak[id_max] == max_val:
                        output["green"].append(peak)
                        comment_text.append(f'The hit at {peak["time"]} for {cpname} in well'
                                            f' {getUserReadableWell(peak["well"])} was selected as it had the largest'
                                            f' {id_max}.')
                        peak_added = True
                    else:
                        output["discarded"].append(peak)
                        comment_text.append(f'The hit at {peak["time"]} for {cpname} in well'
                                            f' {getUserReadableWell(peak["well"])} was discarded as'
                                            f' it did not have the largest {id_max} of all hits found.')
            else:
                output["green"].append(well[0])


        mean_time = sum([i["time"] for i in output["green"]]) / len(output["green"])
           
    return [output, comment_text, cluster_bands]
        
def removeDupAssigns(compoundDF, internalSTD, SMs, products, by_products):
    """
    Checks each compound to ensure that no peak has been assigned 
    to two different compounds. If this has happened, the internalSTD
    (if present) takes first priority, limiting reactant is second priority,
    product third priority and finally a by-product is lowest priority.

    :param compoundDF: Pandas dataframe
    :param internalSTD: A string representing the name of the internal standard
    :param SMs: A list of starting material names
    :param products: A list of product names
    :param by_products: A list of by_product names
    
    :return: compoundDF as Pandas dataframe
    """
    def isSamePeak(hit1, row):
        overlap = False
        
        for peak in row["hits"]["green"]:
            if hit1["time"] == peak["time"] and hit1["well"] == peak["well"]:
                overlap = True
        return overlap
    
    items_to_remove = {}
    for index, row in compoundDF.iterrows():

        for bindex, brow in compoundDF.iterrows():
            
            if index != bindex:
                overlap = [i for i, j in enumerate(brow["hits"]["green"]) if isSamePeak(j, row)]
                #If the name is the internalSTD, we always remove the overlap from the brow
                if row["name"] == internalSTD:
                    for i in overlap:
                        if bindex in items_to_remove:
                            items_to_remove[bindex].append(i)
                        else:
                            items_to_remove[bindex] = [i]
                #Peaks should be assigned to the starting material if there is a conflict,
                #unless that peak belongs to the internal standard      
                elif row["name"] in SMs and brow["name"] != internalSTD:
                    for i in overlap:
                        if bindex in items_to_remove:
                            items_to_remove[bindex].append(i)
                        else:
                            items_to_remove[bindex] = [i]
                
                #Peaks should be assigned to the product if there is a conflict between product 
                #and by-product, as long as both compounds are expected in the well the peak 
                #belongs to.
                elif row["name"] in products and brow["name"] in by_products:
                    for i in overlap:

                        if bindex in items_to_remove:
                            items_to_remove[bindex].append(i)
                        else:
                            items_to_remove[bindex] = [i]
    
    #If a peak is scheduled for removal, the program should check to see if 
    #a previously discarded peak in the same cluster and well
    #should be promoted to replace it. 
    
    for index, to_remove in items_to_remove.items():
        
        hit_copy = compoundDF.loc[index, "hits"]

        new_green = []
        new_discarded = hit_copy["discarded"][:]
        new_dbc = hit_copy["discarded_by_cluster"][:]
        new_comments = compoundDF.loc[index, "comments"][:]

        for i, j in enumerate(hit_copy["green"]):
            if i not in to_remove:
                new_green.append(j)
            else:
            
                new_discarded.append(j)
                new_comments.append(f'<strong>The peak in well {getUserReadableWell(j["well"])} at {j["time"]} minutes was deselected '
                                    'because it was already assigned to the SM, internalSTD or by-product.</strong>')
                
                #find a replacement peak if possible to do so from the same cluster and same well
                cluster = j["cluster"]
                well = j["well"]
                poss_replace = [peak for peak in hit_copy["discarded"] if peak["cluster"] == cluster and peak["well"] == well]
                if len(poss_replace) > 0:
                    id_max = options.mass_or_area
                    best_peak = max(poss_replace, key = lambda x: x[id_max])
                    new_green.append(best_peak)
                    new_comments.append(f'<strong>The peak in well {getUserReadableWell(well)} at {best_peak["time"]} minutes was promoted to replace the '
                                        f'peak that was deselected, and was selected for promotion as it had the largest '
                                        f'{id_max}.</strong>')                    

                
        compoundDF.at[index, "hits"] = {
            "green": new_green,
            "discarded": new_discarded,
            "discarded_by_cluster": new_dbc
        }
        compoundDF.at[index, "comments"] = new_comments 
    return compoundDF

def generateOutputTable(compoundDF, internalSTD, SMs, products, by_products, total_area_abs):
    """
    Reformats the validated hits into a pandas table ready for visualisation and export. 
    
    :param compoundDF: The pandas datatable containing all compounds with their respective hits.
    :param internalSTD: the name of the internalSTD
    :param SMs: a list of indices for the starting materials
    :param products: a list of indices for the products
    :param by_products: a list of indices for the by-products
    :param total_area_abs: A float corresponding to the sum of all peak_area_absolutes
    
    :return: A Pandas table named outputTable
    
    """
    
    outputTable = {}
    #generate the structure of the table so that it is independant of the 
    #hits that are found

    for i in range(1, options.plate_row_no * options.plate_col_no + 1):
        well_id = getUserReadableWell(i)
        goingIn = {
            "well_no": i,
            "Well": well_id,
            "SMILES": "",
            "canonSMILES": "",
            "SMarea": 0,
            "SMareaAbs": 0,
            "Parea": 0,
            "PareaAbs": 0,
            "STDarea": 0,
            "STDareaAbs": 0,
            "Uarea": 100,
            "UareaAbs": total_area_abs[i],
            "corrSMarea": 0,
            "corrParea": 0,
            "corrSTDarea": 0,
            "P/SM+P": 0,
            "P/STD":0
            }
        for by_prod in by_products:            
            name = f'{by_prod}area'
            corr_name = f'corr{by_prod}area'
            abs_name = f'{by_prod}areaAbs'
            name_std = f'{by_prod}/STD'

            goingIn[name] = 0
            goingIn[abs_name] = 0
            goingIn[corr_name] = 0
            goingIn[name_std] = 0   

        outputTable[i] = goingIn  

    for index, row in compoundDF.iterrows():
        #Add the product SMILES to the outputTable
        if row["name"] in products:
            for location in row["locations"]:
                outputTable[location]["canonSMILES"] = index
                outputTable[location]["SMILES"] = row["g_smiles"]
        
        
        if len(row["hits"]["green"]) != 0:
            max_area = max([hit["area"] for hit in row["hits"]["green"]])
            
            for hit in row["hits"]["green"]:
                            
                if row["name"] in SMs:
                    outputTable[hit["well"]]["SMarea"] = hit["area"]
                    outputTable[hit["well"]]["SMareaAbs"] = hit["areaAbs"]
                    outputTable[hit["well"]]["corrSMarea"] = hit["area"] / max_area
                    outputTable[hit["well"]]["Uarea"] = outputTable[hit["well"]]["Uarea"] - hit["area"]
                    outputTable[hit["well"]]["UareaAbs"] = outputTable[hit["well"]]["UareaAbs"] - hit["areaAbs"]
                elif row["name"] == internalSTD:
                    outputTable[hit["well"]]["STDarea"] = hit["area"]
                    outputTable[hit["well"]]["STDareaAbs"] = hit["areaAbs"]
                    outputTable[hit["well"]]["corrSTDarea"] = hit["area"] / max_area
                    outputTable[hit["well"]]["Uarea"] = outputTable[hit["well"]]["Uarea"] - hit["area"]
                    outputTable[hit["well"]]["UareaAbs"] = outputTable[hit["well"]]["UareaAbs"] - hit["areaAbs"]
                elif row["name"] in products:
                    outputTable[hit["well"]]["Parea"] = hit["area"]
                    outputTable[hit["well"]]["PareaAbs"] = hit["areaAbs"]
                    outputTable[hit["well"]]["corrParea"] = hit["area"] / max_area
                    outputTable[hit["well"]]["Uarea"] = outputTable[hit["well"]]["Uarea"] - hit["area"]
                    outputTable[hit["well"]]["UareaAbs"] = outputTable[hit["well"]]["UareaAbs"] - hit["areaAbs"]

                    

                elif row["name"] in by_products:
                    name = f'{row["name"]}area'
                    corr_name = f'corr{row["name"]}area'
                    abs_name = f'{row["name"]}areaAbs'
                    outputTable[hit["well"]][name] = hit["area"]
                    outputTable[hit["well"]][abs_name] = hit["areaAbs"]
                    outputTable[hit["well"]][corr_name] = hit["area"] / max_area
                    outputTable[hit["well"]]["Uarea"] = outputTable[hit["well"]]["Uarea"] - hit["area"]
                    outputTable[hit["well"]]["UareaAbs"] = outputTable[hit["well"]]["UareaAbs"] - hit["areaAbs"]
    
    for i in range(1, options.plate_row_no * options.plate_col_no + 1):

        well = outputTable[i]
        if well["Uarea"] < 0: 
            well["Uarea"] = 0
        if well["PareaAbs"] != 0 or well["SMareaAbs"] != 0:
            well["P/SM+P"] = round(well["PareaAbs"] / (well["SMareaAbs"] + well["PareaAbs"]), 2)

        if well["STDareaAbs"] != 0:
            
            if well["PareaAbs"] != 0:
                well["P/STD"] = round(well["PareaAbs"] / well["STDareaAbs"], 2)
            for by_prod in by_products:
                name_std = f'{by_prod}/STD'
                name = f'{by_prod}areaAbs'
                well[name_std] = round(well[name] / well["STDareaAbs"], 2)

    #generate pandas dataframe, and sort it on the well_no (index)
    dataframe = pd.DataFrame.from_dict(outputTable, orient="index")
    dataframe.sort_index(inplace=True)
    
    #Generate composite metrics like P/STD, and add them to the dataframe. 
    corrP_SM = {}
    corrP_STD = {}
    for index, row in compoundDF.iterrows():
        if row["name"] in products:
            new_slice = dataframe.loc[(dataframe["well_no"].isin(row["locations"])), 
                                      ["well_no", "P/SM+P", "P/STD"]]
            
            max_val_P_SM = new_slice["P/SM+P"].max()
            max_val_P_STD = new_slice["P/STD"].max()
            
            for sindex, srow in new_slice.iterrows():
                if max_val_P_SM != 0:
                    corrP_SM[srow["well_no"]] = srow["P/SM+P"] / max_val_P_SM
                else:
                    corrP_SM[srow["well_no"]] = 0
                if max_val_P_STD != 0:
                    corrP_STD[srow["well_no"]] = srow["P/STD"] / max_val_P_STD
                else:
                    corrP_STD[srow["well_no"]] = 0 

        elif row["name"] in by_products:
            name_std = f'{by_prod}/STD'
            new_slice = dataframe.loc[(dataframe["well_no"].isin(row["locations"])), 
                                      ["well_no", name_std]]
            max_val_P_STD = new_slice[name_std].max()
            new_column = {}
                 
            for sindex, srow in new_slice.iterrows():
                if max_val_P_STD != 0:
                    new_column[srow["well_no"]] = srow[name_std] / max_val_P_STD
                else:
                    new_column[srow["well_no"]] = 0
            going_in = []
            for i in range(options.plate_row_no * options.plate_col_no):
                if i in new_column:
                    going_in.append(new_column[i])
                else:
                    going_in.append(0)
            dataframe[f'corr{by_prod}/STD'] = going_in
    
    #Convert the above data to a list so that it can be added
    #as a column to the dataframe. 
    going_in_sm = []
    going_in_std = []

    for i in range(options.plate_row_no * options.plate_col_no):
        if (i+1) in corrP_SM:
            going_in_sm.append(corrP_SM[i+1])
            going_in_std.append(corrP_STD[i+1])
        else:
            going_in_sm.append(0)
            going_in_std.append(0)
    
    dataframe["corrP/SM+P"] = going_in_sm    
    dataframe["corrP/STD"] = going_in_std

    return dataframe

def findOverlap(dataTable, well, time):
    """
    Provided with the dataTable and a retention time+well of interest, 
    this function will look for any overlapping peaks, and return the retention 
    times of those peaks to the user. 
    
    :param dataTable: A dictionary of list of dicts, indexed by well
    :param well: An integer to describe the well number
    :param time: The retention time, as a float, of the peak of interest
    
    :return: A list of overlaps
    """
    
    overlaps = []
    
    #Sort the peaks in that well by their retention time
    peaks = sorted(dataTable[well], key = lambda x: x["time"])
    
    #Select the peak index for the compound of interest
    peak_index = [index for index, i in enumerate(peaks) if i["time"] == time]
    
    if len(peak_index) > 0:
        index = peak_index[0]
        
        #Find if the peak prior has the same end time as this peak's start time
        if index != 0: 
            if peaks[index-1]["pEnd"] == peaks[index]["pStart"]:
                overlaps.append(peaks[index-1]["time"])
        #Find if the peak after has the same start time as this peak's end time
        if index != len(peaks)-1:
            if peaks[index]["pEnd"] == peaks[index+1]["pStart"]:
                overlaps.append(peaks[index+1]["time"])
        
    return overlaps

def findPotentialConflicts(compoundDF):
    """
    Find any products which are in danger of overlapping with other compounds that
    are expected in the same well. 
    This function is useful as the program may not detect SM (etc) for that peak/well
    because the mass_conf is too low, but the user should still be warned that there 
    may be a problem. 

    :param compoundDF: Pandas datatable for compounds
    :return: A text output for that compound. 
    """

    output = []

    for index, row in compoundDF.iterrows():
        
        new_slice = compoundDF.loc[(((compoundDF["time"] - row["time"]).abs() < 0.02)
                                    & (compoundDF["name"] != row["name"])
                                    & (compoundDF["time"] != 0))]
        close_compounds = []
        for bindex, brow in new_slice.iterrows():
            overlap = len([i for i in row["locations"] if i in brow["locations"]])
            if overlap > 0:
                close_compounds.append(brow["name"])
        if len(close_compounds) > 0:
            text = "<strong>The compound was found to have a similar retention time as "
            for cpindex in range(len(close_compounds)):
                if cpindex != len(close_compounds)-1:
                    text = f'{text}{close_compounds[cpindex]}, '
                elif len(close_compounds) > 1:
                    text = f'{text}and {close_compounds[cpindex]}. '
                else:
                    text = f'{text}{close_compounds[cpindex]}. '
            text = text + "The user should check this result manually.</strong>"
            output.append(text)
        else:
            text = "No potential conflicts found."
            output.append(text)

    return output        

def plotHitValidationGraph(cpname, validatedHits, save_dir, cluster_bands):
    """
    Plots all the hit peaks in a scatter graph of peaktime vs
    well, colour coded by whether the hit was included or discarded
    from the final output. 
    
    :param cpname: a string for the name of the compound of interest
    :param validatedHits: a dict, where each header contains a list of dicts
    :param save_dir: File directory for where to save the matplotlib figure
    :param cluster_bands: A list of the average retention times for each cluster.
    
    :return: jpg of the hit validation graph saved to output directory 
    """  
    
    #Find the largest peak area, and the minimum/maximum retention times
    #of the hits to be plotted
    allHits = validatedHits["green"] + validatedHits["discarded"]
    for cluster in validatedHits["discarded_by_cluster"]:
        for key, value in cluster.items():
                allHits = allHits + value
    
    max_size = max(allHits, key = lambda x: x["area"])["area"] if len(allHits) > 0 else 0
    max_time = max(allHits, key = lambda x: x["time"])["time"] if len(allHits) > 0 else 0
    min_time = min(allHits, key = lambda x: x["time"])["time"] if len(allHits) > 0 else 0

    
    #Used to calculate the size of each hit on the graph
    if max_size > 0:
        size_ratio = 10 / max_size

        #A list of possible markers used for plotting each hit, based on which
        #cluster that hit was assigned to (see validateHits fn)
        t_markers = ["o", "v", "^", "<", ">", "1", "s", "p", "P", "+", 
                     "x", "D", "2", "3", "4", "8", "h", "H", "d", "_"]

        #Determine the right scale for the y-axis now to avoid scaling issues 
        #between different points in the scatter graph            
        y_lims = [0, 0]
        if max_time - min_time < 0.1:
            mid_point = (max_time + min_time) / 2
            y_lims[0] = mid_point - 0.05
            y_lims[1] = mid_point + 0.05
        else:
            y_lims[0] = min_time
            y_lims[1] = max_time
        
        fig, ax = plt.subplots()  

        #Set the axis ranges   
        total_wells = options.plate_col_no * options.plate_row_no
        ax.set_xlim(-total_wells*0.2, total_wells*1.1)
        ax.set_ylim(y_lims[0] - (y_lims[1]-y_lims[0])/10, y_lims[1] + (y_lims[1]-y_lims[0])/10)
        
        #Plot all hits that were used in the analysis
        for hit in validatedHits["green"]:
            if hit["cluster"] > len(t_markers)-1:
                thismarker = "o"
            else:
                thismarker = t_markers[hit["cluster"]]
            ax.scatter(hit["well"], hit["time"], s = math.ceil(hit["area"] * size_ratio)**2, 
                        marker = thismarker, color="black")
        
        #Plot all hits that were discarded in favour of a better hit
        for hit in validatedHits["discarded"]:
            if hit["cluster"] > len(t_markers)-1:
                thismarker = "o"
            else:
                thismarker = t_markers[hit["cluster"]]
            ax.scatter(hit["well"], hit["time"], s = math.ceil(hit["area"] * size_ratio)**2, 
                        marker = thismarker, color="red")
        
        #Plot all hits that were discarded by cluster
        for cluster in validatedHits["discarded_by_cluster"]:   
            for key, value in cluster.items():
                for hit in value:
                    if hit["cluster"] > len(t_markers)-1:
                        thismarker = "o"
                    else:
                        thismarker = t_markers[hit["cluster"]]
                    ax.scatter(hit["well"], hit["time"], s = math.ceil(hit["area"] * size_ratio)**2, 
                                marker = thismarker, color="red")
        
        #Annotate the graph with the average retention time of each cluster
        for i in range(len(cluster_bands)):
            x = -total_wells*0.17
            y = cluster_bands[i]
            ax.annotate(f'C{i}: {cluster_bands[i]} min.', [x, y])

        #Label the graph and axes, and save to the output directory.   
        plt.title(cpname)
        plt.xlabel("Well")
        plt.ylabel("Retention Time /min")
        
        plt.savefig(f'{save_dir}graphs/hits-{cpname}.jpg', format="jpg")
        plt.close()

    


def plotChroma(cpname, wellno, trace, pStart, pEnd, annotate_peaks, save_dir, 
               ms_plus, ms_minus, mass1):
    """
    Plots the LCMS trace with labels for compounds found for a specific well,
    and highlights a particular peak of interest, providing m/z data for that
    peak. 

    :param cpname: a string of the compound name
    :param wellno: an integer representing the well
    :param trace: a list [x-values, y-values] to plot the Uv chromatogram 
    :param pStart: a float for the time where a specific peak begins
    :param pEnd: a float for the time where a specific peak ends. 
    :param annotate_peaks: a list of dictionaries for peaks to annotate
    :param save_dir: a string for the output directory
    :param ms_plus: a list [x-values, y-values] for MS+ spectrometric data
    :param ms_minus: a list [x-values, y-values] for MS- spectrometric data
    :param mass1: the isotopic mass of a compound, to which +1/-1 should be added
        to get to an expected observed mass (typically parent isotopic mass)
    
    :return: jpg of the chromatogram saved to output directory
    """
    fig, (a0, a1, a2) = plt.subplots(3, 1, gridspec_kw={'height_ratios': [2, 2, 6]})
    
    #Function to plot the mass spectrometric data
    def plotMS(axes, data, title):
        x = []
        y = []
        for i in data:
            x.append(i[0])
            y.append(i[1])
        
        axes.bar(x, y, width = 2)
        axes.set_title(title)
        axes.set_xlim(0, 1000)
        axes.set_ylim(0, 200)
        
        last_annotation = 0
        for peak in data:
            if peak[1] > 20:
                if math.isclose(peak[0], last_annotation, abs_tol = 3):
                    axes.annotate(peak[0], [peak[0]+15, peak[1]+20], ha="center", 
                                  va="bottom", rotation=90, size = 8)
                    axes.arrow(peak[0]+10, peak[1]+18, -10, -8)
                else:
                    axes.annotate(peak[0], [peak[0]+1, peak[1]], ha="center", 
                                  va="bottom", rotation=90, size = 8)
                last_annotation = peak[0]
    
    #Plot both MS- and MS+    
    plotMS(a0, ms_minus, "MS-")
    plotMS(a1, ms_plus, "MS+")

    #Plot the UV chromatogram
    a2.plot(trace[0], trace[1])

    #label the graph and axes
    label = getUserReadableWell(wellno)
    fig.suptitle(f'{cpname} ({mass1}): Well {label}')
    a2.set_xlabel("Time /min")
    a2.set_ylabel("AUs")

    #Find the index of the start and end times of the peak
    #that should be highlighted. 
    if pStart != 0 and pEnd != 0:
        x_index_start = min(enumerate(trace[0]), 
                        key = lambda x: abs(x[1]-pStart))[0]
        x_index_end = min(enumerate(trace[0]), 
                        key = lambda x: abs(x[1]-pEnd))[0]

        diff = trace[1][x_index_end] - trace[1][x_index_start]
        
        #Generate a second curve which can be used to specifically fill the hit
        #well. 
        second_curve = []
        for index, value in enumerate(trace[0]):
            #Calculate where the baseline value should be based on the relative
            #heights at the start and end of the peak. 
            if index >= x_index_start and index <= x_index_end:
                new_val = (trace[1][x_index_start] 
                    + diff*((index - x_index_start)/(x_index_end - x_index_start)))
                if new_val > trace[1][index]:
                    second_curve.append(trace[1][index])
                else:
                    second_curve.append(new_val)
            #For all x-values that don't fill inside the hit peak, 
            #the second curve should match the LCMS trace so that these areas don't
            #get filled. 
            else:
                second_curve.append(trace[1][index])
            
    #Annotate the peaks that were matched to a compound using an arrow
    #and dynamic positioning for clarity
    annotate_peaks.sort(key = lambda x: x["time"])
    last_x_position = 0
    max_x = max(trace[0])
    for index, i in enumerate(annotate_peaks):

        x_index = min(enumerate(trace[0]), 
                      key = lambda x: abs(x[1]-i["time"]))[0]
        if index < len(annotate_peaks) / 2:
            h_align = "right"
        else:
            h_align = "left"
        if index == 0:
            offset = -(max_x / 40)
        elif index == len(annotate_peaks)-1:
            offset = (max_x / 40)
        else:
            offset = 0
        arrow_props = dict(arrowstyle="->", connectionstyle="arc3")
        if math.isclose(last_x_position, i["time"]+offset, abs_tol = (max_x / 20)) and index != len(annotate_peaks)-1:
            if h_align == "left":
                h_align = "right"
            else:
                h_align = "left"

        if index % 2 == 0:
            plt.annotate(f'{i["cpname"]}: {i["time"]}', xy = (i["time"], trace[1][x_index]), 
                        xytext = (i["time"]+offset, 120), horizontalalignment=h_align,
                        arrowprops = arrow_props)
        else:
            plt.annotate(f'{i["cpname"]}: {i["time"]}', xy = (i["time"], trace[1][x_index]), 
                        xytext = (i["time"]+offset, 110), horizontalalignment=h_align,
                        arrowprops = arrow_props)
        last_x_position = i["time"]+offset
    min_y = min(trace[1])
    
    #Set axis limits
    a2.set_ylim(min_y - 10, 130)
    a2.set_xlim(-0.1, max_x+0.1)

    #Fill to highlight the peak of interest 
    if pStart != 0 and pEnd != 0:
        a2.fill_between(trace[0], second_curve, trace[1], color="red")

    #Set the number of tick points for each axis. 
    a2.locator_params(axis='y', nbins=6)
    a2.locator_params(axis='x', nbins=20)

    #Save graph to output directory. 
    plt.savefig(f'{save_dir}graphs/chroma-{cpname}-best.jpg', format="jpg")
    plt.close()

def plotPieCharts(zvalue, outputTable, save_dir, by_products):
    """
    Plots a set of pie charts for the full plate
    using the full dataset. The size of the pie chart is dependant on
    the value in the datatable for the column specified by zvalue
    
    :param zvalue: a string corresponding to the desired output metric (e.g. P/STD)
    :param outputTable: a pandas datatable
    :param save_dir: a string for the output directory
    :param by_products: a list of names of byproducts
    
    :return: jpg of the piecharts saved to output directory
    """

    #declare color palette
    palette = ["black", "brown", "red", "sienna",
                "peru", "orange", "gold", "olive", "lawngreen", "darkgreen", "lime", "aqua", 
                "steelblue", "slategray", "navy"]
    
    
    def buildPies(chart_type):
        """
        Define function to build a trellised pie chart. 

        :param chart_type: String defining whether this set of pie charts has fixed or variable 
                            diameters
        
        :return: jpg of trellised pie chart saved to output folder
        """

        #declare new subplots
        fig, axs = plt.subplots(options.plate_row_no, options.plate_col_no)
        
        #get maximum value, so diameter of all pies can be calculated in relation
        #to the best performing well. 
        max_val = max([row[zvalue] for index, row in outputTable.iterrows()])
        
        #format data and generate pie charts
        pies_baked = []
        if max_val != 0:
            
            for index, row in outputTable.iterrows():
                pies_baked.append(index)
                rowVal = math.floor((index-1) / options.plate_col_no)
                colVal = int((index-1) % options.plate_col_no)

                if chart_type == "fixed_width":
                    chart_size = 0.95
                else:
                    chart_size =  (row[zvalue] / max_val) * 0.95

                color_tracker = 0
                colors = []
                sum_byprod = sum([row[f'{by_prod}area'] for by_prod in by_products])
                
                total = row["SMarea"] + row["Parea"] + row["STDarea"] + sum_byprod
                #the LCMS machine can calculate a total area of >100% due to rounding errors
                #this needs to be corrected before the piechart is built
                if total > 100:
                    total = 100
                data = []
                if chart_size > 0: 
                    if total != 100:
                        data.append(100-total)
                        colors.append("grey")

                    if row["SMarea"] != 0:
                        data.append(math.floor(row["SMarea"]))
                        colors.append("cyan")

                    if row["Parea"] != 0:
                        data.append(math.floor(row["Parea"]))
                        colors.append("yellow")

                    if row["STDarea"] != 0:
                        data.append(math.floor(row["STDarea"]))
                        colors.append("green")
                    
                    for by_prod in by_products:
                        if math.floor(row[f'{by_prod}area']) != 0:
                            data.append(math.floor(row[f'{by_prod}area']))
                            colors.append(palette[color_tracker])
                        color_tracker = color_tracker + 1
                else:
                    #Set chart size to 0.01 to avoid a non-zero radius
                    #which causes matplotlib to fail. 
                    chart_size = 0.01
                    
                
                if options.plate_row_no == 1:
                    axs[colVal].pie(data, 
                            textprops={"size": "smaller"}, 
                            colors = colors,  
                            radius=chart_size,
                            normalize = True)
                elif options.plate_col_no == 1:
                    axs[rowVal].pie(data, 
                            textprops={"size": "smaller"}, 
                            colors = colors,  
                            radius=chart_size,
                            normalize = True)
                else:
                    axs[rowVal, colVal].pie(data, 
                            textprops={"size": "smaller"}, 
                            colors = colors,  
                            radius=chart_size,
                            normalize = True)
                
        #Remove any unused sections of the trellised pie chart
        #graph so that the visualisation is neater.         
        for i in range(1, options.plate_col_no * (options.plate_row_no) + 1):
            rowVal = math.floor((i-1) / options.plate_col_no)
            colVal = int((i-1) % options.plate_col_no)
            if i not in pies_baked:
                
                #matplotlib axes drop a dimension 
                #when that dimension is length = 1. 
                #Adjust code to match
                if options.plate_row_no == 1:
                    plt.delaxes(axs[colVal])
                elif options.plate_col_no == 1:
                    plt.delaxes(axs[rowVal])
                else:
                    plt.delaxes(axs[rowVal, colVal])
           

        #Add a key to the graph 
        lines = [Line2D([0], [0], color="grey", lw=4), 
                Line2D([0], [0], color="cyan", lw=4),
                Line2D([0], [0], color="yellow", lw=4),
                Line2D([0], [0], color="green", lw=4),]
        labels = ["Untagged", "Reactant", "Product", "InternalSTD"]
        for i in range(len(by_products)):
            lines.append(Line2D([0], [0], color=palette[i], lw=4))
            labels.append(by_products[i])
        if options.plate_row_no == 1:
            axs[options.plate_col_no-1].legend(lines,
                    labels, loc="lower center",
                    bbox_to_anchor=(1,1), ncol=math.ceil(len(labels)/3))
        elif options.plate_col_no == 1:
            axs[options.plate_row_no-1].legend(lines,
                    labels, loc="upper left",
                    bbox_to_anchor=(1,1), ncol=math.ceil(len(labels)/3))
        else:
            axs[options.plate_row_no-1, options.plate_col_no-1].legend(lines,
                    labels, loc="upper right",
                    bbox_to_anchor=(1,0), ncol=math.ceil(len(labels)/3))
            
        #Add titles to the graphs
        if chart_type == "fixed_width":
            fig.suptitle("Fixed Diameter Trellised Pie Charts", y=0.9)
        else:
            fig.suptitle(f'Trellised Pie Charts Sized by {zvalue}', y=0.9)
        
        #Save graph to output directory
        plt.savefig(f'{save_dir}graphs/piecharts_{chart_type}.jpg', format="jpg")
        plt.close()
    
    buildPies("fixed_width")
    buildPies("variable_width")

def plotHeatmaps(outputTable, save_dir):
    """
    Plots and saves heatmaps for the full dataset
    
    :param outputTable: a pandas datatable
    :param save_dir: a string for the output directory
    
    :return: jpg of the heatmap saved to output directory
    """

    
    zvalues = {
        "Parea": "Parea", 
        "conversion": "P/SM+P", 
        "ratio_to_IS": "P/STD",
        "corrParea": "corrParea", 
        "corrected_conversion": "corrP/SM+P",
        "corrected_ratio_to_IS": "corrP/STD"
    }
    for key, zvalue in zvalues.items():
        
        returnVal = []
        labels = []
        #convert output table to 2-D list that can be used to plot heatmap. 
        for index, row in outputTable.iterrows():
            rowVal = int(math.floor((index-1) / options.plate_col_no))
            colVal = int((index-1) % options.plate_col_no)
            if colVal == 0:
                returnVal.append([])
                labels.append([])
            if row["SampleID"] == "No Data.":
                labels[rowVal].append("-")
                returnVal[rowVal].append(0)
            else:
                if row[zvalue] > 1:
                    labels[rowVal].append(math.floor(row[zvalue]))
                elif row[zvalue] >= 100:
                    labels[rowVal].append(100)
                elif row[zvalue] == 0:
                    labels[rowVal].append("0")
                else:
                    labels[rowVal].append(round(row[zvalue],2))

                returnVal[rowVal].append(row[zvalue])

        #Configure heatmap
        pdTable = pd.DataFrame(returnVal)
        xLabels = [i for i in range(1, options.plate_col_no+1)]
        yLabels = [chr(ord('@')+i) for i in range(1, options.plate_row_no+1)]
        
        ax = sns.heatmap(pdTable, xticklabels=xLabels, yticklabels=yLabels, cmap = "viridis",
                        annot = labels, cbar_kws={"label": zvalue}, fmt="")
        ax.xaxis.set_ticks_position("top")
        
        #Save heatmap to output directory
        plt.savefig(f'{save_dir}graphs/heatmap_{key}.jpg', format="jpg")
        plt.close()

def plotHistogram(dataframe, save_dir):
    """
    Takes the outputTable dataframe and generates a histogram for product %area
    to be saved in the output folder. 

    :param dataframe: A Pandas dataframe (AKA outputTable)
    :param save_dir: a string for the output directory
    
    :return: jpg of the histogram saved to output directory
    """

    #Configure histogram. Use colours for clear differentiation
    #between poor, medium and good results for that plate. 

    new_slice = dataframe.loc[dataframe["SampleID"] != "No Data."]
    fig, ax = plt.subplots()

    N, bins, patches = ax.hist(new_slice["Parea"], 50)
    for i in range(0, 16):
        patches[i].set_facecolor("blue")
    for i in range(16, 31):
        patches[i].set_facecolor("green")
    for i in range(31, 50):
        patches[i].set_facecolor("yellow")
    
    #Label axes and graph
    plt.xlabel("Product Percentage Area")
    plt.ylabel("Count")
    plt.title("Product Percentage Area over Plate")

    #Save graph to output directory
    plt.savefig(f'{save_dir}graphs/histogram.jpg', format="jpg")
    plt.close()

def plotDonut(dataframe, save_dir):
    """
    Takes the outputTable dataframe and generates a donut chart for product %area
    to be saved in the output folder.

    :param dataframe: A Pandas dataframe (AKA outputTable)
    :param save_dir: a string for the output directory
    
    :return: Saved jpg of the donut chart
    
    """
    
    #Bucket data into whether the plate gave high, medium, low product
    #purity, or whether no product was observed
    Parea_list = dataframe["Parea"].to_list()
    num_none_formed = len([i for i in Parea_list if i == 0])
    num_trace_formed = len([i for i in Parea_list if i > 0 and i <= 20])
    num_some_formed = len([i for i in Parea_list if i > 20 and i <=50])
    num_lots_formed = len([i for i in Parea_list if i > 50])

    #Format data and labels
    values = [num_none_formed, num_trace_formed, num_some_formed, num_lots_formed]
    colours = ["black", "blue", "green", "yellow"]
    labels = [f'None formed: {num_none_formed}', f'Poor: {num_trace_formed}',
        f'Medium: {num_some_formed}', f'Good: {num_lots_formed}']
    
    #Configure pie chart
    fig, ax = plt.subplots()
    ax.pie(values, colors = colours, labels = labels)

    #Create and add white circle to make it a donut graph
    centre_circle = plt.Circle((0,0), 0.7, fc="white")
    fig.gca().add_artist(centre_circle)

    #Label graph and save it to output directory
    plt.title("Product Percentage Area over Plate")
    plt.savefig(f'{save_dir}graphs/donut.jpg', format="jpg")
    plt.close()


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


def buildHTML(save_dir, compoundDF, all_compounds, analysis_name, times = {}):
    """
    Build a HTML output file using jinja2 and a html_template
    that is stored in the directory "templates". 
    
    :param save_dir: A string designating the output directory
    :param compoundDF: Pandas datatable containing all information on 
                        the compounds used for analysis
    :param all_compounds: a list of all compound names
    :param times: Optional parameter of a list of floats related to processing time
                    for each step of the analysis. 
    :return: HTML file saved to save_dir 
    """


    #Add key data for each compound out of the compoundDF
    env = Environment(
        loader = FileSystemLoader("templates")
    )

    template = env.get_template("html_template.html")

    cptablerows = []
    for index, row in compoundDF.iterrows():
        cptablerows.append(row)
    

    with open(f'{save_dir}/html_output.html', "w") as fo:
        fo.write(template.render(
            cpnames = all_compounds,
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
                       
                        min_peak_area = 0.5,
                        min_massconf_threshold = 10,
                        min_uv_threshold = 20,

                        min_no_of_wells = 5,

                        uv_match_threshold = 0.5,
                        uv_cluster_threshold = 0.5,
                        massconf_threshold = 0.5,
                        cluster_size_threshold = 0.8,
                           
                        plate_col_no = 12,
                        plate_row_no = 8,

                        points_per_trace = 500,
                        
                        mass_or_area = "mass_conf",
                        plot_type = "Parea", 
                        calc_higherions = "True",
                        calc_boc = "True",
                        
                        gen_csv = "True",
                        gen_zip = "False",

                        output = "output",

                        )
    
    parser.add_argument("input_rpt", help = "Input .rpt file for analysis")
    
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
                               "from the following options: Parea, corrParea, P/SM+P, P/STD, corrP/STD, corrP/SM+P.\n")
    
    parser.add_argument("-chi", "--calc_higherions", action="store", type=str, dest = "calc_higherions",
                        help = "Look for [M+2H]2+ and [M+3H]3+ to find hits and calculate the mass confidence"
                        " of a peak, True/False")
                        
    parser.add_argument("-cboc", "--calc_boc", action="store", type=str, dest = "calc_boc",
                        help = "Look for Boc degradation ions, True/False")
                        
    parser.add_argument("-g", "--generate_csv", action="store", type=str, dest = "gen_csv", 
                        help = "Choose to generate and save a CSV, True/False.\n")
    
    parser.add_argument("-z", "--generate_zip", action="store", type=str, dest = "gen_zip", 
                        help = "Choose to generate and save a zip file, True/False.\n")
    
    parser.add_argument("-n", "--name", action="store", type=str, dest = "analysis_name", 
                        help = "Choose a name for the analysis.\n")


    #Set options to global and parse arguments        
    times = {}
    global options
    options = parser.parse_args()
    
    #Set standard matplotlib graph size
    plt.rcParams["figure.figsize"] = (12, 6)
    
    #Ensure all required input options have a valid address
    root_names = [options.input_csv, options.input_rpt, options.output]
    
    for name in root_names:
        if name.startswith("/"):
            pass
        else:
            if name.startswith("./"):
                name.replace("./", os.getcwd())
            else:
                name = os.getcwd() + "/" + name    
    
    #Create the save directory if one isn't already present.
    global save_dir
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

    #Generate a generic HTML output file so that Pipeline Pilot displays
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

    #Check to ensure the root names of the rpt and csv files are correct, i.e. the files exist
    for name in [root_names[0], root_names[1]]:
        if not os.access(name,os.R_OK) or not os.path.isfile(name):
            logging.error(f'Input file {name} does not exist. Please use an appropriate root name.')
            with open(f'{save_dir}html_output.html', 'w') as f:
                f.write(f'Input file {name} does not exist. Please use an appropriate root name.')
                f.close()
            sys.exit(2)

    #Check to make sure the input files are in the .csv and .rpt format. 
    if options.input_csv[-4:].lower() != ".csv":
        logging.error(f'Plate map is not in the .csv format. Please use a .csv file format.')
        with open(f'{save_dir}html_output.html', 'w') as f:
            f.write(f'Plate map is not in the .csv format. Please use a .csv file format.')
            f.close()
        sys.exit(2)

    if options.input_rpt[-4:].lower() != ".rpt":
        logging.error(f'LCMS data is not in the .rpt format. Please use a .rpt file format.')
        with open(f'{save_dir}html_output.html', 'w') as f:
            f.write(f'LCMS data is not in the .rpt format. Please use a .rpt file format.')
            f.close()
        sys.exit(2)

    times["Initialise Script"] = time.perf_counter() - time_start 

    #Import the structure data from the comma-separated values platemap file provided
    #by the user to initiate the compoundDF. 
    pre_import = time.perf_counter()            
    [compoundDF, internalSTD, SMs, products, by_products] = importStructures(root_names[0])

    #Check that an internal standard was provided if the plot_type expects one
    if options.plot_type in ["P/STD", "corrP/STD"] and internalSTD == "":
        logging.error("No internal standard was provided. Halting program.")
        with open(f'{save_dir}html_output.html', 'w') as f:
            f.write("No internal standard was provided for an analysis type which required one.")
            f.close()
        sys.exit(2)

    #append all the compounds into a master list to make a 
    #pandas categorical for sorting

    all_compounds = [product for product in products]
    all_compounds = all_compounds + [byprod for byprod in by_products]
    all_compounds = all_compounds + [SM for SM in SMs]
    all_compounds.append(internalSTD)
    try: 
        compoundDF["name"] = pd.Categorical(compoundDF["name"], categories = all_compounds, ordered=True)
        compoundDF.sort_values(["name"], inplace=True)
    except:
        logging.error("The substrates could not be categorised. Check to ensure that duplicate compounds have not been assigned different labels "
                    "or that two by-products have not been given the same label.")
        with open(f'{save_dir}html_output.html', 'w') as f:
            f.write("The substrates could not be categorised. Check to ensure that duplicate compounds have not been assigned different labels "
                    "or that two by-products have not been given the same label.")
            f.close()
        sys.exit(2)

    #Log the time taken to complete these actions
    times["Import Compound File"] = time.perf_counter() - pre_import
    logging.info(f'{len(compoundDF.index)} compounds were imported.')

    
    #Convert the .rpt file into a useful format. The output is a dictionary,
    #indexed by the well on the plate, and each value is a list. 
    #each value in the list is a dictionary corresponding to each 
    #peak in the well. 
    #chroma is a list of lists containing the information required to 
    #generate the LCMS trace. 
    
    pre_rpt = time.perf_counter()
    [dataTable, chroma, sample_IDs] = convertRPT2Dict(root_names[1])
    times["Import RPT File"] = time.perf_counter() - pre_rpt 
    logging.info(f'{len(dataTable.items())} wells were found in the rpt.')

    #Determine the total areaAbs for each well, so that these values can be added to the outputTable
    total_area_abs = {}
    for i in range(1, options.plate_row_no * options.plate_col_no + 1):
        total_area_abs[i] = 0
    for index, well in dataTable.items():
        for peak in well:
            if index in total_area_abs:
                total_area_abs[index] = total_area_abs[index] + peak["areaAbs"]

    #For each compound, find all hits using the dataTable, validate them, and append 
    #them to a list ready for insertion into the compoundDF pandas dataframe
    pre_hit_and_val = time.perf_counter()
    validated_hits = []
    comment_text = []
    cluster_bands = []
    for index, row in compoundDF.iterrows():
        hits = findHits(row, dataTable)
        
        logging.debug(f'Hit finding was completed for {row["name"]}.')
        if len(hits) > 0:
            #Validate the hits (or just choose the best well based on
            #options.mass_or_area) using the name of the compound
            #the hits, and the expected retention time if one was provided
            #by the user. 
            validHits = validateHits(row["name"], hits, row["rt"])
            
            #Append the data to appropriate columns ready to add to the 
            #pandas datatable for compounds. 
            validated_hits.append(validHits[0])
            comment_text.append(validHits[1])
            cluster_bands.append(validHits[2])
            logging.debug(f'Validation was completed for {row["name"]}.')
        else:
            #if no hits were found, append the appropriate data
            #structure so that no error is incurred in the future. 
            placeholder = {
                    "green": [],
                    "discarded": [], 
                    "discarded_by_cluster": []
                    }
            validated_hits.append(placeholder)

            cluster_bands.append([])
            new_comments = [f'<strong>No hits were found for {row["name"]}.</strong>']
            new_comments.append(f'No validation was performed for {row["name"]} '
                                'as no hits were found.')
            comment_text.append(new_comments)

    

    #append the hits for each compound as a new column to the dataframe
    compoundDF["hits"] = validated_hits
    compoundDF["comments"] = comment_text
    compoundDF["cluster_bands"] = cluster_bands  
    times["Find Hits and Validate"] = time.perf_counter() - pre_hit_and_val
    
    #Check no peak has been assigned twice, and reassign if necessary
    pre_rem_dup = time.perf_counter()
    compoundDF = removeDupAssigns(compoundDF, internalSTD, SMs, products, by_products)
    logging.info(f'Duplicate assignments were removed.')
    times["Remove Duplicate Assignments"] = time.perf_counter() - pre_rem_dup

    #Generate the output table using validated hits
    pre_output_table = time.perf_counter()
    outputTable = generateOutputTable(compoundDF, internalSTD, SMs, products, by_products, total_area_abs)
    
    #Insert blanks in sample_IDs to match length of datatable
    for i in range(1, options.plate_row_no * options.plate_col_no +1):
        if i not in sample_IDs:
            sample_IDs[i] = "No Data."
    sorted_IDs = [value for key, value in sorted(sample_IDs.items(), key = lambda x: x[0])]
    outputTable.insert(1, "SampleID", sorted_IDs)
    logging.info(f'The output table was generated.')
    times["Generate Output Table"] = time.perf_counter() - pre_output_table
    
    
    #if the plot_type requires the internalSTD to be present, 
    #check to ensure a non-zero area of internalSTD is present in all wells
    warn_for_std = False
    if internalSTD != "":
        min_std = outputTable["STDarea"].min()
        if min_std == 0:
            warn_for_std = True
    
    #Plot graphs of retentiontime vs well for all hits, for each compound
    #During the same loop, gather information to push into the compoundDF. 
    pre_plotting = time.perf_counter()
    best_wells = []
    best_masses_plus = []
    best_masses_minus = []
    best_peak_areas = []
    best_times = [] 
    overlaps = []   
    
    for index, row in compoundDF.iterrows():
         
        #Plot hit validation graph for the compound
        plotHitValidationGraph(row["name"], row["hits"], save_dir, row["cluster_bands"])
        
        #Append a warning to the internalSTD if necessary
        if row["name"] == internalSTD and warn_for_std == True:
            compoundDF.at[index, "comments"].append("<strong>The internal standard was not found in all wells. "
                                    "Recommend re-running analysis using Parea or corrParea analysis type.</strong>")
                                    
        if len(row["hits"]["green"]) > 0:
            
            #find which well gives the best conditions:
            #For the product, this should reflect the metric the user selected.
            #e.g. Parea or corrP/STD 
            if row["name"] in products:
                new_slice = outputTable.loc[outputTable["well_no"].isin(row["locations"])]
                max_index = new_slice[options.plot_type].idxmax()
                best_well_id = new_slice["well_no"][max_index]
                best_well_poss = [i for i in row["hits"]["green"] if i["well"] == best_well_id]
                
                if len(best_well_poss) > 0:
                    best_well = best_well_poss[0]
                else:
                    best_well = max(row["hits"]["green"], key = lambda x: x["area"])
                    
            else:
                best_well = max(row["hits"]["green"], key = lambda x: x["area"])
            
            annotate_peaks = []
            for bindex, brow in compoundDF.iterrows():
                data = [i for i in brow["hits"]["green"] if i["well"] == best_well["well"]]
                if len(data) > 0:
                    annotate_peaks.append({
                            "cpname": brow["name"], 
                            "time": data[0]["time"]
                            })
            
            #Get the ms data for the relevant peak in the best_well
            ms_plus = [i["MS+"] for i in dataTable[best_well["well"]] 
                                if i["time"] == best_well["time"]][0]
            ms_minus = [i["MS-"] for i in dataTable[best_well["well"]] 
                                if i["time"] == best_well["time"]][0]
            try:  
                plotChroma(row["name"], best_well["well"], chroma[best_well["well"]], 
                        best_well["pStart"], best_well["pEnd"], 
                        annotate_peaks, save_dir, ms_plus, ms_minus, row["mass1"])
                logging.debug(f'A chromatogram plotted for {row["name"]} for well '
                                f'{best_well["well"]}.')
            except:
                logging.info("No chromatogram could be plotted. Data not found.")

            #Convert the well no. into a user readable format, then append it 
            #to an array for upload to the compoundDF. 
            best_wells.append(getUserReadableWell(best_well["well"]))
            
            #Prepare to add the best peak area, and a warning as to whether
            #the hit peak overlaps with another peak, to the compoundDF. 
            best_peak_areas.append(best_well["area"])
            best_masses_plus.append(best_well["mass+"])
            best_masses_minus.append(best_well["mass-"])
            best_times.append(round(best_well["time"], 2))
            is_overlap = findOverlap(dataTable, best_well["well"],
                                                    best_well["time"])
            if len(is_overlap) > 0:
                overlaps.append("<strong>Peak overlap detected!</strong>")
            else:
                overlaps.append("No peak overlap detected.")

            
        else:
            #Append appropriate placeholders to the data structure so that
            #there are no null values in the compoundDF. 
            if len(row["hits"]["discarded"]) == 1:

                annotate_peaks = []
                best_well = row["hits"]["discarded"][0]
                for bindex, brow in compoundDF.iterrows():
                    data = [i for i in brow["hits"]["green"] if i["well"] == best_well["well"]]
                    if len(data) > 0:
                        annotate_peaks.append({
                                "cpname": brow["name"], 
                                "time": data[0]["time"]
                                })
                
                #Get the ms data for the relevant peak in the best_well
                ms_plus = [i["MS+"] for i in dataTable[best_well["well"]] 
                                    if i["time"] == best_well["time"]][0]
                ms_minus = [i["MS-"] for i in dataTable[best_well["well"]] 
                                    if i["time"] == best_well["time"]][0]   
                try:
                    plotChroma(row["name"], best_well["well"], chroma[best_well["well"]], 
                            best_well["pStart"], best_well["pEnd"], 
                            annotate_peaks, save_dir, ms_plus, ms_minus, row["mass1"])
                except:
                    logging.info("No chromatogram could be plotted. Data not found.")
            
            elif len(row["locations"]) == 1:

                annotate_peaks = []
                only_well = row["locations"][0]
                for bindex, brow in compoundDF.iterrows():
                    data = [i for i in brow["hits"]["green"] if i["well"] == only_well]
                    if len(data) > 0:
                        annotate_peaks.append({
                                "cpname": brow["name"], 
                                "time": data[0]["time"]
                                })
                
                #Get the ms data for the relevant peak in the best_well
                ms_plus = []
                ms_minus = []
                try:  
                    plotChroma(row["name"], only_well, chroma[only_well], 0, 0, 
                            annotate_peaks, save_dir, ms_plus, ms_minus, row["mass1"])
                except:
                    logging.info("No chromatogram could be plotted. Data not found.")



            overlaps.append("N/A")
            best_wells.append("None found.")
            best_peak_areas.append(0)
            best_masses_plus.append("-")
            best_masses_minus.append("-")
            best_times.append(0)
            
    #Merge the new data into the compoundDF. 
    compoundDF["best_well"] = best_wells
    compoundDF["best_purity"] = best_peak_areas
    compoundDF["overlaps"] = overlaps
    compoundDF["mass+"] = best_masses_plus
    compoundDF["mass-"] = best_masses_minus
    compoundDF["time"] = best_times
    times["Plot Peaks"] = time.perf_counter() - pre_plotting
    
    #check all compounds for similar retention times to another compound
    #which may also appear in one of the same wells, to flag to the user
    #that a compound may be hiding underneath the peak of another. This 
    #additional check is useful in case the program didn't find both compounds
    #for a peak because the m/z for one compound was too small. 
    pre_checks = time.perf_counter()
    pot_conflicts = findPotentialConflicts(compoundDF)
    compoundDF["conflicts"] = pot_conflicts
    times["Find Potential Conlicts"] = time.perf_counter() - pre_checks
    logging.info('Potential conflicts were searched for all compounds.')

    #Return a heatmap to the user
    pre_heatmap = time.perf_counter()
    plotHeatmaps(outputTable, save_dir)  
    times["Generate Heatmap"] = time.perf_counter() - pre_heatmap
    logging.info(f'Heatmaps were generated successfully.')

    #Return a series of piecharts to the user, as long as it's not larger than
    #a 96 well plate. For larger plates, the piecharts become too small to be
    #useful. 
    if options.plate_row_no * options.plate_col_no < 97:
        pre_pie = time.perf_counter()
        plotPieCharts(options.plot_type, outputTable, save_dir, by_products)
        logging.info(f'A set of pie-charts was generated using {options.plot_type} '
                    f'as the index.')
        times["Generate Piechart"] = time.perf_counter() - pre_pie

    #Plot a histogram and donut chart of Parea as a measure of plate success. 
    pre_donut = time.perf_counter()
    plotHistogram(outputTable, save_dir)
    plotDonut(outputTable, save_dir)
    times["Generate Histogram and Donut"] = time.perf_counter() - pre_donut
    logging.info(f'A histogram and donut chaty was generated.')

    #Generate a set of PNG files to depict each compound
    
    for index, row in compoundDF.iterrows():
        generateMol(row["g_smiles"], row["name"], save_dir)

    #Generate the HTML output. 
    times["Total time"] = time.perf_counter() - pre_donut
    buildHTML(save_dir, compoundDF, all_compounds, options.analysis_name, times = times)
    logging.info('The HTML output was generated.')

    #Generate an csv of the output table.
    if options.gen_csv == "True":
        csv = outputTable.to_csv(f'{save_dir}outputTable.csv', index = False)
        newslice = compoundDF.loc[:, ["name", "g_smiles", "mass1", "mass2", "mass3",
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
    
    

    total_time = time.perf_counter() - time_start
    logging.info(f'The analysis was completed in {total_time} seconds.')
    print(f'The analysis was completed in {round(total_time, 2)} seconds.')

if __name__ == "__main__":
    try:
        print("")
        print("Running analysis....")
        main()
    except Exception:
        logging.exception("A fatal exception occurred. Contact administrator.")
