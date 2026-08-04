#External Libraries
import pandas as pd
from rdkit import Chem

#Internal Libaries
import plate_tools

class Output:

    """
    Takes the assigned peaks from a cpTable (see cpAssignment.py) and lcData (see, e.g. getWatersData.py) 
    and produces a Pandas Dataframe summarising the results for each well. Facilitates the creation of 
    heatmaps. 

    The following are descriptions of the Pandas dataframes that will be generated for this object:

    self.df = pd.DataFrame(columns =['well_no', 'Well', 'canonSMILES', 'SMarea', 'SMareaAbs', 'Parea', 'PareaAbs', 
    'STDarea', 'STDareaAbs', 'Uarea', 'UareaAbs', 'corrSMarea', 'corrParea', corrSTDarea', 'P/SM+P', 'P/STD'])
    self.aTable = pd.DataFrame(columns =['well', 'well_ID', 'Reaction_ID', 'analyte_smiles', 'analyte_class', 
                                'data_entry', 'data_type', 'data_units', 'data_description'])
    """
    
    def __init__(self, sample_IDs, plate_col_no, plate_row_no):
        """
        Initialise the Output object. 
        :param sample_IDs: a list of strings
        :param plate_col_no: integer, describing number of columns in plate
        :param plate_row_no: integer, describing number of rows in plate

        :return: The initialised Output object. 
        """
        self.plate_col_no = plate_col_no
        self.plate_row_no = plate_row_no
        self.sample_IDs = sample_IDs
        
        #Insert blanks in sample_IDs to match length of datatable
        for i in range(1, plate_row_no * plate_col_no +1):
            if i not in self.sample_IDs:
                self.sample_IDs[i] = "No Data."
        
    def generateOutputTable(self, cpTable, lcData):
        """
        Reformats the validated hits into a pandas table ready for visualisation and export. 

        :param cpTable: Pandas DataFrame, describing compounds and their assigned peaks
        :param lcData: Pandas DataFrame, describing the LC peak data. 

        :return: A Pandas table named outputTable

        """
        
        def addCompoundToTable(row):
            #Add the product SMILES to the outputTable
            if row["type"] == "Product":
                for location in row["locations"]:
                    outputTable[location]["canonSMILES"] = row["canonSMILES"]


            if len(row["final_result"]["green"]) != 0:
                max_area = max([lcData.at[i, "area"] for i in row["final_result"]["green"]])

                for i in row["final_result"]["green"]:
                    if row["type"] == "Reactant":
                        outputTable[lcData.at[i, "well"]]["SMarea"] = lcData.at[i, "area"]
                        outputTable[lcData.at[i, "well"]]["SMareaAbs"] = lcData.at[i, "areaAbs"]
                        outputTable[lcData.at[i, "well"]]["corrSMarea"] = lcData.at[i, "area"] / max_area
                        outputTable[lcData.at[i, "well"]]["Uarea"] = outputTable[lcData.at[i, "well"]]["Uarea"] - lcData.at[i, "area"]
                        outputTable[lcData.at[i, "well"]]["UareaAbs"] = outputTable[lcData.at[i, "well"]]["UareaAbs"] - lcData.at[i, "areaAbs"]
                    
                    elif row["type"] == "InternalSTD":
                        outputTable[lcData.at[i, "well"]]["STDarea"] = lcData.at[i, "area"]
                        outputTable[lcData.at[i, "well"]]["STDareaAbs"] = lcData.at[i, "areaAbs"]
                        outputTable[lcData.at[i, "well"]]["corrSTDarea"] = lcData.at[i, "area"] / max_area
                        outputTable[lcData.at[i, "well"]]["Uarea"] = outputTable[lcData.at[i, "well"]]["Uarea"] - lcData.at[i, "area"]
                        outputTable[lcData.at[i, "well"]]["UareaAbs"] = outputTable[lcData.at[i, "well"]]["UareaAbs"] - lcData.at[i, "areaAbs"]
                    
                    elif row["type"] == "Product":
                        outputTable[lcData.at[i, "well"]]["Parea"] = lcData.at[i, "area"]
                        outputTable[lcData.at[i, "well"]]["PareaAbs"] = lcData.at[i, "areaAbs"]
                        outputTable[lcData.at[i, "well"]]["corrParea"] = lcData.at[i, "area"] / max_area
                        outputTable[lcData.at[i, "well"]]["Uarea"] = outputTable[lcData.at[i, "well"]]["Uarea"] - lcData.at[i, "area"]
                        outputTable[lcData.at[i, "well"]]["UareaAbs"] = outputTable[lcData.at[i, "well"]]["UareaAbs"] - lcData.at[i, "areaAbs"]
                        
                    elif row["name"] in byproducts:
                        name = f'{row["name"]}area'
                        corr_name = f'corr{row["name"]}area'
                        abs_name = f'{row["name"]}areaAbs'
                        outputTable[lcData.at[i, "well"]][name] = lcData.at[i, "area"]
                        outputTable[lcData.at[i, "well"]][abs_name] = lcData.at[i, "areaAbs"]
                        outputTable[lcData.at[i, "well"]][corr_name] = lcData.at[i, "area"] / max_area
                        outputTable[lcData.at[i, "well"]]["Uarea"] = outputTable[lcData.at[i, "well"]]["Uarea"] - lcData.at[i, "area"]
                        outputTable[lcData.at[i, "well"]]["UareaAbs"] = outputTable[lcData.at[i, "well"]]["UareaAbs"] - lcData.at[i, "areaAbs"]
        
        def calculateNormalisedMetrics(row):
            #Calculate the corrP/STD value for each row
            
            max_P_STD = self.df.loc[self.df["canonSMILES"] == row["canonSMILES"], "P/STD"].max()
            if max_P_STD != 0:
                row["corrP/STD"] = row["P/STD"]/ max_P_STD
            else:
                row["corrP/STD"] = 0
            
            #Calculate the corrP/SM+P value for each row
            max_P_SM = self.df.loc[self.df["canonSMILES"] == row["canonSMILES"], "P/SM+P"].max()
            if max_P_SM != 0:
                row["corrP/SM+P"] = row["P/SM+P"]/ max_P_SM
            else:
                row["corrP/SM+P"] = 0
            
            return row
                                                         
        
        outputTable = {}
        
        products = list(cpTable.loc[cpTable["type"] == "Product", "name"])
        internalSTD = list(cpTable.loc[cpTable["type"] == "InternalSTD", "name"])
        SMs = list(cpTable.loc[cpTable["type"] == "Reactant", "name"])
        byproducts = list(cpTable.loc[cpTable["type"] == "Byproduct", "name"])
        #generate the structure of the table so that it is independant of the 
        #hits that are found

        for i in range(1, self.plate_row_no * self.plate_col_no + 1):
            well_id = plate_tools.getUserReadableWell(i, self.plate_col_no)
            total_area_abs = lcData.loc[lcData["well"] == i][["areaAbs"]].sum().values[0]

            goingIn = {
                "well_no": i,
                "Well": well_id,
                "canonSMILES": "",
                "SMarea": 0,
                "SMareaAbs": 0,
                "Parea": 0,
                "PareaAbs": 0,
                "STDarea": 0,
                "STDareaAbs": 0,
                "Uarea": 100,
                "UareaAbs": total_area_abs,
                "corrSMarea": 0,
                "corrParea": 0,
                "corrSTDarea": 0,
                "P/SM+P": 0,
                "P/STD":0
                }
            
            
            for by_prod in byproducts:            
                name = f'{by_prod}area'
                corr_name = f'corr{by_prod}area'
                abs_name = f'{by_prod}areaAbs'
                name_std = f'{by_prod}/STD'

                goingIn[name] = 0
                goingIn[abs_name] = 0
                goingIn[corr_name] = 0
                goingIn[name_std] = 0   

            outputTable[i] = goingIn  

        #Fill in the output table using the given cpTable and lcData
        cpTable.apply(addCompoundToTable, axis = 1)
        
        #reset any non-sensical values that arise as a result of miniscule rounding errors
        #for the unidentified area (Uarea, UareaAbs)
        for i in range(1, self.plate_row_no * self.plate_col_no + 1):

            well = outputTable[i]
            
            if round(well["Uarea"]) <= 0 or round(well["UareaAbs"]) <= 0: 
                well["Uarea"] = 0
                well["UareaAbs"] = 0

        #generate pandas dataframe, and sort it on the well_no (index)
        self.df = pd.DataFrame.from_dict(outputTable, orient="index")
        self.df.sort_index(inplace=True)
        
        #calculate the hybrid metrics
        self.df["P/SM+P"] = self.df.apply(lambda row: round(row["PareaAbs"] / (row["SMareaAbs"] + row["PareaAbs"]), 2) 
                                             if row["PareaAbs"] != 0 else 0, axis = 1)
        
        self.df["P/STD"] = self.df.apply(lambda row: round(row["PareaAbs"] / row["STDareaAbs"], 2) 
                                             if row["STDareaAbs"] != 0 else 0, axis = 1)
        for byprod in byproducts:
            name_std = f'{by_prod}/STD'
            name = f'{by_prod}areaAbs'
            self.df[name_std] = self.df.apply(lambda row: round(row[name] / row["STDareaAbs"], 2) 
                                             if row["STDareaAbs"] != 0 else 0, axis = 1)
            
        #calculate the normalised metrics
        
        self.df = self.df.apply(calculateNormalisedMetrics, axis = 1)

        #Add the sample IDs to the dataframe
        self.df["SampleID"] = self.sample_IDs

    def genAnalyticalTable(self, platemapDF, cpTable, lcData, sample_IDs, products, SMs, by_products, standards, plate_col_no):
        """ 
        Generate an analytical table that includes all input and output information
        for a plate, including reagent amounts, IDs, plate temperature, irradiation
        power, etc. 
        There is a row for each piece of information, for each reagent, product or
        plate parameter, for each well in the plate. 
        The column headers are based on the assumption that the co-repository, "PyParse_designer",
        was used to generate the platemap. Other methods of generating the required platemap 
        may need to edit the column headers used below. 
        
        :param platemapDF: Pandas datatable containing the full platemap
        :param cpTable: Pandas datatable containing the full output by compound
        :param lcData: Pandas DataFrame containing LC peak data
        :param sample_IDs: Dictionary of all sample IDs, indexed by well number
        :param products/SMs/by_products/internalSTD: List of names of compounds, by each compound type

        :output: Pandas DataFrame added to the Output object. 
        """

        class analyteRow:
            def __init__(self, well, wellID, sample_ID, smiles, analyte_class, data_type, data_value, data_desc, data_unit):
                self.well = well
                self.well_ID = wellID
                self.Reaction_ID = sample_ID
                self.analyte_smiles = smiles
                self.analyte_class = analyte_class
                self.data_entry = data_value
                self.data_type = data_type
                self.data_units = data_unit
                self.data_description = data_desc
                
        aTable = []

        for index, compound in cpTable.iterrows():
            
            if compound["name"] in products:
                a_class = "desired product"
            elif compound["name"] in by_products:
                a_class = "byproduct"
            elif compound["name"] in SMs:
                a_class = "limiting reactant"
            elif compound["name"] in standards:
                a_class = "internalstd"
            else:
                a_class = "impurity"

            for well in compound["locations"]:
                well_name = plate_tools.getUserReadableWell(well, plate_col_no)
                valid_hit = [i for i in compound["final_result"]["green"] if well == lcData.at[i, "well"]]

                if len(valid_hit) > 0:
                    valid_hit = valid_hit[0]
                    MS_plus = "-"
                    MS_minus = "-"
                    
                    peakID = -1
                    for i in compound["final_result"]["green"]:
                        if lcData.at[i, "well"] == well:
                            peakID = lcData.at[i, "peakID"]
    
                    
                    if peakID > -1 and a_class != "impurity":
                        for hit in compound["hits"]:
                            if hit["well"] == well and hit["peakID"] == peakID:
                                MS_plus = hit["MS_plus"]
                                MS_minus = hit["MS_minus"]
                    
                    aTable.append(analyteRow(well, well_name, sample_IDs[well], index, a_class, "Retention Time", lcData.at[valid_hit, "time"], "Time at which compound elutes", "min").__dict__)
                    aTable.append(analyteRow(well, well_name, sample_IDs[well], index, a_class, "ES+", MS_plus, "Positive m/z", "Da/Q").__dict__)
                    aTable.append(analyteRow(well, well_name, sample_IDs[well], index, a_class, "ES-", MS_minus, "Negative m/z", "Da/Q").__dict__)
                    aTable.append(analyteRow(well, well_name, sample_IDs[well], index, a_class, "Percentage Area", lcData.at[valid_hit, "area"], "LCMS UV Peak percentage area", "%").__dict__)
                    aTable.append(analyteRow(well, well_name, sample_IDs[well], index, a_class, "Area Absolute", lcData.at[valid_hit, "areaAbs"], "LCMS UV Peak absolute area", "AU").__dict__)
                #If no hit was found, capture in the table that the product was expected but not observed. 
                else: 
                    aTable.append(analyteRow(well, well_name, sample_IDs[well], index, a_class, "Retention Time", 0, "Time at which compound elutes", "min").__dict__)
                    aTable.append(analyteRow(well, well_name, sample_IDs[well], index, a_class, "ES+", 0, "Positive m/z", "Da/Q").__dict__)
                    aTable.append(analyteRow(well, well_name, sample_IDs[well], index, a_class, "ES-", 0, "Negative m/z", "Da/Q").__dict__)
                    aTable.append(analyteRow(well, well_name, sample_IDs[well], index, a_class, "Percentage Area", 0, "LCMS UV Peak percentage area", "%").__dict__)
                    aTable.append(analyteRow(well, well_name, sample_IDs[well], index, a_class, "Area Absolute", 0, "LCMS UV Peak absolute area", "AU").__dict__)
            
        #Add the remaining data from the platemapDF, containing info on IDs, amounts, quantity used, etc

        entry_types = ["desired product", "limiting reactant", "excess reactant", "internalstd", "base", "catalyst1", "catalyst2", "additive", "solvent1", "solvent2"]
        #Up to 10 columns each are allowed for byproducts and excess reactants. Note that more than 10 different compounds
        #are permitted in this category, as each column can allow multiple compounds. 
        for i in range(1, 11):
            entry_types.append(f'byproduct{i}')
            entry_types.append(f'excess reactant{i}')
        for _, row in platemapDF.iterrows():
            well_as_num = plate_tools.convertWellToNum(row["well"], plate_col_no)
            
            for entry in entry_types:
                if f'{entry} smiles' in row:
                    if row[entry+" smiles"] != "":
                        if entry == "solvent1" or entry == "solvent2":
                            amount_unit = "uL"
                        else:
                            amount_unit = "umol"

                        #Ensure that catalysts and co-catalysts are all given the same analyte class
                        if entry == "catalyst1" or entry == "catalyst2":
                            a_class = "catalyst"
                        elif "byproduct" in entry:
                            a_class = "byproduct"
                        elif "excess reactant" in entry:
                            a_class = "excess reactant"
                        else:
                            a_class = entry

                        #generate a canonical smiles for use as an index
                        try:
                            mol = Chem.MolFromSmiles(row[entry+" smiles"].strip())
                            c_smiles = Chem.MolToSmiles(mol)

                        except:
                            c_smiles = row[entry+" smiles"]
                        
                        if f'{entry} amount' in row:
                            if row[entry +" amount"] != "":
                                aTable.append(analyteRow(well_as_num, row["well"], sample_IDs[well_as_num], c_smiles, a_class, "Expected Amount", row[entry +" amount"], 
                                        "Amount of analyte expected in plate.", amount_unit).__dict__)

                        if f'{entry} id' in row:
                            if row[entry +" id"] != "":
                                aTable.append(analyteRow(well_as_num, row["well"], sample_IDs[well_as_num], c_smiles, a_class, "Analyte ID", row[entry +" id"], 
                                                "ID of analyte in well.", "N/A").__dict__)
                        if f'{entry} name' in row:
                            if row[entry +" name"] != "":
                                aTable.append(analyteRow(well_as_num, row["well"], sample_IDs[well_as_num], c_smiles, a_class, "Analyte Name", row[entry +" name"], 
                                                "Name of analyte in well.", "N/A").__dict__)

            if "time" in row:
                aTable.append(analyteRow(well_as_num, row["well"], sample_IDs[well_as_num], "Time", "plate parameter", "Time", row["time"], 
                                "Reaction time prior to analysis.", "h").__dict__)

            if "temperature" in row:
                aTable.append(analyteRow(well_as_num, row["well"], sample_IDs[well_as_num], "Temperature", "plate parameter", "Temperature", row["temperature"], 
                                "Reaction temperature prior to analysis", "C").__dict__)

            if "irradiation_power" in row:
                if(row["irradiation_power"] != ""):
                    aTable.append(analyteRow(well_as_num, row["well"], sample_IDs[well_as_num], "Irradiation Power", "plate parameter", "Irradiation Power", row["irradiation_power"], 
                                    "Irradiation power applied", "mW per well").__dict__)
            if "irradiation_wavelength" in row:
                if(row["irradiation_wavelength"] != ""):
                    aTable.append(analyteRow(well_as_num, row["well"], sample_IDs[well_as_num], "Irradiation Wavelength", "plate parameter", "Irradiation Wavelength", row["irradiation_wavelength"], 
                                    "Irradiation wavelength applied", "nm").__dict__)




        #Build the Pandas dataframe and sort it inplace
        self.aTable = pd.DataFrame(aTable)
        self.aTable.sort_values(["well", "analyte_smiles"], inplace = True)
        del self.aTable["well"]


