import math
import pandas as pd
import plate_tools

class Output:
    def __init__(self, sample_IDs, plate_col_no, plate_row_no):
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

        :param compoundDF: The pandas datatable containing all compounds with their respective hits.
        :param internalSTD: the name of the internalSTD
        :param SMs: a list of indices for the starting materials
        :param products: a list of indices for the products
        :param by_products: a list of indices for the by-products
        :param total_area_abs: A float corresponding to the sum of all peak_area_absolutes

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
            