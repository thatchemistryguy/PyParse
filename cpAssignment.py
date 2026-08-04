
#External Libraries
import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
import logging
from statistics import mean
import math

#Internal Libraries
import plate_tools

class Assignment:

    """
    Processes the platemap defined by the user and creates a Pandas DataFrame with a row for each compound
    along with the wells it is expected in. 
    
    This object is designed to work in partnership and parallel to an Object for the LCMS data, e.g. getWatersData.rawWatersData

    In this class are functions used to find hits by comparison of the compound with the data in the LCMS data Object,
    refine them, then store a "final_result" as a new column. This final_result is a list of all the peaks that have been 
    assigned to the compound in that row, as the index for the corresponding row in the lcData DataFrame, which is provided
    post-initialisation (e.g. rawWatersData.rawDADTable). 

    The following are descriptions of the Pandas dataframes that will be generated for this object:

    self.cpTable = pd.DataFrame(columns=['canonSMILES', 'type', 'locations', 'name', 'rt', 'comments', 'mass1', 'mass2', 'mass3', 'hits', 
    'clusters', 'culster_bands', 'refined_clusters', 'discarded_clusters', 'final_result'])
    """

    def __init__(self, filename, plate_col_no, plate_row_no):   

        """
        Initialises the Assignment object. Requires the platemap to be provided, as well as the number of 
        rows+cols in the plate. 

        :param filename: a string describing the file location of the .csv platemap. 
        :param plate_row_no: an integer describing the number of rows in the plate, as specified by the user
        :param plate_col_no: an integer describing the number of columns in the plate, as specified by the user
        
        :return: A python object, containing a DataFrame created from the .csv platemap. 
        """

        #read csv file into dataframe
        #replace empty cells with an empty string
        #convert all column names to lower case and remove whitespace
        self.inputCSV = pd.read_csv(filename)
        self.inputCSV.fillna("", inplace=True)
        self.inputCSV.columns = self.inputCSV.columns.str.strip().str.lower()
        
        self.plate_col_no = plate_col_no
        self.plate_row_no = plate_row_no
        
    
    def convertWellToNum(self, wells):
        """
        Function to convert a list human-readable well like B5 into machine format like 
        the number 11. 

        :param wells: A list of strings, each describing one well. 

        :return: A list of integers, where each integer corresponds to a well number. 
        """
        result = []
        
        for well in wells:
            row = well[0]
            column = well[1:]
            #if the format of the row/column conforms to expectations
            if isinstance(int(column[0]), int):
                result.append(int((ord(row) - 65) * self.plate_col_no + int(column)))
            #Plates with more than 26 rows are unsupported at present. 
            else:
                logging.info("The well specified implies an unsupported plate.")
                sys.exit(2)
        return result
    
    
    def getCanonSmiles(self, smiles):
        """
        Fn to convert smiles into canonicalised smiles

        :param smiles: A string, representing a SMILES molecule
        :return: A string, representing a canonicalised SMILES. 
        """

        mol = Chem.MolFromSmiles(smiles.strip())
        return Chem.MolToSmiles(mol)
            
    def generateCPTable(self):
        """
        Takes the self.inputCSV dataframe (the platemap the user provided) and 
        performs an unpivot operation with additional santisation to create a Pandas
        DataFrame where each compound is it's own row with a list of expected locations
        and expected mass ions. 

        :return: DataFrame of compounds. 
        self.cpTable = pd.DataFrame(columns=['canonSMILES', 'type', 'locations', 'name', 'rt', 'comments', 'mass1', 'mass2', 'mass3', 'hits', 
        'clusters', 'culster_bands', 'refined_clusters', 'discarded_clusters', ])
        """

        
        compound_list = []

        #The following describe the standardised column names the user should have used in the 
        #platemap, with the associated compound type next to them. 
        type_dic = {
            "desired product smiles": "Product",
            "limiting reactant smiles":"Reactant",
            "internalstd smiles": "InternalSTD"
        }

        #In the absence of the user providing a name for the compound, this counter allows
        #auto-naming using incremental values. 
        counter = {
            "desired product smiles": 1,
            "limiting reactant smiles": 1,
            "byproduct": 1
        }
        
        #Canonicalise all the incoming smiles in key columns in case the user hadn't done so already
        for col in self.inputCSV:
            if col in type_dic or ("byproduct" in col and "smiles" in col):
                self.inputCSV[col] = self.inputCSV[col].apply(self.getCanonSmiles)
        
        #For each column in the platemap, extract information based on standardised column names. 
        for col in self.inputCSV:
            if col in type_dic or ("byproduct" in col and "smiles" in col):
                cpname_column = f'{col.split(" smiles")[0]} name'
                cprt_column = f'{col.split(" smiles")[0]} rt'
                
                #group the input CSV by canonical smiles in that column, aggregated the wells into a list
                groupeddf = self.inputCSV.groupby(col, as_index=False)[["well"]].agg(lambda x: list(x))
                
                #Iterate through each of those grouped entries
                for index, row in groupeddf.iterrows():
                    cpindex = row[col]
                    cptype = type_dic[col] if col in type_dic else "Byproduct"
                    #get a name for the compound if one was provided
                    name = ""
                    if cpname_column in self.inputCSV.columns:
                        name_series = self.inputCSV.loc[self.inputCSV[col] == row[col]][cpname_column]
                        potential_names = [x for x in name_series if x != ""]
                        if len(potential_names) != 0:
                            name = potential_names[0]
                    #if a name could not be generated, create a generic one using a simple counter
                    if name == "":
                        if col == "internalstd smiles":
                            name = "InternalSTD"
                        elif col in counter:
                            name = f'{type_dic[col]}{counter[col]}'
                            counter[col] = counter[col] + 1
                        elif "byproduct" in col:
                            name = f'Byproduct{counter["byproduct"]}'
                            counter["byproduct"] = counter["byproduct"] + 1
                    
                    #get a rentention time for the compound if one was provided
                    rt = 0
                    if cprt_column in self.inputCSV.columns:
                        rt_series = self.inputCSV.loc[self.inputCSV[col] == row[col]][cprt_column]
                        potential_rt = [x for x in rt_series if x != ""]
                        if len(potential_rt) != 0:
                            rt = potential_rt[0]

                    if cpindex != "":        
                        new_entry = {
                            "canonSMILES": cpindex,
                            "type": cptype,
                            "locations": row["well"],
                            "name": name,
                            "rt": rt,
                            "comments": []
                        }
                        compound_list.append(new_entry)
        
        #Create DataFrame from the compound list. 
        self.cpTable = pd.DataFrame(compound_list)
        
        #set the index to be the canonicalised smiles
        self.cpTable.index = list(self.cpTable["canonSMILES"])
        
        #convert the "A1" style well IDs in the "locations" column into a integer, 
        #to allow matching to a well in the LCMS rawData (e.g. getWatersData.rawWatersData)
        self.cpTable["locations"] = self.cpTable["locations"].apply(self.convertWellToNum)
    

    
    def generateEMs(self, calc_boc):

        """
        Generate Exact Molecular Weights for each compound in the cpTable DataFrame. 
        May generate up to three different masses that will be used for searching for hits, 
        accounting for halide isotopes and Boc protecting group degradation patterns as they
        fly through the mass spec. 
        
        :param calc_boc: string, either "True" or "False"

        :return: Updated cpTable with the new columns "mass1", "mass2" and "mass3". 
        """
        
        def getMW(smiles):
            mol = Chem.MolFromSmiles(smiles)
            return round(Descriptors.ExactMolWt(mol), 2)
        
        def transform_and_getMW(smiles, smirks, stage):
            try:
                mol = Chem.MolFromSmiles(smiles)
                rxn1 = AllChem.ReactionFromSmarts(smirks)
                new_mol1 = rxn1.RunReactants((mol, ))[0][0]
                #Sanitise the molecule to make sure that a sensible molecule was produced. 
                Chem.SanitizeMol(new_mol1)
                return round(Descriptors.ExactMolWt(new_mol1), 2)
            except:
                if ("Cl" in smiles or "Br" in smiles) and stage == "mass2":
                    mol = Chem.MolFromSmiles(smiles)
                    return round(Descriptors.ExactMolWt(mol)+2, 2)
                else:
                    return 0
            
        self.cpTable["mass1"] = self.cpTable["canonSMILES"].apply(lambda smiles: getMW(smiles))
        
        if calc_boc == "True":
            
            smirks1 = "[NX3,n:1][C:2](=[O:3])[O:4][C]([CH3])([CH3])[CH3]>>[*:1][C:2](=[O:3])[O:4]"
            self.cpTable["mass2"] = self.cpTable["canonSMILES"].apply(lambda smiles: transform_and_getMW(smiles, smirks1, "mass2"))
            
            smirks2 = "[NX3,n:1][C](=[O])[O][C]([CH3])([CH3])[CH3]>>[*:1][H]"
            self.cpTable["mass3"] = self.cpTable["canonSMILES"].apply(lambda smiles: transform_and_getMW(smiles, smirks2, "mass3"))
            
    def clusterHits(self, hits, lcData, time_abs_tol = 0.025):
        
        df = lcData.loc[lcData.index.isin(hits)].sort_values("time")

        clusters = []
        for index in df.index:
            if len(clusters) == 0:
                clusters.append([index])
            else:
                clusterFound = False
                for cluster in clusters:
                    mean_rt = mean([df.at[i, "time"] for i in cluster])

                    if abs(mean_rt - df.at[index, "time"]) < time_abs_tol:
                        cluster.append(index)
                        clusterFound = True
                        break
                if not clusterFound:
                    clusters.append([index])

        clusterbands = []
        for cluster in clusters:
            mean_time = lcData.loc[lcData.index.isin(cluster)]["time"].mean()
            clusterbands.append([round(mean_time, 5), cluster])
            
        return [clusters, clusterbands]
        
    def findHits(self, msData, lcData, mass_abs_tol = 0.5, min_massconf_threshold = 10, 
                                    time_abs_tol = 0.025, calc_higherions = "True"):
        """
        This function operates on the full cpTable finding hits for each compound (row)
        in said table. It does this by either finding peaks in the lcData DataFrame that have 
        a retention time close to that specified by the user in the platemap (optional) or
        by finding m/z values in the msData that are close to the exact mass expected. 

        :param msData: Pandas DataFrame, containing mass spec data
        :param lcData: Pandas DataFrame, containing LC peak data
        :param mass_abs_tol: float
        :param min_massconf_threshold: float
        :param time_abs_tol: float
        :param calc_higherions: string, either "True" or "False"

        :return: An updated cpTable with a new column, "hits"
        """

        def findHitsForRow(row):

            def getMassWithHighestConf(mass_row, df):
                """
                From all the matching MS+ and MS- ions, fetch the one with the largest intensity.    
                """
                MS_plus = list(df.loc[(df.index.isin(MS_hits)) & (df["MStype"] == "+") &
                                        (df["MSintensity"] == mass_row["max_intensity"]), "MSvalue"])
                MS_minus = list(df.loc[(df.index.isin(MS_hits)) & (df["MStype"] == "-") &
                                        (df["MSintensity"] == mass_row["max_intensity"]), "MSvalue"])
                mass_row["MS_plus"] = MS_plus[0] if len(MS_plus) > 0 else "-"
                mass_row["MS_minus"] = MS_minus[0] if len(MS_minus) > 0 else "-"
                
                return mass_row
        
            #Filter the msData DataFrame to only those locations we're interested in
            #for this compound. 
            df = msData.loc[msData["well"].isin(row["locations"])]
            
            MS_hits = []
            
            for mass in [row["mass1"], row["mass2"], row["mass3"]]:
                if calc_higherions == "True":
                    hits = df.loc[(df["MStype"] == "+") & 
                                    ((abs(df["MSvalue"] - (mass + 1.01)) <= mass_abs_tol) |
                                        (abs(df["MSvalue"] - (mass + 2.02)/2) <= mass_abs_tol) |
                                        (abs(df["MSvalue"] - (mass + 3.03)/3) <= mass_abs_tol))]
                else:
                    hits = df.loc[(df["MStype"] == "+") & 
                                    (abs(df["MSvalue"] - (mass + 1.01)) <= mass_abs_tol)]
                
                MS_hits = MS_hits + list(hits.index.values)
                
                hits = df.loc[(df["MStype"] == "-") & 
                                    (abs(df["MSvalue"] - (mass - 1.01)) <= mass_abs_tol)]
                MS_hits = MS_hits + list(hits.index.values)
            
            grouped = (df[df.index.isin(MS_hits)].groupby(["well", "peakID"], as_index = False)
                        .agg(mass_conf = ("perc_intensity", "sum"), max_intensity = ("MSintensity", "max")))
            grouped = grouped.apply(getMassWithHighestConf, args=(df,), axis = 1)
            
            grouped = grouped.loc[grouped["mass_conf"] >= min_massconf_threshold]
            grouped = grouped.drop("max_intensity", axis=1)
            
            #fetch any peaks where the time matches the retention time provided
            
            new_slice = lcData.loc[lcData["well"].isin(row["locations"]) & 
                                    (abs(lcData["time"] - row["rt"]) < time_abs_tol), ["well", "peakID"]]
            #For these additional hits, set the requisite columns so that the two sets of hits can be merged into a 
            #single dataframe. 
            new_slice["mass_conf"] = 0
            new_slice["MS_plus"] = "-"
            new_slice["MS_minus"] = "-"

            if len(new_slice.index) != 0:
                row["hits"] = pd.concat([grouped, new_slice], axis=0).to_dict("records")
            else:
                row["hits"] = grouped.to_dict("records")                     

            return row

        #Apply the above function to all rows in the cpTable. 
        self.cpTable = self.cpTable.apply(findHitsForRow, axis=1)
        
        
    def validateHits(self, lcData, uvData, time_abs_tol = 0.025, massconf_threshold = 0.5, uv_abs_tol = 10,
                    uv_cluster_threshold = 0.5, uv_match_threshold = 0.5, cluster_size_threshold = 0.8, min_no_of_wells = 5,
                    validate = "True", mass_or_area = "mass_conf"):
        """
        This function operates on the full cpTable validating the hits for each compound 
        and ultimately delivering a new column called "final_result". Other columns added 
        are "clusters", "refined_clusters", "discarded_clusters", "comments" and "clusters_indexed_by_well".
        The "final_result" column contains the hits that should be used to build the output heatmap. 

        :param lcData: Pandas DataFrame, containing LC peak data
        :param uvData: Pandas DataFrame, containing uv absorbance data
        :param time_abs_tol: float
        :param min_massconf_threshold: float
        :param uv_abs_tol: integer
        :param uv_cluster_threshold: float
        :param uv_match_threshold: float
        :param cluster_size_threshold: float, between 0 and 1
        :param min_no_of_wells: integer, >1
        :param validate: string, either "True" or "False"
        :param mass_or_area: string, either "mass_conf" or "area"


        :return: An updated cpTable with several new columns, most importantly "final_result".
        """
              
        def clusterRow(row):
            """
            Groups the hits so that ones with similar retention times are in the same group. 

            :param row: A row from a Pandas DataFrame, that must contain columns called "hits"
                                    and "comments"
                                    
            :return: The amended row, where the hits have been grouped and these groups placed into a 
                    new column called "clusters". A second column called "cluster_bands" is also added which
                    details the mean retention time for each cluster. 
            """    
            hits = []
            for hit in row["hits"]:
                hits = (hits + list(lcData[(lcData["well"] == hit["well"]) & 
                                                                (lcData["peakID"] == hit["peakID"])].index.values)) 

            cluster_result = self.clusterHits(hits, lcData, time_abs_tol)

            row["clusters"] = cluster_result[0]
            row["cluster_bands"] = cluster_result[1]

            return row
                       
        def selectCluster_ifrt(row):

            #If the user has specified a retention time, we should select only the cluster 
            #that is closest to that retention time, and within the specified time_abs_tol
            if row["rt"] != 0:
                suitable_clusters = [index for index, i in enumerate(row["cluster_bands"]) if abs(i - row["rt"]) < time_abs_tol]

                #If there is more than one cluster close to the specified retention time
                #take only the cluster which is closest 
                if len(suitable_clusters) > 1:
                    diffs = [row["cluster_bands"][index]-row["rt"] for index in suitable_clusters]
                    index_min = min(range(len(diffs)), key=diffs.__getitem__)
                    row["clusters"] = [row["clusters"][suitable_clusters[index_min]]]
                    #Update the cluster bands to only include the correct label
                    cluster_bands = [row["cluster_bands"][suitable_clusters[index_min]]]

                    row["comments"].append("<strong>Multiple clusters were found close the specified"
                                        " retention time.</strong>")
                    row["comments"].append(f'<strong>Cluster {index_min} was selected as it was closest'
                                        ' to the specified retention time.</strong>')
                elif len(suitable_clusters) == 1:
                    row["clusters"] = [row["clusters"][suitable_clusters[0]]]
                    row["cluster_bands"] = [row["cluster_bands"][suitable_clusters[0]]]
                    row["comments"].append("<strong>A single cluster was found close the specified"
                                        " retention time and this was selected for analysis.</strong>")
                else:
                    row["comments"].append("<strong>No cluster was found near to the specified "
                                        "retention time. Proceeding with analysis using all "
                                        f'{len(row["clusters"])} clusters.</strong>')
            return row
            
        def refineClustersByTime(row):
            """
            Takes in input cluster of all the hit peaks, 
            and refines them by finding a mid-value for the retention
            time based on which hit has the greatest number of nearest neighbours. 
            Sorts the best hits into "green", uncertain ones into "orange" 
            and those where another peak closer to the mid-value was found
            in the same well into "discarded". 

            :param row: A row from a Pandas DataFrame, that must contain a column called "clusters"

            :return: The amended row, where the clusters have been refined into green, orange and discarded categories. 
            """

            [clusters, comments, expected_rt] = [row["clusters"], row["comments"], row["rt"]]
            refined_clusters = []

            for cluster in clusters:
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
                        mid_value = lcData.at[i, "time"]
                        mid_values.append([mid_value, len([lcData.at[j, "time"] for j in cluster 
                                                           if abs(lcData.at[j, "time"] - mid_value) < time_abs_tol/4])])
                    mid_value = max(mid_values, key = lambda x: x[1])[0]

                #sort the peaks by the well they occupy
                peaks_by_wells = {}
                for i in cluster:
                    if lcData.at[i, "well"] not in peaks_by_wells:
                        peaks_by_wells[lcData.at[i, "well"]] = []
                    peaks_by_wells[lcData.at[i, "well"]].append(i)

                #For each well, select the peak that's closest to the mid-value in cases
                #where there was more than one hit in that cluster in one well
                #Use peakAdded to ensure that only a single peak per well is added to green, in the 
                #unlikely case that there are two peaks in the same well with the same retention time
                #(i.e. LCMS machine processing error)
                for index, well in peaks_by_wells.items():
                    if len(well) > 1:
                        min_diff = min([abs(lcData.at[i, "time"]-mid_value) for i in well])
                        peakAdded = False
                        for i in well:

                            if abs(lcData.at[i, "time"]-mid_value) == min_diff and not peakAdded:
                                refined_cluster["green"].append(i)
                                peakAdded = True
                            else:
                                refined_cluster["discarded"].append(i)
                                row["comments"].append(f'Peak at {lcData.at[i, "time"]} '
                                            f'in well {plate_tools.getUserReadableWell(index, self.plate_col_no)} was discarded '
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
                
                for i in refined_cluster["green"]:
                    if abs(lcData.at[i, "time"] - mid_value) < time_abs_tol / 2:
                        ref2_cluster["green"].append(i)
                    else:
                        ref2_cluster["orange"].append(i)
                        comments.append(f'Peak at {lcData.at[i, "time"]} in '
                                       f'well {plate_tools.getUserReadableWell(lcData.at[i, "well"], self.plate_col_no)} was '
                                       'marked as tentative as it was found to be too '
                                       'far from the mid-value of the cluster.')
                refined_clusters.append(ref2_cluster)

            row["clusters"] = refined_clusters
            return row
        
        def refineClustersByMassConf(row):
            """
            Takes in input cluster of all the hit peaks, 
            and refines them by ensuring all peaks have a similar mass confidence
            to the cluster's mean. Those which do are left in "green"; 
            those which don't are moved to the "orange" category.
            The "mass confidence" is just the relative intensity of the ion for that eluting peak. 

            :param row: A row from a Pandas DataFrame, that must contain columns called "clusters" and "comments"

            :return: The amended row, where the clusters have been refined and comments added to the "comments" column
            """
            refined_clusters = []
            for cluster in row["clusters"]:
                refined_cluster = {
                    "green": [],
                    "orange": cluster["orange"],
                    "discarded": cluster["discarded"],
                }
                if len(cluster["green"]) > 0:
                    #Get a dictionary of the mass_conf for each peak, indexed by the index present in lcData
                    mass_conf_dict = {}
                    for i in cluster["green"]:
                        well = lcData.at[i, "well"]
                        peakID = lcData.at[i, "peakID"]
                        for j in row["hits"]:
                            if j["well"] == well and j["peakID"] == peakID:
                                mass_conf_dict[i] = j["mass_conf"]

                    #find what the mean mass confidence is of all peaks
                    #currently under the "green" category
                    total = sum(mass_conf_dict[x] for x in cluster["green"])
                    mean_mass_conf = total / len(cluster["green"])

                    for i in cluster["green"]:
                        if mass_conf_dict[i] < mean_mass_conf * massconf_threshold:
                            refined_cluster["orange"].append(i)
                            row["comments"].append(f'Peak at {lcData.at[i, "time"]} in well '
                                            f'{plate_tools.getUserReadableWell(lcData.at[i, "well"], self.plate_col_no)} '
                                            'was marked tentative due to poor mass '
                                            'confidence in comparison to the rest of the cluster.')
                        else:
                            refined_cluster["green"].append(i)
                refined_clusters.append(refined_cluster)
            row["clusters"] = refined_clusters
            return row
        
        def refineClustersByUV(row):
            """
            Takes in input cluster of all the hit peaks, 
            and refines the cluster by ensuring all peaks have a similar set
            of UV maxima. Those which do are left in "green", those which
            don't are moved to the "orange" category

            :param row: A row from a Pandas DataFrame, that must contain columns called "clusters" and "comments"
            
            :return: The amended row, where the clusters have been refined and comments added to the "comments" column
            """
            refined_clusters = []
            
            for cluster in row["clusters"]:
                refined_cluster = {
                    "green": [],
                    "orange": cluster["orange"],
                    "discarded": cluster["discarded"],
                }

                if len(cluster["green"]) > 0:
                    

                    UVclusters = []

                    for i in cluster["green"]:
                        for UV in uvData.loc[(uvData["well"] == lcData.at[i, "well"]) & (uvData["peakID"] == lcData.at[i, "peakID"])]["UVvalue"]:
                            if len(UVclusters) == 0:
                                UVclusters.append([UV])

                            else:
                                clusterFound = False
                                for UVcluster in UVclusters:
                                    if abs(UVcluster[-1] - UV) < uv_abs_tol:
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
                            if len(UVcluster) >= lengthOfMostCommon * uv_cluster_threshold:
                                meanUV.append(mean(UVcluster))

                        for i in cluster["green"]:
                            UVvalues = (uvData.loc[(uvData["well"] == lcData.at[i, "well"]) 
                                               & (uvData["peakID"] == lcData.at[i, "peakID"])]["UVvalue"].values)
                            intersection = []
                            for meanvalue in meanUV:
                                for UV in UVvalues:
                                    if abs(meanvalue - UV) < uv_abs_tol:
                                        intersection.append(meanvalue)

                            if len(intersection) >= len(meanUV) * uv_match_threshold:
                                refined_cluster["green"].append(i)
                            else:
                                refined_cluster["orange"].append(i)
                                row["comments"].append(f'Peak at {lcData.at[i, "time"]} in well '
                                                f'{plate_tools.getUserReadableWell(lcData.at[i, "time"], self.plate_col_no)} '
                                                'was marked tentative due to mismatch in UV maxima with the rest of the cluster.')
                    #If UVclusters is an empty set, pass through the original hits without further
                    #validation. 
                    else:
                        row["comments"].append('UV validation not permitted where cluster contains peaks with no UV data. '
                            'UV validation was not performed.')
                        refined_cluster["green"] = cluster["green"]
                else:
                    row["comments"].append('UV data was not found for the plate, so UV validation was not performed.')
                    refined_cluster["green"] = cluster["green"]
                refined_clusters.append(refined_cluster)
            
            row["clusters"] = refined_clusters
            return row
        
        def selectClusterByMassConf(row):
            """
            If more than one cluster was found for the compound, 
            this function is called to try to select a single cluster based on 
            which cluster has the highest mean massConf. If more than one cluster
            has a close-to-highest-mean massconf, take them all. 

            :param row: A row from a Pandas DataFrame, that must contain columns called "clusters" and "comments"
            
            :return: The amended row, where selected clusters have been added to a new column called "refined_clusters"  and
                    deselected clusters added to a new column called "discarded_clusters". Additional comments also added. 
            """
            
            refined_clusters = []
            discarded_clusters = []

            means = []
            #find the mean massconf for each cluster
            for cluster in row["clusters"]:
                if len(cluster["green"]) > 0:
                    #Get a dictionary of the mass_conf for each peak, indexed by the index present in lcData
                    mass_conf_dict = {}
                    for i in cluster["green"]:
                        well = lcData.at[i, "well"]
                        peakID = lcData.at[i, "peakID"]
                        for j in row["hits"]:
                            if j["well"] == well and j["peakID"] == peakID:
                                mass_conf_dict[i] = j["mass_conf"]
                    
                    #Calculate the mean mass_conf of the cluster
                    mean_mass_conf = sum([mass_conf_dict[i] for i in cluster["green"]]) / len(cluster["green"])
                    means.append(mean_mass_conf)
                else:
                    means.append(0)
                    
            #find the maximum mean value to compare all clusters against
            max_mean = max(means)
            
            #Filter the clusters by those which have a mean mass confidence that is at least the specified
            #percentage of the maximum observed mean mass confidence (set by massconf_threshold, typically 80%)
            for i in range(len(means)):
                if max_mean * massconf_threshold < means[i] or (max_mean == 0 and means[i] == 0):
                    refined_clusters.append(row["clusters"][i])
                else:
                    discarded_clusters.append(row["clusters"][i])
                    
            row["refined_clusters"] = refined_clusters
            row["discarded_clusters"] = discarded_clusters
            row["comments"].append(f'{len(discarded_clusters)} clusters were discarded after refinement'
                                 ' by mass confidence. ')
            return row

        def selectClusterBySize(row):
            """
            If more than one cluster was found for the compound, 
            this function is called to try to select a single cluster based on 
            which cluster is the largest. If more than one cluster
            has a close-to-largest size, take them all. 

            :param row: A row from a Pandas DataFrame, that must contain columns called "refined_clusters", "discarded_clusters"
                        and "comments"
                        
            :return: The amended row, where selected clusters have been added to a new column called "refined_clusters"  and
                    deselected clusters added to a new column called "discarded_clusters". Additional comments also added. 
            """

            refined_clusters2 = []
            discarded_clusters = row["discarded_clusters"]

            lengths = [len(cluster["green"]) + len(cluster["orange"]) for cluster in row["refined_clusters"]]
            max_length = max(lengths)

            for cluster in row["refined_clusters"]:
                if len(cluster["green"]) + len(cluster["orange"]) > max_length * cluster_size_threshold:
                    refined_clusters2.append(cluster)
                else:
                    discarded_clusters.append(cluster)

            #If there is still more than one cluster, reset the process using only those peaks
            #that haven't been marked as suspicious (orange)
            if len(refined_clusters2) > 1:

                lengths = [len(cluster["green"]) for cluster in row["refined_clusters"]]
                max_length = max(lengths)

                #As long as at least one cluster had a green hit, filter the clusters by size
                if max_length != 0:
                    refined_clusters2 = []
                    discarded_clusters = row["discarded_clusters"]

                    for cluster in row["refined_clusters"]:
                        if len(cluster["green"]) > max_length * cluster_size_threshold:
                            refined_clusters2.append(cluster)
                        else:
                            discarded_clusters.append(cluster)
            
            row["comments"].append(f'{len(row["refined_clusters"]) - len(refined_clusters2)} clusters were discarded after refinement'
                                 ' by size comparison. ')

            if len(refined_clusters2) > 1:
                row["comments"].append(f'<strong>High Priority: More than one cluster was found for {row["name"]}. ' 
                                            'User MUST perform their own analysis.</strong>')

            row["refined_clusters"] = refined_clusters2
            row["discarded_clusters"] = discarded_clusters
        
            return row
        
        
            
        def refine_and_select(row): 
            """
            This function runs all the validation functions (refineClustersByTime, refineClustersByMassConf, 
            refineClustersByUV, selectClusterByMassConf, selectClusterBySize). It will do so if there are 
            sufficient wells to perform refinement. If not, this function reverts to the simplest strategy of 
            assuming all hits are "green". 

            :param row: A row from a Pandas DataFrame, that must contain columns called "hits", "clusters" and "comments"
                                    
            :return: The amended row, where hits have been refined, selected clusters have been added to a new column called "refined_clusters"  and
                     deselected clusters added to a new column called "discarded_clusters". Additional comments also added. 
            """
            if len(row["hits"]) > min_no_of_wells and validate == "True":
                row = refineClustersByTime(row)
                row = refineClustersByMassConf(row)
                row = refineClustersByUV(row)
                row = selectClusterByMassConf(row)
                row = selectClusterBySize(row)
            else:
                refined_clusters = []
                for cluster in row["clusters"]:
                    refined_cluster = {
                        "green": [],
                        "orange": [], 
                        "discarded": []
                    }
                    for i in cluster:
                        refined_cluster["green"].append(i)
                    refined_clusters.append(refined_cluster)
                row["refined_clusters"] = refined_clusters
                row["discarded_clusters"] = []
                
                if validate != "True":
                    row["comments"].append("Validation was not performed as requested by the user.")
                else:
                    row["comments"].append(f'Validation was not performed for {row["name"]} as '
                                        'there were insufficient hits.')
            return row

        def indexClusterByWells(row):
            """
            This function sorts all the hits from the refined clusters, whether they're
            green, orange, or discarded into lists that are indexed by the well those hits
            come from. This allows easier decision analysis by the "finalResult" function to
            select just a single hit for each well. 

            :param row: A row from a Pandas DataFrame, that must contain columns called "refined_clusters"

            :return: The amended row, where hits in the refined_clusters have been reindexed and placed into
                    a new column called "clusters_indexed_by_well". 
            """
            cluster_by_well = {}
            for cluster in row["refined_clusters"]:
                for i in cluster["green"]:
                    if lcData.at[i, "well"] not in cluster_by_well:
                        cluster_by_well[lcData.at[i, "well"]] = {
                                "green": [],
                                "orange": [],
                                "discarded": [],
                                }

                    cluster_by_well[lcData.at[i, "well"]]["green"].append(i)

                for i in cluster["orange"]:
                    if lcData.at[i, "well"] not in cluster_by_well:
                        cluster_by_well[lcData.at[i, "well"]] = {
                                "green": [],
                                "orange": [],
                                "discarded": [],
                                }

                    cluster_by_well[lcData.at[i, "well"]]["orange"].append(i)

                for i in cluster["discarded"]:
                    if lcData.at[i, "well"] not in cluster_by_well:
                        cluster_by_well[lcData.at[i, "well"]] = {
                                "green": [],
                                "orange": [],
                                "discarded": [],
                                }

                    cluster_by_well[lcData.at[i, "well"]]["discarded"].append(i)
                    
            row["clusters_indexed_by_well"] = cluster_by_well
            return row
        
        def finalResult(row):

            """
            This final function in the validateHits sueprfunction sorts through the hits that have been clusters,
            refined and reindexed by well, and gives an output that can be used ready for creating heatmaps and 
            other visualisations. 
            The "final_result" column is simply a list of all the LC peaks that should be used to create visualisations, 
            provided as just the DataFrame index for that LC peak in the lcData DataFrame that was provided. 

            :param row: A Pandas DataFrame row that contains a column called "clusters_indexed_by_well"
            :return: The amended row where a new column called "final_result" has been added. Comments also added to 
                    "comments" column. 

            """
            final_result = {
                "green": [],
                "discarded": []    
            }
            
            for well in row["clusters_indexed_by_well"].values():
                #Move all discarded peaks immediately to the analogous result in "final result"
                final_result["discarded"] += well["discarded"]
                
                #Depending on how many green or orange peaks were found in this well, proceed accordingly, 
                #where green peaks are selected in preference to orange ones. 
                if len(well["green"]) > 0:
                    if len(well["green"]) == 1:
                        final_result["green"].append(well["green"][0])
                    else:
                        
                        #Get a list of peaks that are deemed suitable, either those that have the largest mass confidence
                        #or, in the case where the peaks are selected purely by their area, a dummy list. 
                        #We're going 
                        if mass_or_area == "mass_conf":
                            #Get a dictionary of the mass_conf for each peak, indexed by the index present in lcData
                            mass_conf_dict = {}
                            for i in well["green"]:
                                peak_well = lcData.at[i, "well"]
                                peakID = lcData.at[i, "peakID"]
                                for j in row["hits"]:
                                    if j["well"] == peak_well and j["peakID"] == peakID:
                                        mass_conf_dict[i] = j["mass_conf"]
                                        
                            max_mass_conf = max(mass_conf_dict.values())
                                                               
                            suitable_indexes = [i for i in well["green"] if mass_conf_dict[i] == max_mass_conf]
                        else:
                            suitable_indexes = well["green"]
                            
                        sorted_list = sorted(well["green"], key = lambda x: lcData.at[x, "area"], reverse = True)

                        peakAdded = False
                        for i in sorted_list:
                            if i in suitable_indexes and peakAdded == False:
                                final_result["green"].append(i)
                                row["comments"].append(f'Peak with largest {mass_or_area} selected in '
                                                       'preference to others available for well '
                                                       f'{plate_tools.getUserReadableWell(lcData.at[i, "well"], self.plate_col_no)}.')
                                peakAdded = True
                            else:
                                final_result["discarded"].append(i)
                                row["comments"].append(f'The peak at {lcData.at[i, "time"]} in well '
                                                       f' {plate_tools.getUserReadableWell(lcData.at[i, "well"], self.plate_col_no)} '
                                                       f'was discarded as it had a smaller {mass_or_area} than an '
                                                       'otherwise equally likely alternative.')
                                                            
                    #If we had a green peak to select, we can safely move all the orange peaks into the discarded pile
                    if len(well["orange"]) > 0:
                        final_result["discarded"] += well["orange"]
                        for i in well["orange"]:
                            row["comments"].append(f'The peak at {lcData.at[i, "time"]} in well '
                                        f'{plate_tools.getUserReadableWell(lcData.at[i, "well"], self.plate_col_no)} '
                                        f'for {row["name"]} was discarded because a better match was found.')
                
                elif len(well["orange"]) == 1:
                    final_result["green"].append(well["orange"][0])
                    row["comments"].append(f'<strong>The tentative peak at {lcData.at[well["orange"][0], "time"]} in well '
                                        f'{plate_tools.getUserReadableWell(lcData.at[well["orange"][0], "well"], self.plate_col_no)} was '
                                        f'used as there was no better option. User should check this well.</strong>')
                elif len(well["orange"]) > 1:
                    #Get a list of peaks that are deemed suitable, either those that have the largest mass confidence
                    #or, in the case where the peaks are selected purely by their area, a dummy list. 
                    #We're going 
                    if mass_or_area == "mass_conf":
                        #Get a dictionary of the mass_conf for each peak, indexed by the index present in lcData
                        mass_conf_dict = {}
                        for i in well["orange"]:
                            peak_well = lcData.at[i, "well"]
                            peakID = lcData.at[i, "peakID"]
                            for j in row["hits"]:
                                if j["well"] == peak_well and j["peakID"] == peakID:
                                    mass_conf_dict[i] = j["mass_conf"]

                        max_mass_conf = max(mass_conf_dict.values())

                        suitable_indexes = [i for i in well["orange"] if mass_conf_dict[i] == max_mass_conf]
                    else:
                        suitable_indexes = well["orange"]

                    sorted_list = sorted(well["orange"], key = lambda x: lcData.at[x, "area"], reverse = True)

                    peakAdded = False
                    for i in sorted_list:
                        if i in suitable_indexes and peakAdded == False:
                            final_result["green"].append(i)
                            row["comments"].append(f'Tentatively assigned peak with largest {mass_or_area} selected in '
                                                   'preference to other tentative assignments available for well '
                                                   f'{plate_tools.getUserReadableWell(lcData.at[i, "well"], self.plate_col_no)}.')
                            peakAdded = True
                        else:
                            final_result["discarded"].append(i)
                            row["comments"].append(f'<strong>The peak at {lcData.at[i, "time"]} in well '
                                                   f' {plate_tools.getUserReadableWell(lcData.at[i, "well"], self.plate_col_no)} '
                                                   f'was discarded as it had a smaller {mass_or_area} than an '
                                                   'otherwise equally likely alternative.</strong>')
                    
            row["final_result"] = final_result
            return row
        


        #Each set of functions is called in turn to validate the hits for each row
        #in the cpTable DataFrame. 
            
        self.cpTable = self.cpTable.apply(clusterRow, axis = 1)
        self.cpTable = self.cpTable.apply(selectCluster_ifrt, axis = 1)
        
        self.cpTable = self.cpTable.apply(refine_and_select, axis = 1)
        
        self.cpTable = self.cpTable.apply(indexClusterByWells, axis = 1)
        
        self.cpTable = self.cpTable.apply(finalResult, axis = 1)

        
    
    
    def removeDupAssigns(self, lcData, mass_or_area):
        """
        Checks each compound to ensure that no peak has been assigned 
        to two different compounds. If this has happened, the internalSTD
        (if present) takes first priority, limiting reactant is second priority,
        product third priority and finally a by-product is lowest priority.
        If a compound gives up a hit to a higher priority compound, a next-best option
        is sought from the list of hits that were discarded in the final result. 


        :param lcData: Pandas dataframe, containing columns called "type", "final_result" and "comments". 
        :param mass_or_area: A string, either "mass_conf" or "area". 

        :return: compoundDF as Pandas dataframe
        """
        
        def checkForOverlap(row):
            #Get a relevant slice of the cpTable to compare this row against
            indexes = []
            if row["type"] == "Byproduct":
                indexes = self.cpTable.loc[self.cpTable["type"].isin(["Product", "InternalSTD", "Reactant"])].index.values
            elif row["type"] == "Product":
                indexes = self.cpTable.loc[self.cpTable["type"].isin(["InternalSTD", "Reactant"])].index.values
            elif row["type"] == "Reactant":
                indexes = self.cpTable.loc[self.cpTable["type"].isin(["InternalSTD"])].index.values
            
            if len(indexes) > 0:
                new_slice = self.cpTable.loc[self.cpTable.index.isin(indexes)]
                
                to_keep = []
                for i, j in enumerate(row["final_result"]["green"]):
                    overlap = False
                    for assigned_peaks in list(new_slice["final_result"]):
                        if j in assigned_peaks["green"]:
                            overlap = True
                    
                    #If overlap was found, this hit is appended to this list of discarded hits         
                    if overlap:
                        row["comments"].append(f'<strong>The peak in well {plate_tools.getUserReadableWell(lcData.at[j, "well"], self.plate_col_no)} '
                                               f'at {lcData.at[j, "time"]} minutes was deselected '
                                               'because it was already assigned to the SM, internalSTD or by-product.</strong>')
                        possibilities = []
                        for replacement in row["final_result"]["discarded"]:
                            if lcData.at[j, "well"] == lcData.at[replacement, "well"]:
                                possibilities.append(replacement)
                                
                        if len(possibilities) > 0:
                            if mass_or_area == "mass_conf":
                                mass_conf_dict = {}
                                for lcData_index in possibilities:
                                    well = lcData.at[lcData_index, "well"]
                                    peakID = lcData.at[lcData_index, "peakID"]
                                    for hit in row["hits"]:
                                        if hit["well"] == well and hit["peakID"] == peakID:
                                            mass_conf_dict[i] = hit["mass_conf"]

                                #Get the best replacement option based on the mass confidence of that option
                                best_option = max(mass_conf_dict, key=mass_conf_dict.get)
                            else:
                                options = lcData.loc[lcData.index.isin(possibilities)]
                                options = options.loc[options["area"] == options["area"].max()]
                                best_option = list(options["area"])[0]
                            
                            to_keep.append(best_option)
                            row["comments"].append(f'<strong>The peak in well {plate_tools.getUserReadableWell(lcData.at[best_option, "well"], self.plate_col_no)} '
                                                   f'at {lcData.at[best_option, "time"]} minutes was promoted to replace the '
                                                   f'peak that was deselected, and was selected for promotion as it had the largest '
                                                   f'{mass_or_area}.</strong>')
                        
                        #finally, append the now unwanted peak into the discard pile
                        row["final_result"]["discarded"].append(j)   
                            
                    else:
                        to_keep.append(j)
                
                row["final_result"]["green"] = to_keep
            return row
        
        self.cpTable = self.cpTable.apply(checkForOverlap, axis = 1)
    
    def findPotentialConflicts(self):
        """
        Find any products which are in danger of overlapping with other compounds that
        are expected in the same well. 
        This function is useful as the program may not detect SM (etc) for that peak/well
        because the mass_conf is too low, but the user should still be warned that there 
        may be a problem. 

        :param compoundDF: Pandas datatable for compounds
        :return: A text output for that compound. 
        """
        
        def getConflicts(row):

            text = ""


            new_slice = self.cpTable.loc[((( self.cpTable["time"] - row["time"]).abs() < 0.02)
                                        & ( self.cpTable["name"] != row["name"])
                                        & ( self.cpTable["time"] != 0))]
            text = "No potential conflicts found."
            if len(new_slice.index) > 0:                            
                new_slice["intersection"] = new_slice["locations"].apply(lambda x: len(set(x).intersection(set(row["locations"]))))

                close_compounds = new_slice.loc[new_slice["intersection"] > 0].apply(lambda brow: brow["name"], axis = 1)
                                                

                if len(close_compounds) > 0:
                    text = "<strong>The compound was found to have a similar retention time as "
                    for i, name in enumerate(close_compounds):
                        if i != len(close_compounds)-1:
                            text = f'{text}{name}, '
                        elif len(close_compounds) > 1:
                            text = f'{text}and {name}. '
                        else:
                            text = f'{text}{name}. '
                    text = text + "The user should check this result manually.</strong>"

            row["conflicts"] = text
            return row
        if "time" in self.cpTable:
            self.cpTable = self.cpTable.apply(getConflicts, axis = 1)
        else:
            print("Necessary column, 'time', not present in dataframe. Run 'setBestWell' and 'setBestTime' first.")
    
    
    def setBestWell(self, outputTable, lcData, sampleTable, plot_type):

        def bestWellPerRow(row):


            if row["type"] == "Product":
                
                #Filter to the rows where the compound is found, then get a list of wells that have the 
                #maximum observed value for the given plot type (note: not guaranteed to return a single well!)
                new_slice = outputTable.loc[outputTable["well_no"].isin(row["locations"])]
                best_well_poss = list(new_slice.loc[(new_slice[plot_type] == new_slice[plot_type].max()) & 
                                    (new_slice[plot_type] > 0),  "well_no"])

                if len(best_well_poss) > 0:
                    best_well = best_well_poss[0]
                
                #If no "best well" is found (perhaps because the internal standard wasn't found in that well, 
                #set the best well to be the one where the highest Parea was observed. 
                else:
                    best_well_poss2 = list(new_slice.loc[(new_slice["Parea"] == new_slice["Parea"].max()) & 
                                                (new_slice["Parea"] > 0), "well_no"])
                    
                    if len(best_well_poss2) > 0:
                        best_well = best_well_poss2[0]
                    else:
                        best_well = -1

            else:
                new_slice = lcData.loc[lcData.index.isin(row["final_result"]["green"])]                 
                new_slice = new_slice.loc[new_slice["area"] == new_slice["area"].max(), "well"]
                    
                if len(list(new_slice)) > 0:
                    best_well = list(new_slice)[0]
                else:
                    best_well = -1
                
            row["best_wellno"] = best_well
            row["best_well"] = plate_tools.getUserReadableWell(best_well, self.plate_col_no) if best_well > -1 else "None found."

            best_sample_filename = list(sampleTable.loc[sampleTable["well"] == best_well, "filename"])
            if len(best_sample_filename) > 0:
                row["best_sample_filename"] = best_sample_filename[0]
            else:
                row["best_sample_filename"] = "None found."
            return row

        self.cpTable = self.cpTable.apply(bestWellPerRow, axis = 1)
    
    def setBestMS(self, lcData):
        
        def bestMSPerRow(row):
            if row["type"] != "Impurity":
                best_well = row["best_wellno"]
                peakID = -1
                for i in row["final_result"]["green"]:
                    if lcData.at[i, "well"] == best_well:
                        peakID = lcData.at[i, "peakID"]

                MS_plus = "-"
                MS_minus = "-"
                if peakID > -1:
                    for hit in row["hits"]:
                        if hit["well"] == best_well and hit["peakID"] == peakID:
                            MS_plus = hit["MS_plus"]
                            MS_minus = hit["MS_minus"]
                row["mass+"] = MS_plus
                row["mass-"] = MS_minus
            return row
        
        self.cpTable = self.cpTable.apply(bestMSPerRow, axis = 1)

    def setBestTime(self, lcData):

        def bestTimePerRow(row): 
            if len(row["final_result"]["green"]) >0:
                best_well = row["best_wellno"]
                best_time = lcData.loc[(lcData.index.isin(row["final_result"]["green"])) 
                                                    & (lcData["well"] == row["best_wellno"]), "time"].values[0] 
                
                row["time"] = best_time
            else:
                row["time"] = 0
            return row
        
        self.cpTable = self.cpTable.apply(bestTimePerRow, axis = 1)
    
    def setBestPurity(self, lcData):
        def bestPurityPerRow(row): 
            if len(row["final_result"]["green"]) >0:
                best_well = row["best_wellno"]
                best_purity = lcData.loc[(lcData.index.isin(row["final_result"]["green"])) 
                                        & (lcData["well"] == row["best_wellno"]), "area"].values[0] 
                row["best_purity"] = best_purity
            else:
                row["best_purity"] = 0
            
            return row
        
        self.cpTable = self.cpTable.apply(bestPurityPerRow, axis = 1)


    def findOverlap(self, lcData):
        """
        This function will determine whether the assigned peak in the "best well"
        happens to overlap with any other peaks in the lc for that well, based on how the pStart and pEnd
        times reported. 
        
        :param lcData: Pandas dataframe containing lcData, taken from a PyParse import (e.g. getWatersData.py -> rawWatersData.rawDADTable) 
        """
        def overlapPerRow(row):
            
            new_slice = lcData.loc[(lcData.index.isin(row["final_result"]["green"])) & (lcData["well"] == row["best_wellno"])]
            if len(new_slice.index) > 0:
                overlapping_peaks = lcData.loc[((lcData["pStart"] == new_slice.iloc[0]["pEnd"]) | (lcData["pEnd"] == new_slice.iloc[0]["pStart"])) 
                                                & (lcData["well"] == row["best_wellno"])]

                if len(overlapping_peaks.index) > 0:
                    row["overlaps"] = "<strong>Peak overlap detected!</strong>"
                else:
                    row["overlaps"] = "No peak overlap detected."
            else: 
                row["overlaps"] = "No peak overlap detected."
            return row

        if "best_wellno" in self.cpTable:
            self.cpTable = self.cpTable.apply(overlapPerRow, axis = 1)  
        else: 
            print("Necessary column, 'best_wellno', not present in dataframe. Run 'setBestWell' first.")

    def findImpurities(self, lcData, msData, time_abs_tol = 0.025, min_no_of_wells = 5, mass_abs_tol = 0.5):

        """
        Find unassigned peaks that have similar retention times and added them as compounds
        to the cpTable. 

        :param lcData: Pandas Dataframe
        :param msData: Pandas Dataframe
        :param time_abs_tol: float (optional)
        :param min_no_of_wells: integer (optional)
        :param mass_abs_tol: float (optional)

        :return: Updated cpTable with additional rows, one per common impurity
        """
            
        #Define function to get a list of wells from a cluster, 
        #without any duplicates
        def unique_wells(cluster):
            unique_wells = lcData.loc[lcData.index.isin(cluster), "well"].nunique()
            return unique_wells
        
        #Define sub-function to find common ions for a cluster of peaks
        def findCommonIons(valid_combos, MStype):
            
            #get all rows in msData that have the right combonations of both well and peakID
            result = (msData.loc[(msData["MStype"] == MStype) & 
                                (msData.set_index(["well", "peakID"]).index.isin(valid_combos))]).sort_values("time")


            #We need to find which m/z were commonly observed across this set of peaks. 
            #In this sorted data frame, find the difference between each successive row, 
            #mark any difference where the it is larger than mass_abs_tol as True,
            #then apply cumulative sum that counts the number of "Trues" seen so far
            #which is effectively a grouping index. 
            groupings = (result.groupby(result["MSvalue"].diff()
                                        .gt(mass_abs_tol)
                                        .cumsum())
                                        .apply(lambda row: list(row.index)))


            #Filter the ions down to those which appear in at least 80% of the peaks    
            filtered = [group for group in groupings if len(group) > len(cluster)*0.8]

            #Filter to just the top 5 most commonly observed peaks, to guarantee a maximum
            #number of 5 mass ions reported. 
            filtered = sorted(groupings, key = lambda x: len(x))[0:5]
            
            #get the mean value for each grouping, then return as a list
            return_ions = []
            for group in filtered:
                return_ions.append(msData.loc[msData.index.isin(group),"MSvalue"].mean())

            return return_ions
        
        #Get a list of all the peaks that have been assigned to a compound
        assigned_peaks = []
        self.cpTable.apply(lambda row: assigned_peaks.extend(row["final_result"]["green"]), axis = 1)
        
        #And therefore find a listof all the peaks that have not been assigned!
        unassigned_peaks = lcData.loc[lcData.index.isin(assigned_peaks) == False].index
        
        [clusters, clusterbands] = self.clusterHits(unassigned_peaks, lcData, time_abs_tol)
        
        #filter to only those clusters where there are sufficient number of unique wells
        clusters = [i for i in clusters if unique_wells(i) > min_no_of_wells]
        
        impurities = []
        
        for index, cluster in enumerate(clusters):
            
            new_slice = lcData.loc[lcData.index.isin(cluster)]
            new_slice = new_slice.loc[new_slice["area"] == new_slice["area"].max(), ["well", "peakID", "area"]]
            mean_rt = round(lcData.loc[lcData.index.isin(cluster), "time"].mean(), 4)
            
            valid_combos = [(lcData.at[peak, "well"], lcData.at[peak, "peakID"]) for peak in cluster]
            mass_plus = findCommonIons(valid_combos, "+")
            mass_minus = findCommonIons(valid_combos, "-")
            
            comments = []
        
            #Criteria to determine when a comment is added:
            #   -Cluster typically occurs in a particular column
            #   -Cluster typically occurs in a particular row
            #   -Cluster typically occurs for a particular compound in the platemap
            #   -Cluster is independant of position (i.e. whole plate) - this overrides all of the above
            #   -describe how many wells contained this impurity
            observed_unique_wells = unique_wells(cluster)
            comments.append(f'This impurity was observed in {observed_unique_wells} wells.')

            if observed_unique_wells == self.plate_row_no * self.plate_col_no:
                comments.append("This impurity was observed in every well of the plate.")
            elif observed_unique_wells > 0.5 * self.plate_row_no * self.plate_col_no:
                comments.append("This impurity was observed across the majority of the plate.")
            else:
                #Build a matrix of where the compound was observed, then iterate through each row/column in turn?
                columns = {}
                rows = {}
                for i in cluster:
                    wellno = lcData.at[i, "well"]
                    row = math.floor((wellno-1)/self.plate_col_no) + 1
                    column = ((wellno-1) % self.plate_col_no) + 1
                    if row in rows:
                        rows[row] = rows[row] + 1
                    else:
                        rows[row] = 1
                    if column in columns:
                        columns[column] = columns[column] + 1
                    else:
                        columns[column] = 1

                for cindex, column in columns.items():
                    if column > 0.8 * self.plate_row_no:
                        comments.append(f'Impurity is frequently observed in column {cindex}.')
                for rindex, row in rows.items():
                    if row > 0.8 * self.plate_col_no:
                        comments.append(f'Impurity is frequently observed in row {chr(ord("@")+(rindex)+1)}.')
                if len(columns.keys()) < 0.5 * self.plate_col_no:
                    readable_cols = [str(i) for i in sorted([int(j) for j in columns.keys()])]
                    comments.append(f'Impurity was only observed in columns {", ".join(readable_cols)}.')
                if len(rows.keys()) < 0.5 * self.plate_row_no:
                    readable_rows = sorted([chr(ord("@")+(x)) for x in rows.keys()])
                    comments.append(f'Impurity was only observed in rows {", ".join(readable_rows)}.')
                    
            #Append this data to a new dictionary, to ultimately build into a dataframe that will be concatonated 
            #with self.cpTable
            entry = {}
            #get the columns that should be filled 
            for column in self.cpTable.columns:
                entry[column] = ""
                
            entry["name"] = f'Impurity{index}'
            entry["locations"] = list(lcData.loc[lcData.index.isin(cluster), "well"].unique())
            entry["canonSMILES"] = "[*]"
            entry["stdname"] = f'Impurity{index}'
            entry["type"] = "Impurity"
            [entry["mass1"], entry["mass2"], entry["mass3"]] = [0,0,0]
            entry["time"] = mean_rt
            entry["mass-"] = mass_minus
            entry["mass+"] = mass_plus
            entry["final_result"] = {
                "green": cluster,
                "discarded": []
            }
            entry["cluster_bands"] = clusterbands
            entry["comments"] = comments
            impurities.append(entry)
        
        
        
        if len(clusters) > 0:
            #Turn the list of impurities into a dataframe
        
            df = pd.DataFrame(impurities)
            df.index = list(df["stdname"])
            df = df.drop("stdname", axis=1)
            self.cpTable = pd.concat([self.cpTable, df], axis=0)
                    
                
            
            
