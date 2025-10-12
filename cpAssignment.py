import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
import logging
from statistics import mean
import math
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.cm as cm
import numpy as np

import plate_tools

class Assignment:
    def __init__(self, filename, plate_col_no, plate_row_no):   
        #read csv file into dataframe
        #replace empty cells with an empty string
        #convert all column names to lower case and remove whitespace
        self.inputCSV = pd.read_csv(filename)
        self.inputCSV.fillna("", inplace=True)
        self.inputCSV.columns = self.inputCSV.columns.str.strip().str.lower()
        
        self.plate_col_no = plate_col_no
        self.plate_row_no = plate_row_no
        
    #Fn to convert a well name like B5 to machine format (11)
    def convertWellToNum(self, wells):
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
    
    #Fn to convert smiles into canonicalised smiles
    def getCanonSmiles(self, smiles):
        mol = Chem.MolFromSmiles(smiles.strip())
        return Chem.MolToSmiles(mol)
            
    def generateCPTable(self):
        compound_list = []
        type_dic = {
            "desired product smiles": "Product",
            "limiting reactant smiles":"Reactant",
            "internalstd smiles": "InternalSTD"
        }
        counter = {
            "desired product smiles": 1,
            "limiting reactant smiles": 1,
            "byproduct": 1
        }
        
        #Canonicalise all the incoming smiles in key columns in case the user hadn't done so already
        for col in self.inputCSV:
            if col in type_dic or ("byproduct" in col and "smiles" in col):
                self.inputCSV[col] = self.inputCSV[col].apply(self.getCanonSmiles)
        
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
                            
                    new_entry = {
                        "canonSMILES": cpindex,
                        "type": cptype,
                        "locations": row["well"],
                        "name": name,
                        "rt": rt,
                        "comments": []
                    }
                    compound_list.append(new_entry)
        
        self.cpTable = pd.DataFrame(compound_list)
        
        #set the index to be the canonicalised smiles
        self.cpTable.index = list(self.cpTable["canonSMILES"])
        
        #convert the "A1" style well IDs into a integer, to allow matching to a well in the rawData
        self.cpTable["locations"] = self.cpTable["locations"].apply(self.convertWellToNum)
    

    
    def generateEMs(self, calc_boc):
        
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
                    return round(Descriptors.ExactMolWt(mol), 2) + 2
                else:
                    return 0
            
        self.cpTable["mass1"] = self.cpTable["canonSMILES"].apply(lambda smiles: getMW(smiles))
        
        if calc_boc == "True":
            
            smirks1 = "[NX3,n:1][C:2](=[O:3])[O:4][C]([CH3])([CH3])[CH3]>>[*:1][C:2](=[O:3])[O:4]"
            self.cpTable["mass2"] = self.cpTable["canonSMILES"].apply(lambda smiles: transform_and_getMW(smiles, smirks1, "mass2"))
            
            smirks2 = "[NX3,n:1][C](=[O])[O][C]([CH3])([CH3])[CH3]>>[*:1][H]"
            self.cpTable["mass3"] = self.cpTable["canonSMILES"].apply(lambda smiles: transform_and_getMW(smiles, smirks2, "mass3"))
            

    def findHits(self, msData, lcData, mass_abs_tol = 0.5, min_massconf_threshold = 10, 
                                      time_abs_tol = 0.025, calc_higherions = "True"):
        
        def getMatches(compound):
            
            def getMassWithHighestConf(row):

                MS_plus = list(df.loc[(df.index.isin(MS_hits)) & (df["MStype"] == "+") &
                                       (df["MSintensity"] == row["max_intensity"]), "MSvalue"])
                MS_minus = list(df.loc[(df.index.isin(MS_hits)) & (df["MStype"] == "-") &
                                       (df["MSintensity"] == row["max_intensity"]), "MSvalue"])
                row["MS_plus"] = MS_plus[0] if len(MS_plus) > 0 else "-"
                row["MS_minus"] = MS_minus[0] if len(MS_minus) > 0 else "-"
                
                return row
            
            df = msData.loc[msData["well"].isin(compound["locations"])]
            
            MS_hits = []
            
            for mass in [compound["mass1"], compound["mass2"], compound["mass3"]]:
                hits = df.loc[(df["MStype"] == "+") & 
                                   ((abs(df["MSvalue"] - (mass + 1.01)) <= mass_abs_tol) |
                                    (abs(df["MSvalue"] - (mass + 2.02)/2) <= mass_abs_tol) |
                                    (abs(df["MSvalue"] - (mass + 3.03)/3) <= mass_abs_tol))]
                MS_hits = MS_hits + list(hits.index.values)
                
                hits = df.loc[(df["MStype"] == "-") & 
                                   (abs(df["MSvalue"] - (mass - 1.01)) <= mass_abs_tol)]
                MS_hits = MS_hits + list(hits.index.values)
            
            grouped = (df[df.index.isin(MS_hits)].groupby(["well", "peakID"], as_index = False)
                        .agg(mass_conf = ("perc_intensity", "sum"), max_intensity = ("MSintensity", "max")))
            grouped = grouped.apply(getMassWithHighestConf, axis = 1)
           
            grouped = grouped.loc[grouped["mass_conf"] >= min_massconf_threshold]
            grouped = grouped.drop("max_intensity", axis=1)
            
            #fetch any peaks where the time matches the retention time provided
            
            new_slice = lcData.loc[lcData["well"].isin(compound["locations"]) & 
                                  (abs(lcData["time"] - compound["rt"]) < time_abs_tol), ["well", "peakID"]]
            #For these additional hits, set the requisite columns so that the two sets of hits can be merged into a 
            #single dataframe. 
            new_slice["mass_conf"] = 0
            new_slice["MS_plus"] = "-"
            new_slice["MS_minus"] = "-"

            return pd.concat([grouped, new_slice], axis=0).to_dict("records")                      

        self.cpTable["hits"] = self.cpTable.apply(getMatches, axis = 1)
    
    def clusterHits(self, hits, lcData, time_abs_tol = 0.025):
        
        df = lcData.loc[lcData.index.isin(hits)]
        df.sort_values("time", inplace = True)

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
        
    def validateHits(self, lcData, msData, uvData, time_abs_tol = 0.025, massconf_threshold = 0.5, uv_abs_tol = 10,
                    uv_cluster_threshold = 0.5, uv_match_threshold = 0.5, cluster_size_threshold = 0.8, min_no_of_wells = 5,
                    validate = "True", mass_or_area = "mass_conf"):
              
        def clusterRow(row):
                
            relevant_indexes = []
            for hit in row["hits"]:
                relevant_indexes = (relevant_indexes + list(lcData[(lcData["well"] == hit["well"]) & 
                                                                (lcData["peakID"] == hit["peakID"])].index.values)) 
            
            outcome = self.clusterHits(relevant_indexes, lcData, time_abs_tol)
            row["clusters"] = outcome[0]
            row["cluster_bands"] = outcome[1]

            return row
                       
        def selectCluster_ifrt(row):
            comments = row["comments"]
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

            :param cluster: list of dictionaries, where each dictionary is a hit
            :param comments: A list of comments for that structure so far.

            :return: List comprising [a dictionary for the refined cluster, list of comments]
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

            :param cluster: a dict, with list of dicts for each header
            :param comments: A list of comments for the compound so far

            :return: List comprising [a dictionary for the refined cluster, list of comments]
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

            :param cluster: a dict, with list of dicts for each header
            :param UVdatafound: boolean for whether the rpt data contains UV data
            :param comments: A list of comments for that structure so far

            :return: List comprising [a dictionary for the refined cluster, list of comments]
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

            :param clusters: a list of dictionaries, with list of dictionaries for each header
            :return refined_clusters: a list of dictionaries, with list of dictionaries for each header

            :return discarded_clusters: a list of dictionaries, with list of dictionaries for each header
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
            return row

        def selectClusterBySize(row):
            """
            If more than one cluster was found for the compound, 
            this function is called to try to select a single cluster based on 
            which cluster is the largest. If more than one cluster
            has a close-to-largest size, take them all. 

            :param clusters: a list of dictionaries, with list of dictionaries for each header
            :return refined_clusters: a list of dictionaries, with list of dictionaries for each header

            :return discarded_clusters: a list of dictionaries, with list of dictionaries for each header
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
            row["refined_clusters"] = refined_clusters2
            row["discarded_clusters"] = discarded_clusters
            return row
        
        def indexClusterByWells(row):
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
            
        def refine_and_select(row): 
            #If there are sufficient wells to perform refine and select a cluster, do so. 
            #Otherwise, simply mark all hits as "green" fill in the necessary table structure. 
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
        
        def finalResult(row):
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

        :param compoundDF: Pandas dataframe
        :param internalSTD: A string representing the name of the internal standard
        :param SMs: A list of starting material names
        :param products: A list of product names
        :param by_products: A list of by_product names

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
            
            close_compounds = new_slice.apply(lambda brow: brow["name"], axis = 1)
                                              

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
            else:
                text = "No potential conflicts found."

            row["conflicts"] = text
            return row
        if "time" in self.cpTable:
            self.cpTable = self.cpTable.apply(getConflicts, axis = 1)
        else:
            print("Necessary column, 'time', not present in dataframe. Run 'setBestWell' and 'setBestTime' first.")
    
    
    def setBestWell(self, outputTable, lcData, plot_type):

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
            best_well = row["best_wellno"]
            best_time = best_purity = lcData.loc[(lcData.index.isin(row["final_result"]["green"])) 
                                                 & (lcData["well"] == row["best_wellno"]), "time"].values[0] 
            
            row["time"] = best_time
            return row
        
        self.cpTable = self.cpTable.apply(bestTimePerRow, axis = 1)
    
    def setBestPurity(self, lcData):
        def bestPurityPerRow(row): 
            best_well = row["best_wellno"]
            best_purity = lcData.loc[(lcData.index.isin(row["final_result"]["green"])) 
                                     & (lcData["well"] == row["best_wellno"]), "area"].values[0] 
            row["best_purity"] = best_purity
            
            return row
        
        self.cpTable = self.cpTable.apply(bestPurityPerRow, axis = 1)

    def findOverlap(self, lcData):
        """
        This function will determine whether the assigned peak in the "best well"
        happens to overlap with any other peaks in the lc for that well, based on how the pStart and pEnd
        times reported. 
        
        :param lcData: Pandas daaframe containing lcData, taken from a PyParse import (e.g. getWatersData.py -> rawWatersData.rawDADTable) 
        """
        def overlapPerRow(row):
            
            new_slice = lcData.loc[(lcData.index.isin(row["final_result"]["green"])) & (lcData["well"] == row["best_wellno"])]
            overlapping_peaks = lcData.loc[((lcData["pStart"] == new_slice.iloc[0]["pEnd"]) | (lcData["pEnd"] == new_slice.iloc[0]["pStart"])) 
                                            & (lcData["well"] == row["best_wellno"])]

            if len(overlapping_peaks.index) > 0:
                row["overlaps"] = "<strong>Peak overlap detected!</strong>"
            else:
                row["overlaps"] = "No peak overlap detected."
            return row

        if "best_wellno" in self.cpTable:
            self.cpTable = self.cpTable.apply(overlapPerRow, axis = 1)  
        else: 
            print("Necessary column, 'best_wellno', not present in dataframe. Run 'setBestWell' first.")

    def findImpurities(self, lcData, msData, save_dir, time_abs_tol = 0.025, min_no_of_wells = 5, mass_abs_tol = 0.5):
        
        #Define function to get a list of wells from a cluster, 
        #without any duplicates
        def unique_wells(cluster):
            unique_wells = lcData.loc[lcData.index.isin(cluster), "well"].nunique()
            return unique_wells
        
        #Define sub-function to find common ions for a cluster of peaks
        def findCommonIons(valid_combos, MStype):
            
            #get all rows in msData that have the right combonations of both well and peakID
            result = (msData.loc[(msData["MStype"] == MStype) & 
                                 (msData.set_index(["well", "peakID"]).index.isin(valid_combos))])

            result.sort_values("time", inplace = True)

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
            entry["type"] = "Impurity"
            [entry["mass1"], entry["mass2"], entry["mass3"]] = [0,0,0]
            entry["mass-"] = mass_minus
            entry["mass+"] = mass_plus
            entry["final_result"] = {
                "green": cluster,
                "discarded": []
            }
            entry["cluster_bands"] = clusterbands
            entry["comments"] = comments
            
            impurities.append(entry)
        
        #Turn the list of impurities into a datafram
        
        df = pd.DataFrame(impurities)
        df.index = list(df["name"])
        
        self.cpTable = pd.concat([self.cpTable, df], axis=0)
            
        #If 1 or more impurities were found, plot a hit validation graph to display these. 
        #Set standard matplotlib graph size
        plt.rcParams["figure.figsize"] = (12, 6)
        
        if len(clusters) > 0:
            fig, ax = plt.subplots()
            total_wells = self.plate_col_no * self.plate_row_no
            ax.set_xlim(-total_wells*0.25, total_wells*1.1)
            #Plot a hit validation graph for all clusters to more easily display the results.
            #Save this to the graphs folder. 
            palette = cm.rainbow(np.linspace(0, 1, len(clusters)))

            for index, cluster in enumerate(clusters):
                x = []
                y = []
                for i in cluster:
                    x.append(lcData.at[i, "well"])
                    y.append(lcData.at[i, "time"])
                plt.scatter(x, y, color = palette[index])

                #Annotate the graph with the average retention time of each cluster

                mean_rt = round(lcData.loc[lcData.index.isin(cluster), "time"].mean(), 2)

                ax.annotate(f'Imp{index}: {mean_rt} min.', [-total_wells*0.22, mean_rt], 
                            verticalalignment='center')     

            #Label the graph and axes, and save to the output directory.   
            plt.title("All Frequent Impurities Found")
            plt.xlabel("Well")
            plt.ylabel("Retention Time /min")
            plt.savefig(f'{save_dir}/graphs/impuritychart.jpg', format = "jpg")    
            plt.close()
            
        
        
