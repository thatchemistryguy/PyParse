import math
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
import pandas as pd

import plate_tools

def plotChroma(row, cpTable, lcData, msData, traceData, save_dir, plate_col_no, filter_by_rt = []):
    
    #cpname, wellno, trace, pStart, pEnd, annotate_peaks, save_dir, 
    #           ms_plus, ms_minus, mass1):
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
    
    
    #Function to plot the mass spectrometric data
    def plotMS(axes, data, title):
        x = data[0]
        y = data[1]
        
        axes.bar(x, y, width = 2)
        axes.set_title(title)
        axes.set_ylim(0, 200)
        
        last_annotation = 0
        highest_mass = 0
        for i, j in enumerate(x):
            if y[i] > 20:
                if math.isclose(j, last_annotation, abs_tol = 3):
                    axes.annotate(j, [j+15, y[i]+20], ha="center", 
                                  va="bottom", rotation=90, size = 8)
                    axes.arrow(j+10, y[i]+18, -10, -8)
                else:
                    axes.annotate(j, [j+1, y[i]], ha="center", 
                                  va="bottom", rotation=90, size = 8)
                last_annotation = j
                if j > highest_mass:
                    highest_mass = j
        
        #Configure mass spec axis to reach at least 1000, but higher if necessary
        #With some additional headroom so that peaks are not obscured by the edge
        #of the plot
        if highest_mass > 950:
            axes.set_xlim(0, math.ceil((highest_mass + 100)/100)*100) 
        else:
            axes.set_xlim(0, 1000)
            
    def getAnnotations(annot_row, well, cpName):
        if (annot_row["type"] != "Impurity" or annot_row["name"] == cpName):

            data = [i for i in annot_row["final_result"]["green"] if lcData.at[i, "well"] == well]
            if len(data) > 0:
                annotations.append({
                    "cpname": annot_row["name"], 
                    "time": lcData.at[data[0], "time"]
                    })
        
    #Get plotting space and edit the default plot dimensions
    plt.rcParams["figure.figsize"] = (12, 6)
    fig, (a0, a1, a2) = plt.subplots(3, 1, gridspec_kw={'height_ratios': [2, 2, 6]})

    #First, determine for the compound given which well it is that should be plotted
    #and which peak should be highlighted with MS data plotted. 
    well_to_use = -1
    if len(row["final_result"]["green"]) > 0:
        well_to_use = row["best_wellno"]
        peakID = -1
        for i in row["final_result"]["green"]:
            if lcData.at[i, "well"] == well_to_use:
                peakID = lcData.at[i, "peakID"]
                
    elif len(row["final_result"]["discarded"]) == 1:
        well_to_use = lcData.at[row["final_result"]["discarded"][0], "well"]
        peakID = -1
        for i in row["final_result"]["green"]:
            if lcData.at[i, "well"] == well_to_use:
                peakID = lcData.at[i, "peakID"]
    elif len(row["locations"]) == 1:
        well_to_use = row["locations"][0]
        peakID = -1
        
    #As long as a well was found, proceed to plot an annotated chromatogram for this compound.     
    if well_to_use > -1:
        
        #Get annotations based on any peaks in this well that were assigned to a compound. 
        annotations = []
        cpTable.apply(getAnnotations, args=(well_to_use, row["name"]), axis = 1)
        
        #If a specific peak was selected, get the relevant MS_plus and MS_minus rows for plotting. 
        if peakID > -1:
            plus_slice = msData.loc[(msData["well"] == well_to_use) & (msData["peakID"] == peakID) & (msData["MStype"] == "+")]
            ms_plus = [list(plus_slice["MSvalue"]), list(plus_slice["MSintensity"])]

            minus_slice = msData.loc[(msData["well"] == well_to_use) & (msData["peakID"] == peakID) & (msData["MStype"] == "-")]
            ms_minus = [list(minus_slice["MSvalue"]), list(minus_slice["MSintensity"])]

            lc_slice = lcData.loc[(lcData["well"] == well_to_use) & (lcData["peakID"] == peakID)]
            
            #Get the time where the peak starts and when it ends, so that we can use matplotlib 
            #to shade the peak in red. 
            pStart = float(list(lc_slice["pStart"])[0])
            pEnd = float(list(lc_slice["pEnd"])[0])
        
        #If a specific peak wasn't selected, provide the necessary variable structure. 
        else:
            ms_plus = []
            ms_minus = []
            pStart = 0
            pEnd = 0

        #Plot both MS- and MS+    
        plotMS(a0, ms_minus, "MS-")
        plotMS(a1, ms_plus, "MS+")

        #Plot the UV chromatogram
        trace = ([list(traceData.loc[traceData["well"] == well_to_use, "time"]), 
                 list(traceData.loc[traceData["well"] == well_to_use, "height"])])
        
        a2.plot(trace[0], trace[1])

        #label the graph and axes
        label = plate_tools.getUserReadableWell(well_to_use, plate_col_no)
        fig.suptitle(f'{row["name"]} ({row["mass1"]}): Well {label}')
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
        annotations.sort(key = lambda x: x["time"])
        last_x_position = 0
        max_x = max(trace[0])
        for index, i in enumerate(annotations):

            x_index = min(enumerate(trace[0]), 
                          key = lambda x: abs(x[1]-i["time"]))[0]
            if index < len(annotations) / 2:
                h_align = "right"
            else:
                h_align = "left"
            if index == 0:
                offset = -(max_x / 40)
            elif index == len(annotations)-1:
                offset = (max_x / 40)
            else:
                offset = 0
            arrow_props = dict(arrowstyle="->", connectionstyle="arc3")
            if math.isclose(last_x_position, i["time"]+offset, abs_tol = (max_x / 20)) and index != len(annotations)-1:
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

        #plot a hatched region if the filter_by_rt option 
        #has been used
        if len(filter_by_rt) > 0:

            #Set default yvalues to match that of chromatogram trace
            hatched_ymin, hatched_ymax = [],[]
            for i in trace[1]:
                hatched_ymin.append(i)
                hatched_ymax.append(i)

            for new_range in filter_by_rt:
                minx = float(new_range.split("-")[0].strip())
                maxx = float(new_range.split("-")[1].strip())
                for i, j in enumerate(trace[0]):
                    #Amend any values such between specified ranges
                    if j > minx and j < maxx:
                        hatched_ymax[i] = 130
                        hatched_ymin[i] = min_y - 10
            #fill between on graph, using transparency = 0.3
            a2.fill_between(trace[0], hatched_ymin, hatched_ymax, hatch = "/", alpha = 0.3)

        #Fill to highlight the peak of interest 
        if pStart != 0 and pEnd != 0:
            a2.fill_between(trace[0], second_curve, trace[1], color="red")

        #Set the number of tick points for each axis. 
        a2.locator_params(axis='y', nbins=6)
        a2.locator_params(axis='x', nbins=20)

        #Save graph to output directory. 
        plt.savefig(f'{save_dir}graphs/chroma-{row["name"]}-best.jpg', format="jpg")
        plt.close()

def plotHitValidationGraph(row, lcData, save_dir, plate_col_no, plate_row_no):
    """
    Plots all the hit peaks in a scatter graph of peaktime vs
    well, colour coded by whether the hit was included or discarded
    from the final output. 
    
    :param row: pandas dataframe row from a PyParse compound table
    :param validatedHits: pandas dataframe containing the liquid chromatography 
    :param save_dir: File directory for where to save the matplotlib figure
    :param cluster_bands: A list of the average retention times for each cluster.
    
    :return: jpg of the hit validation graph saved to output directory 
    """  
    
    #Find the largest peak area, and the minimum/maximum retention times
    #of the hits to be plotted
    unused_hits = row["final_result"]["discarded"]
    green_hits = row["final_result"]["green"]
    for cluster in row["discarded_clusters"]:
        for key, value in cluster.items():
                unused_hits = unused_hits + value
    
    all_hits = green_hits + unused_hits

    full_slice = lcData.loc[lcData.index.isin(all_hits)]
    max_size = full_slice["area"].max() if len(all_hits) > 0 else 0
    max_time = full_slice["time"].max() if len(all_hits) > 0 else 0
    min_time = full_slice["time"].min() if len(all_hits) > 0 else 0

    
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
        total_wells = plate_col_no * plate_row_no
        ax.set_xlim(-total_wells*0.2, total_wells*1.1)
        ax.set_ylim(y_lims[0] - (y_lims[1]-y_lims[0])/10, y_lims[1] + (y_lims[1]-y_lims[0])/10)
       
        #Plot the green points
        green_slice = lcData.loc[lcData.index.isin(green_hits)]
        ax.scatter(green_slice["well"], green_slice["time"], s = green_slice["area"].apply(lambda x: math.ceil(x * size_ratio)**2),
                   marker = "o", color = "black")
        
        #Plot the unused points
        unused_slice = lcData.loc[lcData.index.isin(unused_hits)]
        ax.scatter(unused_slice["well"], unused_slice["time"], s = unused_slice["area"].apply(lambda x: math.ceil(x * size_ratio)**2),
                   marker = "v", color = "red")

        """
        #Plot all hits that were used in the analysis
        for i in row["final_result"]["green"]:
            #find which cluster this hit came from
            for index, band in enumerate(row["cluster_bands"]):
                if i in band[1]:
                    cluster_id = index
            if cluster_id > len(t_markers)-1:
                thismarker = "o"
            else:
                thismarker = t_markers[cluster_id]
                
            ax.scatter(lcData.at[i, "well"], lcData.at[i, "time"], s = math.ceil(lcData.at[i, "area"] * size_ratio)**2, 
                        marker = thismarker, color="black")
        
        #Plot all hits that were discarded in favour of a better hit
        for i in row["final_result"]["discarded"]:
            #find which cluster this hit came from
            for index, band in enumerate(row["cluster_bands"]):
                if i in band[1]:
                    cluster_id = index
            if cluster_id > len(t_markers)-1:
                thismarker = "o"
            else:
                thismarker = t_markers[cluster_id]
            ax.scatter(lcData.at[i, "well"], lcData.at[i, "time"], s = math.ceil(lcData.at[i, "area"] * size_ratio)**2, 
                        marker = thismarker, color="red")
        
        #Plot all hits that were discarded by cluster
        for cluster in row["discarded_clusters"]:   
            for key, value in cluster.items():
                for i in value:
                    #find which cluster this hit came from
                    for index, band in enumerate(row["cluster_bands"]):
                        if i in band[1]:
                            cluster_id = index
                    if cluster_id > len(t_markers)-1:
                        thismarker = "o"
                    else:
                        thismarker = t_markers[cluster_id]
                    ax.scatter(lcData.at[i, "well"], lcData.at[i, "time"], s = math.ceil(lcData.at[i, "area"] * size_ratio)**2, 
                                marker = thismarker, color="red")
        """
        #Annotate the graph with the average retention time of each cluster
        for i in range(len(row["cluster_bands"])):
            x = -total_wells*0.17
            y = row["cluster_bands"][i][0]
            ax.annotate(f'C{i}: {row["cluster_bands"][i][0]} min.', [x, y])

        #Label the graph and axes, and save to the output directory.   
        plt.title(row["name"])
        plt.xlabel("Well")
        plt.ylabel("Retention Time /min")
        
        plt.savefig(f'{save_dir}graphs/hits-{row["name"]}.jpg', format="jpg")
        plt.close()

def plotHeatmaps(outputTable, save_dir, plate_row_no, plate_col_no):
    """
    Plots and saves heatmaps for the full dataset
    
    :param outputTable: a pandas datatable
    :param save_dir: a string for the output directory
    
    :return: jpg of the heatmap saved to output directory
    """
    #Set standard matplotlib graph size
    plt.rcParams["figure.figsize"] = (12, 6)
    
    zvalues = {
        "SMarea": "SMarea",
        "Parea": "Parea", 
        "conversion": "P/SM+P", 
        "ratio_to_IS": "P/STD",
        "corrSMarea": "corrSMarea",
        "corrParea": "corrParea", 
        "corrected_conversion": "corrP/SM+P",
        "corrected_ratio_to_IS": "corrP/STD"
    }
    for key, zvalue in zvalues.items():
        
        returnVal = []
        labels = []
        #convert output table to 2-D list that can be used to plot heatmap. 
        for index, row in outputTable.iterrows():
            rowVal = int(math.floor((index-1) / plate_col_no))
            colVal = int((index-1) % plate_col_no)
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
        xLabels = [i for i in range(1, plate_col_no+1)]
        yLabels = [chr(ord('@')+i) for i in range(1, plate_row_no+1)]
        
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

def genLocationHeatmaps(row, save_dir, plate_col_no, plate_row_no):
    """
    Generate a heatmap type graph to visualise the 
    expected locations of each compound, based on the platemap that was provided. 

    :param row: A row from Pandas dataframe containing describing a compound
    :param save_dir: string describing the location of the output folder
    :output: Heatmaps saved to the output folder. 
    """
    
    if row["type"] != "Impurity":
        location_data = []
        labels = []
        for i in range(plate_col_no * plate_row_no):
            rowVal = int(math.floor(i / plate_col_no))
            colVal = int(i % plate_col_no)
            
            if colVal == 0:
                location_data.append([])
                labels.append([])
            if i+1 in row["locations"]:
                location_data[rowVal].append(1)
                labels[rowVal].append(u'\u2713')
            else:
                location_data[rowVal].append(0)
                labels[rowVal].append("")
                

        pdTable = pd.DataFrame(location_data)
        xLabels = [i for i in range(1, plate_col_no+1)]
        yLabels = [chr(ord('@')+i) for i in range(1, plate_row_no+1)]
        
        ax = sns.heatmap(pdTable, xticklabels=xLabels, yticklabels=yLabels, cmap = "viridis", annot = labels, fmt="")
        ax.xaxis.set_ticks_position("top")
        
        #Save heatmap to output directory
        plt.savefig(f'{save_dir}graphs/loc_heatmap_{row["name"]}.jpg', format="jpg")
        plt.close()   

def plotPieCharts(zvalue, outputTable, save_dir, by_products, plate_row_no, plate_col_no):
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
        #Set standard matplotlib graph size
        plt.rcParams["figure.figsize"] = (12, 6)
        
        #declare new subplots
        fig, axs = plt.subplots(plate_row_no, plate_col_no)
        
        #get maximum value, so diameter of all pies can be calculated in relation
        #to the best performing well. 
        max_val = outputTable[zvalue].max()
        
        #format data and generate pie charts
        pies_baked = []
        if max_val != 0:
            
            for index, row in outputTable.iterrows():
                pies_baked.append(index)
                rowVal = math.floor((index-1) / plate_col_no)
                colVal = int((index-1) % plate_col_no)

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
                    
                
                if plate_row_no == 1:
                    axs[colVal].pie(data, 
                            textprops={"size": "smaller"}, 
                            colors = colors,  
                            radius=chart_size,
                            normalize = True)
                elif plate_col_no == 1:
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
        for i in range(1, plate_col_no * (plate_row_no) + 1):
            rowVal = math.floor((i-1) / plate_col_no)
            colVal = int((i-1) % plate_col_no)
            if i not in pies_baked:
                
                #matplotlib axes drop a dimension 
                #when that dimension is length = 1. 
                #Adjust code to match
                if plate_row_no == 1:
                    plt.delaxes(axs[colVal])
                elif plate_col_no == 1:
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
        if plate_row_no == 1:
            axs[plate_col_no-1].legend(lines,
                    labels, loc="lower center",
                    bbox_to_anchor=(1,1), ncol=math.ceil(len(labels)/3))
        elif plate_col_no == 1:
            axs[plate_row_no-1].legend(lines,
                    labels, loc="upper left",
                    bbox_to_anchor=(1,1), ncol=math.ceil(len(labels)/3))
        else:
            axs[plate_row_no-1, plate_col_no-1].legend(lines,
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