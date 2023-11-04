"""
The .daml files provided by Shimadzu come in an xml format
where key data pieces (chromatogram, observed m/z ratios 
with their corresponding intensity and UV traces) are 
encoded in base64 format. 

The following script extracts the data relevant to PyParse, and reformats
it into a common structure. 

"""

import xml.etree.ElementTree as ET
import base64
import glob
import math


def init(args):
    global options
    options = args


def getData(input_dir):
    wellData = []
    masterTable = {}
    chroma = {}
    sample_IDs = {}

    input_files = glob.glob(input_dir + '*.daml')

    for filename in input_files:
        #Open the .daml file, which is structured as an xml
        tree = ET.parse(filename)
        root = tree.getroot()

        #If the number of rows or number of columns has not yet been 
        #specified, determine these values from the .daml file
        if options.plate_col_no == 0 or options.plate_row_no == 0:
            vessels = root.find("datafile_attribs").find("Plate").find("Vessels")
            #determine the number of columns/rows by how many times 
            #a well has the y-value = 0 or x-val = 0 respectively

            for vessel in vessels.findall("Vessel"):
                if int(vessel.get("x")) == 0:
                    options.plate_row_no += 1
                if int(vessel.get("y")) == 0:
                    options.plate_col_no += 1

        #get the well number
        wellno = int(root.find("datafile_attribs").find("vial").text)

        #get the sampleID
        sample_IDs[wellno] = str(root.find("datafile_attribs").find("desc").text)

        #Get to the results secton of the .daml file
        results = root.find("results")

        peaks = {}
        #start by getting the peaks observed in the UV spectrum
        raw_UV = results.find("lc_chros").find("chro")

        for peak in raw_UV.findall("peak"):
            #Fetch all data about a peak from the details section
            new_entry = {}
            for detail in peak.get("details").split(","):
                name = detail.strip().split("=")[0]
                value = detail.strip().split("=")[1]
                if name == "Peak#":
                    new_entry["peakID"] = int(value)
                elif name == "RT":
                    new_entry["time"] = float(value.split("mins")[0].strip())
                    new_entry["pStart"] = float(value.split("(")[1].split("-")[0])
                    new_entry["pEnd"] = float(value.split("-")[1].split("mins")[0])
                elif name == "Area":
                    new_entry["areaAbs"] = float(value)
                elif name == "Area%":
                    new_entry["area"] = float(value)
            #Add placeholders for the correct data structure
            new_entry["MS+"] = [[]]
            new_entry["MS-"] = [[]]

            #Add the new entry to the peaks dictionary
            peaks[new_entry["peakID"]] = new_entry

        #get the chromatographic trace for this well

        coded_trace = raw_UV.find("data").text
        decoded_trace = base64.b64decode(coded_trace)
        trace_list = str(decoded_trace).split(";;;")[1].split("'")[0].split(";")

        xval = []
        yval = []
        for i in range(0, len(trace_list), 2):
            xval.append(float(trace_list[i])/60000)
            yval.append(float(trace_list[i+1]))

        #Normalise the heights of the peaks so that they plot correctly

        max_y_val = max(yval)
        normalised_yval = [100*i/max_y_val for i in yval]
        chroma[wellno] = [xval, normalised_yval]

        #next, get the mass spec data for each peak, and match it to the UV peaks

        raw_ms = results.find("spectra")
        for spec in raw_ms.findall("spec"):
            peakID = int(spec.find("peak_num").text)
            if peakID in peaks.keys():
                msdata = spec.find("data").text
                decoded = base64.b64decode(msdata)

                #determine whether this data is for MS+ or MS-
                ms_type = "MS+" if str(decoded).find("MS(+)") > -1 else "MS-"

                #Get the list of masses with their intensities, and group them
                #pairwise ([massion, intensity])
                decoded_list = str(decoded).split(";;;")[1].split("'")[0].split(";")

                masses = []
                total = 0
                for i in range(0, len(decoded_list), 2):
                    datapair = [float(decoded_list[i])/10000, float(decoded_list[i+1])]
                    masses.append(datapair)
                
                #normalise the intensities 
                max_intensity = max(masses, key = lambda x: x[1])[1]
                normalised_masses = [[i[0], round(i[1]*100/max_intensity, 3)] for i in masses]
                
                #Remove any masses which are so small as to be negligable
                total = sum([i[1] for i in normalised_masses])
                refined_masses = [i for i in normalised_masses if math.floor(i[1]*100/total) > 0]

                if ms_type == "MS+":
                    peaks[peakID]["MS+"] = refined_masses
                else:
                    peaks[peakID]["MS-"] = refined_masses

        #finally, fetch the UV absorption maxima values

        raw_UV_ads = results.find("pda_spectra")
        for pda_spec in raw_UV_ads.findall("pda_spec"):
            #Get the peakID so the pda data is matched to the correct peaks
            peakID = int(pda_spec.find("peak_num").text)
            if peakID in peaks.keys():
                
                UV_maxima = []

                #Decode the base64 data, and group the data into pairs of x/y values
                UV_pairs = []
                pda_data = base64.b64decode(pda_spec.find("data").text)
                decoded_list = str(pda_data).split(";;;")[1].split("'")[0].split(";")
                for i in range(0, len(decoded_list), 2):
                        datapair = [float(decoded_list[i])/100, float(decoded_list[i+1])]
                        UV_pairs.append(datapair)

                #Tag the very first value as a maxima if required
                if UV_pairs[0][1] > UV_pairs[1][1] and UV_pairs[0][1] > options.min_uv_threshold:
                    UV_maxima.append(UV_pairs[0][0])
                
                #Iterate through all other n-1/n/n+1 comparisons to find all other maxima
                for i in range(1, len(UV_pairs)-1):
                    if (UV_pairs[i][1] > UV_pairs[i-1][1] and UV_pairs[i][1] > UV_pairs[i+1][1] 
                        and UV_pairs[i][1] > options.min_uv_threshold):
                        UV_maxima.append(UV_pairs[i][0])
                
                if UV_pairs[-1][1] > UV_pairs[-2][1] and UV_pairs[-1][1] > options.min_uv_threshold:
                    UV_maxima.append(UV_pairs[-1][0])
                #Store these maxima values with the peak 
                peaks[peakID]["UV"] = UV_maxima


        #Append the peaks to the specific well
        masterTable[wellno] = peaks.values()

    return [masterTable, chroma, sample_IDs]