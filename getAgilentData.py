"""
The .rpt file provided by Waters instrument is a text file
containing the output from each detector (UV, ELSD, ES+ and ES-),
described as "functions" sequentially. The functions are described for 
each well in turn. 

The following script extracts the data relevant to PyParse, and reformats
it into a common structure. It finds the appropriate sections of text by 
looking for keywords that describe the start of a particular section. 

"""

import math
import logging
import sys

def init(args):
    global options
    options = args

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
    
    #Normalise the heights of the peaks so that they plot correctly

    max_y_val = max(y_val)
    normalised_y_val = [100*i/max_y_val for i in y_val]

    return [x_val, normalised_y_val]
           
        
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
    #Remove any masses which, as a percentage, round to 0 
    #to remove unnecessary baseline ions
    refined_masses = []
    for i in masses:
        if math.floor((i[1]/total)*100) > 0:
            refined_masses.append([i[0], i[1]])
    
    #Get ioncount, in rare cases where the instrument assigns two mass spectrum to the same peakID
    ioncount = float(spectrum.split("TIC")[1].split("\n")[0].strip())

    return [refined_masses, ioncount]
   
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

def getData(filename):
    """
    Converts the input .rpt file into a dictionary where each index corresponds to a well.
    Each dictionary value contains a list of peaks, where each peak is represented as a dictionary
    with the peak area, m/z and UV (if available) data present. 
    
    :param filename: Address of the input file
    
    :return: List comprising: a dictionary of all peaks in all wells; a dictionary containing each chromatogram;
        a dictionary of the sample IDs for each well; a dictionary of the total_abs_area of each well.
        Each of these dictionaries are indexed by well. 
        
    """
    
    wellData = []
    masterTable = {}
    chroma = {}
    sample_IDs = {}
    total_area_abs = {}

    with open(filename, errors = "ignore") as f:
        fullText = f.read()
        text_split = fullText.split("[SAMPLE]") #Split the file into individual wells
        wellData = text_split[1:] 
        

    def findDataEntry(peaks, peakID):
        """
        The program should match up equivalent peaks from different 
        detectors. Detectors may have subtly different retention times, but 
        the asr file tags equivalent peaks with the same peakID. Use this value 
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
            
        
        
        try:
            #Find from the "Well" header for this sample which well it is in
            welltext = functions[0].split("Well")[1].split("-")[1].split("\n")[0].strip() 
            col_as_int = int(welltext[1:])
            row_as_int = ord(welltext[0].capitalize())-64

            #Calculate the wellno as a single integer
            wellno = (row_as_int - 1) * options.plate_col_no + col_as_int

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
            if "Type\tDAD" in lines[4] and options.detector == "UV":
                chromatograms = function.split("[CHROMATOGRAM]")[1:]
                for chromatogram in chromatograms:
                    c_lines = chromatogram.split("\n")
                    #Use just the 254nm wavelength chromatogram rather than a specific wavelength
                    if "Description	* DAD1C,Sig=254,4  Ref=off [s]" in c_lines[3]:
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
                                        "peakID": peakID,
                                        "well": wellno,
                                        "MS+ioncount": 0,
                                        "MS-ioncount": 0,

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
                                        "well": wellno,
                                        "MS+ioncount": peaks[peak_index]["MS+ioncount"],
                                        "MS-ioncount": peaks[peak_index]["MS-ioncount"],

                                        }
                                #Delete the old entry now the data has been copied over
                                #to use the correct retention time. 
                                del peaks[peak_index]
                                
                            peaks[retTime]["area"] = peakArea
                            peaks[retTime]["areaAbs"] = peakAreaAbs
                            peaks[retTime]["pStart"] = float(peakWidth[0])
                            peaks[retTime]["pEnd"] = float(peakWidth[1])
                            
                            if welltext == "C2":
                                print(welltext, retTime, peaks[retTime]["peakID"], peaks[retTime]["areaAbs"])
                                



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
                            "peakID": peakID,
                            "well": wellno,
                            "MS+ioncount": 0,
                            "MS-ioncount": 0,
                            }
                else:
                    retTime = peak_index
                    
                #For each function, parse all the associated data and store
                #as a float. 
                #For UV, take only the lambaMax
                if "IonMode\tEsi" in lines[3] and "Polarity\tPositive" in lines[5]:
                    [data, ioncount] = getMSData(spectrum)
                    if ioncount > peaks[retTime]["MS+ioncount"]:
                        peaks[retTime]["MS+"] = data
                        peaks[retTime]["MS+ioncount"] = ioncount
                    if welltext == "C2":
                        print(welltext, retTime, peaks[retTime]["peakID"], peaks[retTime]["areaAbs"])
                        print(data)

                      
                if "IonMode\tEsi" in lines[3] and "Polarity\tNegative" in lines[5]:
                    [data,ioncount] = getMSData(spectrum)
                    if ioncount > peaks[retTime]["MS-ioncount"]:
                        peaks[retTime]["MS-"] = data
                        peaks[retTime]["MS-"] = ioncount
                    
                
                #Placeholder if a sample Agilent ASR file with UV absorbance data
                #for each peak is ever obtained. 
                #elif lines[4] == "Type\tDAD ":
                #    data = getUVData(spectrum)
                #    peaks[retTime]["UV"] = data
            
        logging.debug(f'{len(peaks)} peaks found in well {wellno}.')             
    
        #Filter out any peaks within retention times specified by the user
        #Expected uses: remove solvent peaks from consideration
        #Expected format is a list of specified ranges, e.g.
        #[min_time1-max_time1, min_time2-max_time2]
        #where the start of the range appears first and is separated from 
        #the end of the range with a hyphen

        if len(options.filter_by_rt) > 0:
            to_remove = []
            for filter_range in options.filter_by_rt:
                min_time = float(filter_range.split("-")[0].strip())
                max_time = float(filter_range.split("-")[1].strip())

                for index in peaks.keys():
                    if float(index) > min_time and index < max_time:
                        to_remove.append(index)

            for index in to_remove:
                del peaks[index]
        
        #If the reprocess_by_rt parameter was used, recalculate the percentage peak area for each peak
        if len(options.filter_by_rt) > 0:
            total_area_abs[wellno] = 0
            for peak in peaks.values():
                total_area_abs[wellno] = total_area_abs[wellno] + peak["areaAbs"]

        #Recalculate the peak area percentages now that the total absolute area 
        #for that well has been calculated. 

        if len(options.filter_by_rt) > 0: 
            for index in peaks.keys():
                peaks[index]["area"] = round((peaks[index]["areaAbs"]*100) / total_area_abs[wellno], 2)
        
        masterTable[wellno] = peaks.values()  

    return [masterTable, chroma, sample_IDs, total_area_abs]