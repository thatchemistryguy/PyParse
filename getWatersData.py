#External Libraries
import math
import pandas as pd

#Internal Libraries
import plate_tools

class rawWatersData:

    """
    Processes the raw .rpt Waters datafile and creates multiple Pandas Dataframes that will facilitate easy 
    access to all the data contained within. 

    The following are descriptions of the Pandas dataframes that will be generated for this object:

    self.rawDADTable = pd.DataFrame(columns =['filename', 'well', 'peakID', 'time', 'pStart', 'pEnd', 'area', 'areaAbs'])
    self.rawUVTable = pd.DataFrame(columns =['filename', 'well', 'peakID', 'time', 'UVvalue', 'UVintensity'])
    self.rawUVMaximaTable = pd.DataFrame(columns =['filename', 'well', 'peakID', 'time', 'UVvalue'])
    self.rawELSDTable = pd.DataFrame(columns =['filename', 'well', 'peakID', 'time', 'pStart', 'pEnd', 'area', 'areaAbs'])
    self.rawMSTable = pd.DataFrame(columns =['filename', 'well', 'peakID', 'time', 'MSvalue', 'MSintensity', 'MStype', 'perc_intensity'])
    self.rawTraceTable = pd.DataFrame(columns =['filename', 'well', 'time', 'height', 'detector_type'])

    """

    def __init__(self, inputfile, row_no = 0, col_no = 0):
        """
        Initialises the rawWatersData object. Requires the inputfile to be set. 

        :param inputfile: a string describing the file location of the rpt file. 
        :param row_no: (optional) an integer describing the number of rows in the plate, as specified by the user
        :param col_no: (optional) an integer describing the number of columns in the plate, as specified by the user
        
        :return: A python object. 
        """

        self.row_no = row_no
        self.col_no = col_no
        self.username = ""
        self.sample_IDs = {}
        self.method_IDs = {}

        #Open the rpt file, and split the filetext into individual wells
        with open(inputfile, errors = "ignore") as f:
            fullText = f.read()
            self.wellData = fullText.split("[SAMPLE]")[1:] 
    
    def getWellFormat(self):

        """
        Get the plate dimensions from the first sample
        Data sample: #Plate	01TL,XY,SD,1: 8,2:12,3: 90.0...
        where number of rows is 8 and number of columns in 12
        Overwrite the default option, but ensure that the user retains control
        such that empty rows can be removed from the heatmap.

        :return: Sets the number of rows/columns (if the user didn't do this), stores the row/col order,
        stores how each well is described. 
        """

        
        if self.row_no == 0 or self.col_no == 0:
            self.row_no = int(self.wellData[0].split("\n")[17].split(",")[3].split(":")[1])
            self.col_no = int(self.wellData[0].split("\n")[17].split(",")[4].split(":")[1])
            self.plate_cols_for_extraction = self.col_no
        else:
            self.plate_cols_for_extraction = int(self.wellData[0].split("\n")[17].split(",")[4].split(":")[1])
            
            
        #Find from rpt file how each well is specified 
        #Data sample: #Plate	01TL,XY,SD,1: 8,2:12,3: 90.0...
        self.row_col_order = self.wellData[0].split("\n")[17].split(",")[1]
        self.well_type = self.wellData[0].split("\n")[17].split(",")[2]

    def getMetaData(self):
        """
        Get the username of the person who submitted the sample and the date it was run. 

        :return: Stores the username and daterun as part of the object. 
        """
        
        self.username = self.wellData[0].split("UserName")[1].split("\n")[0].strip()
        self.daterun = self.wellData[0].split("Date")[1].split("\n")[0].strip()

        
    def getWell(self, position):
        """
        Takes in the string from the .rpt file describing the well and returns
        an integer. Works whether the well is described as a row/col coordinate or single digit. 

        :param position: Integer describing which index of self.wellData is currently being examined. 
        :return: Integer describing the well number of the sample currently being examined. 

        """

        well_number = -1
        #If the well type is just a Single Digit...
        if self.well_type == "SD":
            #If the well is simply an integer between 1 and infinity
            #Single line function to trim full string to just the well number used
            well_number = int(self.wellData[position].split("Well")[1].split("\n")[0].split(":")[1].strip()) 

            #If the column number specified by the user is different to that found in the rpt file, 
            #then we assume that the user looking to trim off blank columns on the right of the plate. Only wells described by 
            #a single digit need to be modified to take this into account. 
            if self.col_no != self.plate_cols_for_extraction:
                well_number = math.floor(well_number / self.plate_cols_for_extraction)*self.col_no + (well_number % self.plate_cols_for_extraction)

        #If the well type is a combination of letters/numbers...
        else:
            #Find if the well is in "XY" format or "YX" format
            if self.row_col_order == "XY":
                column = self.wellData[position].split("Well")[1].split("\n")[0].split(":")[1].split(",")[0].strip()
                row = self.wellData[position].split("Well")[1].split("\n")[0].split(":")[1].split(",")[1].strip()
            else: 
                column = self.wellData[position].split("Well")[1].split("\n")[0].split(":")[1].split(",")[1].strip()
                row = self.wellData[position].split("Well")[1].split("\n")[0].split(":")[1].split(",")[0].strip()

            #Convert the column to integer, either by direct
            #conversion, or by finding position in the alphabet
            try:
                col_as_int = int(column)
            except:
                col_as_int = ord(column.capitalize()) - 64

            #Convert the row to integer, either by direct
            #conversion, or by finding position in the alphabet
            try: 
                row_as_int = int(row)
            except: 
                row_as_int = ord(row.capitalize()) - 64

            #Calculate the wellno as a single integer
            well_number = (row_as_int - 1) * self.plate_cols_for_extraction + col_as_int
        return well_number
    
    def processSamples(self):
        """
        #Get a dataframe containing all the individual sample names

        :return: Pandas dataframe containing sample names. 
        self.rawSampleTable = pd.DataFrame(columns =['filename', 'sampleID', 'methodID', 'well', 'well_readable'])
        """

        
        sample_list = []
        for i in range(len(self.wellData)):
            sample_list.append({
                "filename": self.wellData[i].split("FileName")[1].split("\n")[0].strip(),
                "sampleID": self.wellData[i].split("SampleID")[1].split("\n")[0].strip(),
                "methodID": self.wellData[i].split("Method\t")[1].split("\n")[0].split("\\")[-1].strip(),
                "well": self.getWell(i), 
                "well_readable": plate_tools.getUserReadableWell(self.getWell(i), self.col_no)  
            })

        self.rawSampleTable = pd.DataFrame(sample_list)


    def processDAD(self):
        """
        Get a dataframe containing all the DAD data (peaks and peak areas). 
        This function also adds a list of sample_IDs and method_IDs to rawWatersData object 
        if this isn't already present. 

        :return: Pandas dataframe containing DAD data, indexed by well and peakID.  
        self.rawDADTable = pd.DataFrame(columns =['filename', 'well', 'peakID', 'time', 'pStart', 'pEnd', 'area', 'areaAbs'])
        
        """

        peak_list = []
        for i in range(len(self.wellData)):
            functions = self.wellData[i].split("[FUNCTION]")
            well_number = self.getWell(i)
            
            #Make sure the sample ID has been added for this well. This process is repeated for the other 
            #"process_something" functions, so that the data is available regardless of which detectors 
            #are of interest. 
            if well_number not in self.sample_IDs:
                well_ID = functions[0].split("SampleID")[1].split("\n")[0].strip()  
                self.sample_IDs[well_number] = well_ID

            if well_number not in self.method_IDs:
                method_ID = functions[0].split("Method")[1].split("\n")[0].strip()  
                self.method_IDs[well_number] = method_ID
                
            for j in range(len(functions[1:])):
                function = functions[1:][j]
                lines = function.split("\n")
                filename = functions[0].split("FileName")[1].split("\n")[0].strip()
                
                #if the current section of file is indeed related to the DAD output...
                if "Type\tDAD" in lines[4]:
                    chromatograms = function.split("[CHROMATOGRAM]")[1:]
                    for chromatogram in chromatograms:
                        c_lines = chromatogram.split("\n")
                        if "Description\tDAD:" in c_lines[3]:
                            spectra = function.split("[SPECTRUM]")[1:] #split by spectrum (i.e. each peak)
                            peaks = chromatogram.split("[PEAK]")[1:]
                            for peak in peaks:
                                new_entry = {
                                    "filename": filename,
                                    "well": well_number,
                                    "peakID": int(peak.split("Peak ID")[1].split("\n")[0].strip()),
                                    "time": float(peak.split("Time")[1].split("\n")[0].strip()),
                                    "pStart": float(peak.split("Peak\t")[1].split("\n")[0].split("\t")[0]),
                                    "pEnd": float(peak.split("Peak\t")[1].split("\n")[0].split("\t")[1]),
                                    "area": float(peak.split("Area %Total")[1].split("\n")[0].strip()),
                                    "areaAbs": float(peak.split("AreaAbs")[1].split("\n")[0].strip()),
                                }
                                peak_list.append(new_entry)

        self.rawDADTable = pd.DataFrame(peak_list)

    def processUV(self):
        """
        Extracts the UV wavelength/intensities for all peaks in all samples and creates a Dataframe from them. 
        Note: Not all LCMSs report the UV data, depending on how they're configured, so this dataframe may be empty. 
            
        :return: Pandas Dataframe for all UV wavelength/intensities for each peak in each sample. 
        self.rawUVTable = pd.DataFrame(columns =['filename', 'well', 'peakID', 'time', 'UVvalue', 'UVintensity'])
        """


        peak_list = []
        for i in range(len(self.wellData)):
            functions = self.wellData[i].split("[FUNCTION]")
            well_number = self.getWell(i)
            
            #Make sure the sample ID has been added for this well. This process is repeated for the other 
            #"process_something" functions, so that the data is available regardless of which detectors 
            #are of interest. 
            if well_number not in self.sample_IDs:
                well_ID = functions[0].split("SampleID")[1].split("\n")[0].strip()  
                self.sample_IDs[well_number] = well_ID

            if well_number not in self.method_IDs:
                method_ID = functions[0].split("Method")[1].split("\n")[0].strip()  
                self.method_IDs[well_number] = method_ID
                
            for j in range(len(functions[1:])):
                function = functions[1:][j]
                lines = function.split("\n")
                filename = functions[0].split("FileName")[1].split("\n")[0].strip()
                if "Type\tDAD" in lines[4]:
                    spectra = function.split("[SPECTRUM]")[1:] #split by spectrum (i.e. each peak)

                    for spectrum in spectra:
                        #Grab the list of intensities for each wavelength observed by the detector
                        #for each peak in each sample. 
                        lineData = spectrum.split(";Mass\t% BPI")[1].split("\n")
                        UVdata = []
                        for line in lineData:
                            
                            if line == "}":#stop for loop if end of UVData section is reached
                                break
                            
                            UVdata = line.split("\t")
                            if len(UVdata) == 2:
                                UVdata.append([float(UVdata[0]), abs(float(UVdata[1]))])

                        #Load this data into a form ready for conversion to DataFrame. 
                        for data in UVdata:
                            new_entry = {
                                "filename": filename,
                                "well": well_number,
                                "peakID": int(spectrum.split("Peak ID")[1].split("\n")[0].strip()),
                                "time": float(spectrum.split("Time")[1].split("\n")[0].strip()),
                                "UVvalue": data[0],
                                "UVintensity": data[1]
                            }
                            peak_list.append(new_entry)

        if len(peak_list) == 0:
            self.rawUVTable = pd.DataFrame(columns =['filename', 'well', 'peakID', 'time', 'UVvalue', 'UVintensity'])
        else:
            self.rawUVTable = pd.DataFrame(peak_list)    
        
    def processUVmaxima(self, min_uv_threshold = 20): 
        """
            Takes the specific region of the rpt file pertaining to 
            UV absorbance spectrum data for a specific peak in a specific
            well, and returns a list containing all the maxima of that spectrum.
            The height of the maxima must be greater than the min_uv_threshold
            specified in options.

            Note: Not all LCMSs report the UV data, depending on how they're configured, so this dataframe may be empty. 
            
            :param min_uv_threshold: (optional) integer describing how intense absorbance at a
            particular wavelength must be before it's classed as a maxima
            
            :return: DataFrame of UV maxima
            self.rawUVMaximaTable = pd.DataFrame(columns =['filename', 'well', 'peakID', 'time', 'UVvalue'])
            """

        peak_list = []
        for i in range(len(self.wellData)):
            functions = self.wellData[i].split("[FUNCTION]")
            well_number = self.getWell(i)
            
            #Make sure the sample ID has been added for this well. This process is repeated for the other 
            #"process_something" functions, so that the data is available regardless of which detectors 
            #are of interest. 
            if well_number not in self.sample_IDs:
                well_ID = functions[0].split("SampleID")[1].split("\n")[0].strip()  
                self.sample_IDs[well_number] = well_ID

            if well_number not in self.method_IDs:
                method_ID = functions[0].split("Method")[1].split("\n")[0].strip()  
                self.method_IDs[well_number] = method_ID
                
            for j in range(len(functions[1:])):
                function = functions[1:][j]
                lines = function.split("\n")
                filename = functions[0].split("FileName")[1].split("\n")[0].strip()

                if "Type\tDAD" in lines[4]:
                    spectra = function.split("[SPECTRUM]")[1:] #split by spectrum (i.e. each peak)

                    for spectrum in spectra:

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
                        if UVy[0] > UVy[1] and UVy[0] > min_uv_threshold:
                            UVmaxima.append(UVx[0])
                        for i in range(1, len(UVy)-1):
                            if UVy[i] > UVy[i-1] and UVy[i] > UVy[i+1] and UVy[i] > min_uv_threshold:
                                UVmaxima.append(UVx[i])
                        if UVy[-1] > UVy[-2] and UVy[-1] > min_uv_threshold:
                            UVmaxima.append(UVx[-1])   

                        #process the data into a form ready for conversion to DataFrame
                        for maxima in UVmaxima:
                            new_entry = {
                                "filename": filename,
                                "well": well_number,
                                "peakID": int(spectrum.split("Peak ID")[1].split("\n")[0].strip()),
                                "time": float(spectrum.split("Time")[1].split("\n")[0].strip()),
                                "UVvalue": maxima,
                            }
                            peak_list.append(new_entry)

        if len(peak_list) == 0:
            self.rawUVMaximaTable = pd.DataFrame(columns =['filename', 'well', 'peakID', 'time', 'UVvalue'])
        else:
            self.rawUVMaximaTable = pd.DataFrame(peak_list)
    
    def processELSD(self):
        """
        Get a dataframe containing all the ELSD data (peaks and peak areas). 
        This function also adds a list of sample_IDs and method_IDs to rawWatersData object 
        if this isn't already present. 

        Note: Not all LCMSs have an ELSD detector, so this may return an empty dataframe. 

        :return: Pandas dataframe containing ELSD data, indexed by well and peakID.  
        self.rawELSDTable = pd.DataFrame(columns =['filename', 'well', 'peakID', 'time', 'pStart', 'pEnd', 'area', 'areaAbs'])
        
        """
        peak_list = []
        for i in range(len(self.wellData)):
            functions = self.wellData[i].split("[FUNCTION]")
            well_number = self.getWell(i)
            
            #Make sure the sample ID has been added for this well. This process is repeated for the other 
            #"process_something" functions, so that the data is available regardless of which detectors 
            #are of interest. 
            if well_number not in self.sample_IDs:
                well_ID = functions[0].split("SampleID")[1].split("\n")[0].strip()  
                self.sample_IDs[well_number] = well_ID

            if well_number not in self.method_IDs:
                method_ID = functions[0].split("Method")[1].split("\n")[0].strip()  
                self.method_IDs[well_number] = method_ID
                
            for j in range(len(functions[1:])):
                function = functions[1:][j]
                lines = function.split("\n")
                filename = functions[0].split("FileName")[1].split("\n")[0].strip()
                #get peakarea for this peak
                if "Description\tANALOG" in lines[3]:
                    chromatograms = function.split("[CHROMATOGRAM]")[1:]
                    for chromatogram in chromatograms:
                        c_lines = chromatogram.split("\n")
                        if "ELSD" in c_lines[3]:
                            spectra = function.split("[SPECTRUM]")[1:] #split by spectrum (i.e. each peak)
                            #chroma[wellno] = getChromatogram(chromatogram.split("[TRACE]")[1]) #get chromatogram for this well
                            peaks = chromatogram.split("[PEAK]")[1:]
                            for peak in peaks:
                                new_entry = {
                                    "filename": filename,
                                    "well": well_number,
                                    "peakID": int(peak.split("Peak ID")[1].split("\n")[0].strip()),
                                    "time": float(peak.split("Time")[1].split("\n")[0].strip()),
                                    "pStart": peak.split("Peak\t")[1].split("\n")[0].split("\t")[0],
                                    "pEnd": peak.split("Peak\t")[1].split("\n")[0].split("\t")[1],
                                    "area": float(peak.split("Area %Total")[1].split("\n")[0].strip()),
                                    "areaAbs": float(peak.split("AreaAbs")[1].split("\n")[0].strip())
                                }
                                peak_list.append(new_entry)
        if len(peak_list) == 0:
            self.rawELSDTable = pd.DataFrame(columns =['filename', 'well', 'peakID', 'time', 'pStart', 'pEnd', 'area', 'areaAbs'])
        else:
            self.rawELSDTable = pd.DataFrame(peak_list)
        
    def processMS(self):
        """
        Get a dataframe containing all the MS data (MS values, intensities and whether its MS+ or MS-). 
        This function also adds a list of sample_IDs and method_IDs to rawWatersData object 
        if this isn't already present. 

        :return: Pandas dataframe containing MS data, indexed by well and peakID.  
        self.rawMSTable = pd.DataFrame(columns =['filename', 'well', 'peakID', 'time', 'MSvalue', 'MSintensity', 'MStype', 'perc_intensity'])
        
        """

        peak_list = []
        for i in range(len(self.wellData)):
            functions = self.wellData[i].split("[FUNCTION]")
            well_number = self.getWell(i)
            
            #Make sure the sample ID has been added for this well. This process is repeated for the other 
            #"process_something" functions, so that the data is available regardless of which detectors 
            #are of interest. 
            if well_number not in self.sample_IDs:
                well_ID = functions[0].split("SampleID")[1].split("\n")[0].strip()  
                self.sample_IDs[well_number] = well_ID
            
            if well_number not in self.method_IDs:
                method_ID = functions[0].split("Method")[1].split("\n")[0].strip()  
                self.method_IDs[well_number] = method_ID

            filename = functions[0].split("FileName")[1].split("\n")[0].strip()    
            for j in range(len(functions[1:])):
                function = functions[1:][j]
                lines = function.split("\n")
                if "IonMode\tES" in lines[3]:
                    
                    spectra = function.split("[SPECTRUM]")[1:] #split by spectrum (i.e. each peak)
                    for spectrum in spectra:
                        if "IonMode\tES+" in lines[3]:
                            MStype = "+"
                        elif "IonMode\tES-" in lines[3]:
                            MStype = "-"

                        #get the list of mass ions observed and their intensities
                        masses = []
                        lineData = spectrum.split(";Mass\t% BPI")[1].split("\n")

                        for line in lineData[1:]:
                            if line == "}": #stop the for loop if end of MSData section is reached
                                break
                        
                            massData = line.split("\t")
                            if len(massData) == 2:
                                floatData = [float(i) for i in massData] #convert all data to float
                                masses.append(floatData)
                                
                        #Remove any masses which round to 0 
                        #to remove unnecessary baseline ions
                        MSdata = []
                        for i in masses:
                            if math.floor(i[1]) > 0:
                                MSdata.append([i[0], i[1]])
                        
                        for ion in MSdata:
                            new_entry = {
                                "filename": filename,
                                "well": well_number,
                                "peakID": int(spectrum.split("Peak ID")[1].split("\n")[0].strip()),
                                "time": float(spectrum.split("Time")[1].split("\n")[0].strip()),
                                "MSvalue": float(ion[0]),
                                "MSintensity": float(ion[1]),
                                "MStype": MStype
                            }
                            peak_list.append(new_entry)

        
        self.rawMSTable = pd.DataFrame(peak_list)
        
        #Calculate the sum of the MSintensities in each peak, then calculate the of each MSintensity to this total
        total_intensities = self.rawMSTable.groupby(["well", "peakID", "MStype"]).agg(total_intensity=('MSintensity', 'sum'))
        self.rawMSTable = self.rawMSTable.join(total_intensities, on=["well", "peakID", "MStype"], rsuffix = "right")
        self.rawMSTable["perc_intensity"] = 100 * self.rawMSTable["MSintensity"] / self.rawMSTable["total_intensity"]
        
    def processTrace(self, points_per_trace = 500):
        
        """
        Takes the region of text in rpt file corresponding to the chromatogram,
        and creates a DataFrame from this information. Contains info from both UV and ELSD detectors. 

        :param points_per_trace: (optional) Integer describing the maximum number of datapoints that 
        should be extracted from the chromatogram. 

        :return: Pandas Dataframe for all datapoints in all chromatograms
        self.rawTraceTable = pd.DataFrame(columns =['filename', 'well', 'time', 'height', 'detector_type'])
        """            
        trace_list = []
        
        for i in range(len(self.wellData)):
            functions = self.wellData[i].split("[FUNCTION]")
            well_number = self.getWell(i)
            
            for function in functions[1:]:
                lines = function.split("\n")
                if "Type\tDAD" in lines[4] or "Description\tANALOG" in lines[3]:
                    chromatograms = function.split("[CHROMATOGRAM]")[1:]
                    filename = functions[0].split("FileName")[1].split("\n")[0].strip()
                    for chromatogram in chromatograms:
                        c_lines = chromatogram.split("\n")
                        
                        #determine if this chromatogram is DAD or ELSD data
                        #note that it depends on how the detector was named during LCMS
                        #installation, so this line may need to be changed for different institutions
                        detector_type = "ELSD" if "ELSD" in c_lines[3] else "DAD"
                        trace_text = chromatogram.split("[TRACE]")[1].split("\n")

                        length = len(trace_text)
                        n_val = math.ceil(length / points_per_trace)

                        for index, value in enumerate(trace_text[1:]):
                            if value == "}": #stop for loop if end of data section is reached
                                break
                            elif index % n_val == 0:

                                if value != "{":
                                    data = value.split("\t")
                                    if len(data) == 2:
                                        goingin = {
                                            "filename": filename,
                                            "well": well_number,
                                            "time": float(data[0]),
                                            "height": float(data[1]), 
                                            "detector_type": detector_type
                                        }
                                        trace_list.append(goingin)

        self.rawTraceTable = pd.DataFrame(trace_list)
        
        
            