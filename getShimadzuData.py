import xml.etree.ElementTree as ET
import base64
import glob

def getData(input_dir):
    wellData = []
    masterTable = {}
    chroma = {}
    sample_IDs = {}

    input_files = glob.glob(input_dir + '/*.daml')
    #print(input_files)
    for filename in input_files:
        #Open the .daml file, which is structured as an xml
        tree = ET.parse(filename)
        root = tree.getroot()
        wellno = int(root.find("datafile_attribs").find("vial").text)
        sample_IDs[wellno] = str(root.find("datafile_attribs").find("desc").text)
        #print(wellno)
        #print(root.tag)
        results = root.find("results")

        peaks = {}
        #start by getting the peaks observed in the UV spectrum
        raw_UV = results.find("lc_chros").find("chro")
        #print(raw_UV)

        for peak in raw_UV.findall("peak"):

            new_entry = {}
            for detail in peak.get("details").split(","):
                name = detail.strip().split("=")[0]
                value = detail.strip().split("=")[1]
                if name == "Peak#":
                    new_entry["peakID"] = value
                elif name == "RT":
                    new_entry["time"] = float(value.split("mins")[0].strip())
                    new_entry["pStart"] = float(value.split("(")[1].split("-")[0])
                    new_entry["pEnd"] = float(value.split("-")[1].split("mins")[0])
                elif name == "Area":
                    new_entry["areaAbs"] = float(value)
                elif name == "Area%":
                    new_entry["area"] = float(value)
                    #Add a placeholder entry for UV data
                    new_entry["UV"] = []
            peaks[new_entry["peakID"]] = new_entry

        #get the chromatographic trace for this well

        coded_trace = raw_UV.find("data").text
        decoded_trace = base64.b64decode(coded_trace)
        trace_list = str(decoded_trace).split(";;;")[1].split("'")[0].split(";")

        xval = []
        yval = []
        for i in range(0, len(trace_list), 2):
            xval.append(float(trace_list[i])/60000)
            yval.append(float(trace_list[i+1])/10000)
        print(filename)
        print(xval)
        print(" ")
        print(yval)
        chroma[wellno] = [xval, yval]

        #next, get the mass spec data for each peak, and match it to the UV peaks

        raw_ms = results.find("spectra")
        for spec in raw_ms.findall("spec"):
            peakID = spec.find("peak_num").text
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
                    total = total + float(decoded_list[i+1])
                    masses.append(datapair)
                
                #normalise the intensities 
                refined_masses = [[i[0], round(i[1]*100/total, 3)] for i in masses]

                if ms_type == "MS+":
                    peaks[peakID]["MS+"] = refined_masses
                else:
                    peaks[peakID]["MS-"] = refined_masses

        



        masterTable[wellno] = peaks.values()

    return [masterTable, chroma, sample_IDs]