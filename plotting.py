#External Libaries
import math


def converttoHeatMapData(outputTable, plate_row_no, plate_col_no):

    """
    Takes the outputTable DataFrame (see genOutput.py) and reformats the 
    data so that it's in a format that's ready for Plotly to conver into a heatmap. 

    :param outputTable: Pandas DataFrame
    :param plate_row_no: Integer, describing the number of rows in the plate
    :param plate_col_no: Integer, describing the number of columns in the plate

    :return: Python dictionary, where each entry is a dictionary containing lists and strings. 
    """

    #Define all the different types of heatmap that are required. 
    #Each entry must correspond to a value in the outputTable. 
    hm_types = {
        "SMarea": "SMarea",
        "Parea": "Parea", 
        "conversion": "P/SM+P", 
        "ratio_to_IS": "P/STD",
        "corrSMarea": "corrSMarea",
        "corrParea": "corrParea", 
        "corrected_conversion": "corrP/SM+P",
        "corrected_ratio_to_IS": "corrP/STD"
    }

    output = {}
    for key, datatype in hm_types.items():
        x_values = []
        y_values = []
        z_values = []

        for i in range(plate_row_no):
            y_values.append(chr(i + 65))
            z_values.append([])
        for i in range(plate_col_no):
            x_values.append(i + 1)

        #Reverse the y-values so that "A" appears at the top of the plotly heatmap
        y_values.reverse()

        for index, row in outputTable.iterrows():
            rowVal = int(math.floor((index-1) / plate_col_no))

            #Plotly requires the data to be in bottom-up
            #format, so reverse the rowVal

            reverse_rowVal = abs(rowVal - plate_row_no + 1)

            if row["SampleID"] == "No Data.":
                z_values[reverse_rowVal].append(0)
            else: 
                if row[datatype] > 1:
                    z_values[reverse_rowVal].append(math.floor(row[datatype]))
                elif row[datatype] >= 100:
                    z_values[reverse_rowVal].append(100)
                elif row[datatype] == 0:
                    z_values[reverse_rowVal].append(0)
                else:
                    z_values[reverse_rowVal].append(round(row[datatype],2))
        
        
        goingin = {
            "name": datatype,
            "x_values": x_values,
            "y_values": y_values, 
            "z_values": z_values,
        }
        output[datatype] = goingin

    return output
