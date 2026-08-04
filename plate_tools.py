#External Libraries
import math
import sys
import logging

def getUserReadableWell(wellno, plate_col_no):
    """
    Converts the well as a number into a user-friendly string,
    e.g. well 11 becomes "B5" for a 4*6 well plate. 

    :param wellno: An integer representing a specific well on the plate
    :param plate_col_no: Integer, the number of columns in the plate
    
    :return: A string representing a specific well on the plate
    """
    
    rowVal = math.floor((wellno-1) / plate_col_no)
    colVal = (wellno) % plate_col_no
    if colVal == 0:
        colVal = plate_col_no
    
    label = f'{chr(ord("@")+(rowVal)+1)}{colVal}'
    return label

def convertWellToNum(well, plate_col_no):
    """
    Converts a human readable well label (like A1) to
    an integer. 

    :param well: string, representing a well in a plate
    :param plate_col_no: integer, number of columns in a plate. 
    """
    row = well[0]
    column = well[1:]
    #if the format of the row/column conforms to expectations
    if isinstance(int(column[0]), int):
        return int((ord(row) - 65) * plate_col_no + int(column))
    #Plates with more than 26 rows are unsupported at present. 
    else:
        logging.info("The well specified implies an unsupported plate.")
        sys.exit(2)