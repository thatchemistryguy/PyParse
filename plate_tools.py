import math

def getUserReadableWell(wellno, plate_col_no):
    """
    Converts the well as a number into a user-friendly string,
    e.g. well 11 becomes "B5" for a 4*6 well plate. 

    :param wellno: An integer representing a specific well on the plate
    
    :return: A string representing a specific well on the plate
    """
    
    rowVal = math.floor((wellno-1) / plate_col_no)
    colVal = (wellno) % plate_col_no
    if colVal == 0:
        colVal = plate_col_no
    
    label = f'{chr(ord("@")+(rowVal)+1)}{colVal}'
    return label