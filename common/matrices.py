#! /usr/bin/python

# matrices.py
# Author: Santhosh Balasubramanian
# Created: March 24, 2013
# Last Modified: March 24, 2013


# - - - - - - - - - - MATRIX CREATION - - - - - - - - - - #


# constructSquareMatrix: Constructs empty square matrix of length <dimension>.
def constructSquareMatrix(dimension):
    # Initialize Rows
    matrix = dimension*[None]
    
    # Fill Columns
    for i in range(dimension):
        matrix[i] = dimension*[None]

    return matrix

# constructSquareIntegerMatrix: Constructs empty integer square matrix.
def constructSquareIntegerMatrix(dimension):
    # Initialize Rows
    matrix = dimension*[0]
    
    # Fill Columns
    for i in range(dimension):
        matrix[i] = dimension*[0]

    return matrix


# - - - - - - - - - - MATRIX I/O - - - - - - - - - - #


# writeMatrixToFile: Writes <matrix> to <fileHandle>.
def writeMatrixToFile(matrix, outFile):

    # Iterate through Rows
    for row in matrix:
        # Iterate through Columns
        numColumns = len(row)
        for i in range(numColumns):
            column = row[i]
            if i == 0:
                outFile.write('%s' % str(column))
            else:
                outFile.write('\t%s' % str(column))
        outFile.write('\n')
    outFile.write('\n')
    
    return

