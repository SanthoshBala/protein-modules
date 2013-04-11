#! /usr/bin/python

# paralellization.py
# Author: Santhosh Balasubramanian
# Created: March 24, 2013
# Last Modified: March 24, 2013


# Python Imports
import math
from multiprocessing import Process

# Global Imports
from settings import *


# - - - - - - - - - - LOCAL SETTINGS - - - - - - - - - - #


NUM_CPU_CORES = 4


# - - - - - - - - - - TISSUE - - - - - - - - - - #


# parallelizeTaskByRecord: Paralellizes <functionHandle> amongst
# 2*<NUM_CPU_CORES> processes by splitting along DB record.
def parallelizeTaskByRecord(functionHandle, numRecords):
    # Determine # Processes and Records Per Process
    numProcesses = 2*NUM_CPU_CORES
    recordsPerProcess = int(math.ceil(numRecords/float(numProcesses)))
    processList = []

    # Create Processes
    for i in range(numProcesses):
        processList.append( Process( target = functionHandle, 
                                     args = ( i*recordsPerProcess,
                                              recordsPerProcess ) ) )

    # Start Processes
    for i in range(numProcesses):
        processList[i].start()

    # Join Processes
    for i in range(numProcesses):
        processList[i].join()

    return

# parallelizeTaskByItem: Parallelizes <functionHandle> by
# dividing a parameter amongst N buckets and setting as arguments.
def parallelizeTaskByItem(functionHandle, numItems):
    # Determine # Processes and # Tissues Per Process
    numProcesses = 2*NUM_CPU_CORES
    itemsPerProcess = int(math.ceil(numItems/float(numProcesses)))

    # Create Processes
    processList = []
    for i in range(numProcesses):
        processList.append( Process(target=functionHandle, 
                                    args = (i*itemsPerProcess, 
                                            (i+1)*itemsPerProcess)) )

    # Start Processes
    for i in range(numProcesses):
        processList[i].start()

    # Join Processes
    for i in range(numProcesses):
        processList[i].join()

    return


# parallelizeTaskByTissue: Paralellizes <functionHandle> amongst 
# 2*<NUM_CPU_CORES> processes. Joins processes before returning.
def parallelizeTaskByTissue(functionHandle, augmented = False):

    if augmented:
        numItems = len(AUGMENTED_TISSUE_LIST)
    else:
        numItems = len(FUNCTIONAL_TISSUE_LIST)

    parallelizeTaskByItem(functionHandle, numItems)

    return

