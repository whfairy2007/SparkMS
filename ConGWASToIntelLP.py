import csv
import os
import time
import threading
from os import listdir
from os.path import isfile, join

class MyUberMultithreadingWarMachine (threading.Thread):
    def __init__(self, threadID, path, filesName, caseOrControl):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.path = path
        self.filesName = filesName
        self.caseOrControl = caseOrControl
    def run(self):
        print("Starting " + self.name)
        modifyCSV(self.name, self.path, self.filesName, self.caseOrControl)
        print("Exiting  " + self.name)

def modifyCSV(threadName, path, filesName, caseOrControl):
    start_time = time.time()
    for fileName in filesName:
        fileLocation = path + fileName
        with open(fileLocation, newline='') as csvfile:
            CSVReader = csv.reader(csvfile, delimiter=',')
            header = next(CSVReader) # I do this to remove the header from the files
            print(header)
            value = "1.0 0.0" if (caseOrControl == "case") else "0.0 1.0"
            weight = "1"

            # Allocating only once the name of the curated file by extracting the first line with payload
            firstLine = next(CSVReader)
            source = firstLine[1]
            destination = firstLine[9] + firstLine[10] + "_" + firstLine[4]
            newRow = [source, value, destination, weight]
            newFileName = path + source + "_curated.csv"
            writeCSV(newFileName, newRow)
            i = 0

            # Now all the other lines
            for row in CSVReader:
                if(i==1000):
                    break
                if (float(row[2]) != 0):
                    #                    source = row[1]
                    i += 1
                    destination = row[9] + row[10] + "_" + row[4]
                    newRow = [source, value, destination, weight]
                    #                    newFileName = path + source + "_curated.csv"
                    writeCSV(newFileName, newRow)
    elapsed_time = time.time() - start_time
    print("the thread %s is done and it took: %f" %(threadName,elapsed_time))
    return 0

def writeCSV(newFileName, newRow):
    writeCSVDirectEdge(newFileName, newRow)
    writeCSVReverseEdge(newFileName, newRow)
    return 0


def writeCSVDirectEdge(newFile, newRow):
    with open(newFile, 'a', newline='') as csvfile:
        lpgwas = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        lpgwas.writerow(newRow)
    return 0

def writeCSVReverseEdge(newFile, newRow):
    newValue = "0.5 0.5"
    reversedEdgeRow = [newRow[2], newValue, newRow[0], newRow[3]]
    with open(newFile, 'a', newline='') as csvfile:
        lpgwas = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        lpgwas.writerow(reversedEdgeRow)
    return 0

def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i+n]

#####################################
# main program                      #
#####################################

#fileDirectory = "C:\\where\\is\\your\\project\ConverterGWAStoIntelLP\\csv_test_files\\"
#fileDirectory = "C:\\Users\\Axel\\Work\\ConverterGWAStoIntelLP\\csv_test_files\\"
fileDirectory = "/home/impuser/Desktop/DataSets/control_temp/"
numberPerChunks = 11
caseOrControl = "control"
if os.path.exists(fileDirectory):
    filesName = [f for f in listdir(fileDirectory) if isfile(join(fileDirectory, f))]
    chunks = chunks(filesName, numberPerChunks)
    chunksList = list(chunks)
    for i in range(len(chunksList)):
        warMachine = MyUberMultithreadingWarMachine(i, fileDirectory, chunksList[i], caseOrControl)
        warMachine.start()
else:
    print("This directory doesn't exist")
