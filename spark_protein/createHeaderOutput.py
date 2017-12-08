import math
import numpy

def createHeaderOutput(FileName,short_name,numFile,nfeature ,listVarMatch,SelShowVar):

    len1 = numFile*nfeature
    allVar = listVarMatch+SelShowVar

    header1 = []
    header2 = []

    for i in range(1,len1+1):
        ind = i%nfeature
        if ind==0:
            ind = nfeature
        if ind==1:
            n = int(math.ceil(i/nfeature))
            filename2 = short_name[n-1]
            header2.append(filename2)
        else:
            header2.append('')

        header1.append(allVar[ind-1])
    header1_1 = numpy.asarray(header1)
    header2_1 = numpy.asarray(header2)

    header = numpy.vstack((header2_1,header1_1))
    return header

    tt=1



