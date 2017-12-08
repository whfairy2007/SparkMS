import numpy

def assignValue1(nfile,cellsz2_x,nfeature,ion_sel):
    alignData_i_j = numpy.zeros((1,nfeature*nfile))
    for lp0 in range(nfile):
        if cellsz2_x[lp0]>0:
            indx = range(lp0*nfeature,(lp0+1)*nfeature)
            alignData_i_j[0,indx] = ion_sel[lp0]

    return alignData_i_j
