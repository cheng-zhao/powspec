#cython: language_level=3
#cython: boundscheck = True
import cython
from cython.parallel import prange, threadid
cimport openmp
from libc.stdlib cimport malloc, free, calloc
from libc.stdio cimport FILE, fprintf, fopen, fclose, printf, fflush, stdout, stderr
import numpy as np

# read_cata.h

cdef extern from "read_cata.h":
    ctypedef struct DATA:
        double x[3]
        double w
    ctypedef struct CATA:
        int num;              # number of data catalogs to be read       
        DATA **data;          # structure for the data catalogs          
        DATA **rand;          # structure for the random catalogs        
        size_t *ndata;        # number of objects in the data catalogs   
        size_t *nrand;        # number of objects in the random catalogs 
        double *wdata;        # total completeness weight of the data    
        double *wrand;        # total completeness weight of the random  
        double *alpha;        # alpha value for each data-random sample  
        double *shot;         # shot noise term for each sample          
        double *norm;         # normalisation term for each sample       
    CATA *cata_init(int ndata) nogil
    void cata_destroy(CATA *cat) nogil


# powspec.h
cdef extern from "powspec.h":

    PK *compute_pk(CATA* cata, bint save_out, int argc, char* argv[]) nogil;

# multipole.h

cdef extern from "multipole.h":
    ctypedef struct PK:
        bint issim;           # indicate if the inputs are from simulations  
        bint log;             # indicate if using logarithm wave number bins 
        bint isauto[2];       # indicate if auto power spectra are required  
        bint iscross;         # indicate if cross power spectrum is required 
        int nl;               # number of multipoles to be computed          
        int nbin;             # number of k bins for the power spectra       
        int nmu;              # number of mu bins for the 2-D power spectra  
        int *poles;           # multipoles to be computed                    
        double los[3];        # line-of-sight vector for simulation boxes    
        double dk;            # width of the wave number bins                
        double *kedge;        # edges of the pre-defined wave number bins    
        double *k;            # centre of the pre-defined wave number bins   
        double *km;           # mean wave number of all grids inside the bin 
        size_t *cnt;          # number of grids inside the bin               
        double *lcnt;         # number of weighted grids for multipoles      
        double **pl[2];       # auto power spectrum multipoles               
        double **xpl;         # cross power spectrum multipoles              
        int nomp;             # number of OpenMP threads                     
        double *pcnt;         # for summing the Fourier space densities      
        double *plcnt;        # for counting the number of multipole modes   
    void powspec_destroy(PK *pk) nogil

    

cdef CATA* numpy_to_cata_auto(double[:,:] positions_a,) nogil:
    cdef int ndata = 1;
    cdef CATA* cat = cata_init(ndata)
    
    cdef size_t i, j, k
    cat.num = 1
    cat.data[0] =  <DATA*> malloc(sizeof(DATA) * positions_a.shape[0])
    cat.ndata[0] = positions_a.shape[0]
    

    

    for i in range(ndata):
        for j in range(positions_a.shape[0]):
            cat.data[i][j].x[0] = positions_a[j,0]
            cat.data[i][j].x[1] = positions_a[j,1]
            cat.data[i][j].x[2] = positions_a[j,2]
            cat.data[i][j].w = positions_a[j,3]
            cat.wdata[i] += cat.data[i][j].w
    
    return cat

    
def compute_auto_powspec(double[:] data_x, #Assumes double precision input/FFTW!
                        double[:] data_y, 
                        double[:] data_z, 
                        double[:] data_w,
                        powspec_conf_file = "powspec.conf",
                        output_file = None):


    save_out = output_file is not None
    if not save_out:
        # Define dummy names for IO so conf does not crash
        test_output = "--auto=test/test.out"
    else:
        test_output = f"--auto={output_file}"
    test_output_bytes = test_output.encode('utf-8') + b'\x00'
    cdef char* test_output_string = test_output_bytes

    # Define name of the configuration file to use
    # TODO: Generate temporary configuration file at fixed location
    #       from options passed to function. See i.e. 
    #       https://github.com/dforero0896/fcfcwrap
    # TODO: (Alternative/harder) override CONF structure
    conf = f"--conf={powspec_conf_file}"
    conf_bytes = conf.encode('utf-8') + b'\x00'
    cdef char* conf_string = conf_bytes

    # Define dummy argc, argv to send to powspec main function
    # This should remain similar once we generate a conf file.
    cdef int argc = 2
    cdef char* argv[2]
    argv[0] = conf_string
    argv[1] = test_output_string

    # Create CATA structure (involves data copying)
    cat = numpy_to_cata_auto(np.c_[data_x, data_y, data_z, data_w])

    cdef PK* pk = compute_pk(cat, <bint> save_out, argc, argv)
    pk_result = {}
    pk_result['n_data_objects'] = cat.ndata[0]
    pk_result['w_data_objects'] = cat.wdata[0]
    pk_result['shot_noise'] = cat.shot[0]
    pk_result['normalisation'] = cat.norm[0]
    pk_result['k'] = np.empty(pk.nbin, dtype = np.double) 
    pk_result['kmin'] = np.empty(pk.nbin, dtype = np.double) 
    pk_result['kmax'] = np.empty(pk.nbin, dtype = np.double) 
    pk_result['kavg'] = np.empty(pk.nbin, dtype = np.double) 
    pk_result['nmodes'] = np.empty(pk.nbin, dtype = np.double) 
    pk_result['multipoles'] = np.empty((pk.nbin, pk.nl), dtype = np.double) 
    


    # Loop assumes computing l=0,2,4
    
    for i in range(pk.nbin):
        pk_result['k'][i] = pk.k[i]
        pk_result['kavg'][i] = pk.km[i]
        pk_result['kmin'][i] = pk.kedge[i]
        pk_result['kmax'][i] = pk.kedge[i+1]
        pk_result['nmodes'][i] = pk.cnt[i]
        for j in range(pk.nl):
            pk_result['multipoles'][i,j] = pk.pl[0][j][i]
        
        #printf("%lf %lf %lf %lf\n", pk.k[i], pk.pl[0][0][i], pk.pl[0][1][i], pk.pl[0][2][i])

    cata_destroy(cat)
    powspec_destroy(pk)
    
    return pk_result


