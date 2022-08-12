#cython: language_level=3
#cython: boundscheck = False
import cython
from cython.parallel import prange, threadid
cimport openmp
from libc.stdlib cimport malloc, free, calloc
from libc.stdio cimport FILE, fprintf, fopen, fclose, printf, fflush, stdout, stderr
import numpy as np
from libc.math cimport sqrt


arg0_bytes = "POWSPEC".encode('utf-8') + b'\x00'
cdef char* arg0_str = arg0_bytes


# load_conf.h

cdef extern from "load_conf.h":
    ctypedef struct CONF:
        char *fconf;          # Name of the configuration file. */
        int ndata;            # Number of data catalogs. */
        char **dfname;        # DATA_CATALOG    */
        char **rfname;        # RAND_CATALOG    */
        int *dftype;          # DATA_FORMAT     */
        int *rftype;          # RAND_FORMAT     */
        long *dskip;          # DATA_SKIP       */
        long *rskip;          # RAND_SKIP       */
        char *dcmt;           # DATA_COMMENT    */
        char *rcmt;           # RAND_COMMENT    */
        char **dfmtr;         # DATA_FORMATTER  */
        char **rfmtr;         # RAND_FORMATTER  */
        char **dpos;          # DATA_POSITION   */
        char **rpos;          # RAND_POSITION   */
        char **dwcomp;        # DATA_WT_COMP    */
        char **rwcomp;        # RAND_WT_COMP    */
        char **dwfkp;         # DATA_WT_FKP     */
        char **rwfkp;         # RAND_WT_FKP     */
        char **dnz;           # DATA_NZ         */
        char **rnz;           # RAND_NZ         */
        char **dsel;          # DATA_SELECTION  */
        char **rsel;          # RAND_SELECTION  */
        bint has_asc[2];      # Indicate whether there is an ASCII file. */
        bint *dcnvt;          # DATA_CONVERT    */
        bint *rcnvt;          # RAND_CONVERT    */
        bint cnvt;            # Indicate if coordinate conversion is required. */
        double omega_m;       # OMEGA_M         */
        double omega_l;       # OMEGA_LAMBDA    */
        double omega_k;       # 1 - OMEGA_M - OMEGA_LAMBDA */
        double eos_w;         # DE_EOS_W        */
        double ecdst;         # CMVDST_ERR      */
        char *fcdst;          # Z_CMVDST_CONV   */
        bint issim;           # CUBIC_SIM       */
        double *los;          # LINE_OF_SIGHT   */
        double *bsize;        # BOX_SIZE        */
        double *bpad;         # BOX_PAD         */
        int gsize;            # GRID_SIZE       */
        int assign;           # PARTICLE_ASSIGN */
        bint intlace;         # GRID_INTERLACE  */
        int *poles;           # MULTIPOLE       */
        int npole;            # Number of multipoles to be computed. */
        double kmin;          # KMIN            */
        double kmax;          # KMAX            */
        bint logscale;        # LOG_SCALE       */
        double kbin;          # BIN_SIZE        */
        char **oauto;         # OUTPUT_AUTO     */
        char *ocross;         # OUTPUT_CROSS    */
        bint isauto[2];       # Indicate whether to compute auto power spectra. */
        bint iscross;         # Indicate whether to compute cross power spectra. */
        bint oheader;         # OUTPUT_HEADER   */
        int ovwrite;          # OVERWRITE       */
        bint verbose;         # VERBOSE         */
    void conf_destroy(CONF *conf) nogil
    CONF *conf_init() nogil
    int check_cosmo_lib(double *omega_m, double *omega_l,
                    double *omega_k, double *eos_w) nogil
        
cdef CONF* init_conf(int ndata, 
                     bint issim, 
                     bint* dcnvt, 
                     bint* rcnvt,
                     double omega_m,
                     double omega_l,
                     double omega_k,
                     double eos_w,
                     double* los,
                     double bsize,
                     double bpad
                     ) nogil:

    cdef CONF* conf = conf_init()
    conf.ndata = ndata
    conf.issim = issim
    conf.cnvt = False
    conf.dcnvt = dcnvt
    conf.rcnvt = rcnvt
    conf.omega_m = omega_m
    conf.omega_l = omega_l
    conf.omega_k = omega_k
    conf.eos_w = eos_w
    conf.los = los

    if not conf.issim:
        if conf.dcnvt:
            for i in range(conf.ndata):
                if conf.dcnvt[i]:
                    conf.cnvt = True
                    break
        if not conf.cnvt:
            if conf.rcnvt:
                for i in range(conf.ndata):
                    if conf.rcnvt[i]:
                        conf.cnvt = True
                        break
    if conf.cnvt:
        check_cosmo_lib(&conf.omega_m, &conf.omega_l, &conf.omega_k,
                        &conf.eos_w)
        conf.ecdst = 1e-8
    cdef double tmp = 0
    if conf.issim:
        if conf.los:
            tmp = sqrt(conf.los[0] * conf.los[0] + conf.los[1] * conf.los[1] + conf.los[2] * conf.los[2])
            conf.los[0] /= tmp
            conf.los[1] /= tmp
            conf.los[2] /= tmp
        else:
            conf.los = <double*> calloc(3, sizeof(double))
            conf.los[2] = 1





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

    PK *compute_pk(CATA* cata, bint save_out, bint has_randoms, int argc, char* argv[]) nogil;

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

def pystring_to_cstring(str string):

    
    str_bytes = string.encode('utf-8') + b'\x00'
    cdef char* cstring = str_bytes
    return cstring

cdef double* numpy_to_cata(double[:,:] positions,
                        CATA* cat,
                        int data_id,
                        bint is_rand,
                        bint is_sim) nogil:  
    cdef size_t i, j, k
    cdef double *sumw2 = <double*> malloc(sizeof(double) * 2)
    sumw2[0] = sumw2[1] = 0
    # sumw2n is sumw2[1]
    if not is_rand:
        cat.data[data_id] = <DATA*> malloc(sizeof(DATA) * positions.shape[0])
        cat.ndata[data_id] = positions.shape[0]
        for j in range(positions.shape[0]):
            cat.data[data_id][j].x[0] = positions[j,0]
            cat.data[data_id][j].x[1] = positions[j,1]
            cat.data[data_id][j].x[2] = positions[j,2]
            cat.wdata[data_id] += positions[j,3]
            if not is_sim:
                cat.data[data_id][j].w = positions[j,3] * positions[j,4]
                sumw2[0] += cat.data[data_id][j].w**2
                sumw2[1] += positions[j,3] * positions[j,4]**2 * positions[j,5]
            else:
                cat.data[data_id][j].w = positions[j,3]
    else:
        cat.rand[data_id] = <DATA*> malloc(sizeof(DATA) * positions.shape[0])
        cat.nrand[data_id] = positions.shape[0]
        for j in range(positions.shape[0]):
            cat.rand[data_id][j].x[0] = positions[j,0]
            cat.rand[data_id][j].x[1] = positions[j,1]
            cat.rand[data_id][j].x[2] = positions[j,2]
            cat.wrand[data_id] += positions[j,3]
            if not is_sim:
                cat.rand[data_id][j].w = positions[j,3] * positions[j,4]
                sumw2[0] += cat.rand[data_id][j].w**2
                sumw2[1] += positions[j,3] * positions[j,4]**2 * positions[j,5]
            else:
                cat.rand[data_id][j].w = positions[j,3]
    return sumw2



    
def compute_auto_box(double[:] data_x, #Assumes double precision input/FFTW!
                        double[:] data_y, 
                        double[:] data_z, 
                        double[:] data_w,
                        powspec_conf_file,
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
    cdef int argc = 3
    cdef char* argv[3]
    argv[0] = arg0_str
    argv[1] = conf_string
    argv[2] = test_output_string

    # Create CATA structure (involves data copying)
    cdef int ndata = 1;
    cdef CATA* cat = cata_init(ndata)
    
    cdef double* sumw2 = numpy_to_cata(np.c_[data_x, data_y, data_z, data_w],
                    cat,
                    0,
                    0,
                    True)

    cdef PK* pk = compute_pk(cat, <bint> save_out, False, argc, argv)
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
    free(sumw2)
    
    return pk_result


def compute_cross_box(double[:] data_1_x, #Assumes double precision input/FFTW!
                      double[:] data_1_y, 
                      double[:] data_1_z, 
                      double[:] data_1_w,
                      double[:] data_2_x, #Assumes double precision input/FFTW!
                      double[:] data_2_y, 
                      double[:] data_2_z, 
                      double[:] data_2_w,
                      powspec_conf_file,
                      output_auto = None,
                      output_cross = None):

    
    save_auto = output_auto is not None
    if not save_auto:
        # Define dummy names for IO so conf does not crash
        auto_output = "--auto=[test/auto_1.out,test/auto_2.out]"
        cross_output = "--cross=test/cross.out"
    else:
        auto_output = f"--auto=[{','.join(output_auto)}]"
        cross_output = f"--cross={output_cross}"
    auto_output_bytes = auto_output.encode('utf-8') + b'\x00'
    cdef char* auto_output_string = auto_output_bytes

    cross_output_bytes = cross_output.encode('utf-8') + b'\x00'
    cdef char* cross_output_string = cross_output_bytes

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
    cdef int argc = 4
    cdef char* argv[4]
    argv[0] = arg0_str
    argv[1] = conf_string
    argv[2] = auto_output_string
    argv[3] = cross_output_string

    # Create CATA structure (involves data copying)
    cdef int ndata = 2;
    cdef CATA* cat = cata_init(ndata)
    
    cdef double* sumw2_a = numpy_to_cata(np.c_[data_1_x, data_1_y, data_1_z, data_1_w],
                    cat,
                    0,
                    False,
                    True)
    cdef double* sumw2_b = numpy_to_cata(np.c_[data_2_x, data_2_y, data_2_z, data_2_w],
                    cat,
                    1,
                    False,
                    True)


    cdef PK* pk = compute_pk(cat, <bint> save_auto, False, argc, argv)
    pk_result = {}
    pk_result['n_data_objects'] = [cat.ndata[i] for i in range(cat.num)]
    pk_result['w_data_objects'] = [cat.wdata[i] for i in range(cat.num)]
    pk_result['shot_noise'] = [cat.shot[i] for i in range(cat.num)]
    pk_result['normalisation'] = [cat.norm[i] for i in range(cat.num)]
    pk_result['k'] = np.empty(pk.nbin, dtype = np.double) 
    pk_result['kmin'] = np.empty(pk.nbin, dtype = np.double) 
    pk_result['kmax'] = np.empty(pk.nbin, dtype = np.double) 
    pk_result['kavg'] = np.empty(pk.nbin, dtype = np.double) 
    pk_result['nmodes'] = np.empty(pk.nbin, dtype = np.double) 
    pk_result['auto_multipoles'] = np.empty((ndata, pk.nbin, pk.nl), dtype = np.double) 
    pk_result['cross_multipoles'] = np.empty((pk.nbin, pk.nl), dtype = np.double) 
    

   
    for i in range(pk.nbin):
        pk_result['k'][i] = pk.k[i]
        pk_result['kavg'][i] = pk.km[i]
        pk_result['kmin'][i] = pk.kedge[i]
        pk_result['kmax'][i] = pk.kedge[i+1]
        pk_result['nmodes'][i] = pk.cnt[i]
        for j in range(pk.nl):
            pk_result['cross_multipoles'][i,j] = pk.xpl[j][i]
            for k in range(ndata):
                pk_result['auto_multipoles'][k,i,j] = pk.pl[k][j][i]
        
        #printf("%lf %lf %lf %lf\n", pk.k[i], pk.pl[0][0][i], pk.pl[0][1][i], pk.pl[0][2][i])

    cata_destroy(cat)
    powspec_destroy(pk)
    free(sumw2_a)
    free(sumw2_b)
    
    return pk_result


def compute_auto_lc(double[:] data_x, #Assumes double precision input/FFTW!
                    double[:] data_y, 
                    double[:] data_z, 
                    double[:] data_wcomp,
                    double[:] data_wfkp,
                    double[:] data_nz,
                    double[:] rand_x, #Assumes double precision input/FFTW!
                    double[:] rand_y, 
                    double[:] rand_z, 
                    double[:] rand_wcomp,
                    double[:] rand_wfkp,
                    double[:] rand_nz,
                    powspec_conf_file,
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
    cdef int argc = 3
    cdef char* argv[3]
    argv[0] = arg0_str
    argv[1] = conf_string
    argv[2] = test_output_string

    # Create CATA structure (involves data copying)
    cdef int ndata = 1;
    cdef int i = 0
    cdef CATA* cat = cata_init(ndata)
    
    cdef double* sumw2_dat = numpy_to_cata(np.c_[data_x, data_y, data_z, data_wcomp, data_wfkp, data_nz],
                    cat,
                    i,
                    0,
                    False)

    cdef double* sumw2_ran = numpy_to_cata(np.c_[rand_x, rand_y, rand_z, rand_wcomp, rand_wfkp, rand_nz],
                    cat,
                    i,
                    1,
                    False)
    
    cat.alpha[i] = cat.wdata[i] / cat.wrand[i]
    # Shot noise from both data and random 
    cat.shot[i] = sumw2_dat[0] + cat.alpha[i] * cat.alpha[i] * sumw2_ran[0]
    # Normalise using the random. 
    if (sumw2_dat[1] == 0): cat.norm[i] = cat.alpha[i] * sumw2_ran[1]
    #Normalise using the data. 
    elif (sumw2_ran[1] == 0): cat.norm[i] = sumw2_dat[1]
    # Check consistency (TODO) and Normalise using the random. 
    else: cat.norm[i] = cat.alpha[i] * sumw2_ran[1]
    
    


    cdef PK* pk = compute_pk(cat, <bint> save_out, True, argc, argv)
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
    free(sumw2_dat)
    free(sumw2_ran)
    
    return pk_result

def compute_cross_lc(double[:] data_1_x, #Assumes double precision input/FFTW!
                     double[:] data_1_y, 
                     double[:] data_1_z, 
                     double[:] data_1_wcomp,
                     double[:] data_1_wfkp,
                     double[:] data_1_nz,
                     double[:] rand_1_x, #Assumes double precision input/FFTW!
                     double[:] rand_1_y, 
                     double[:] rand_1_z, 
                     double[:] rand_1_wcomp,
                     double[:] rand_1_wfkp,
                     double[:] rand_1_nz,
                     double[:] data_2_x, #Assumes double precision input/FFTW!
                     double[:] data_2_y, 
                     double[:] data_2_z, 
                     double[:] data_2_wcomp,
                     double[:] data_2_wfkp,
                     double[:] data_2_nz,
                     double[:] rand_2_x, #Assumes double precision input/FFTW!
                     double[:] rand_2_y, 
                     double[:] rand_2_z, 
                     double[:] rand_2_wcomp,
                     double[:] rand_2_wfkp,
                     double[:] rand_2_nz,
                     powspec_conf_file,
                     output_auto = None,
                     output_cross = None):
 
    save_auto = output_auto is not None
    if not save_auto:
        # Define dummy names for IO so conf does not crash
        auto_output = "--auto=[test/auto_1.out,test/auto_2.out]"
        cross_output = "--cross=test/cross.out"
    else:
        auto_output = f"--auto=[{','.join(output_auto)}]"
        cross_output = f"--cross={output_cross}"
    auto_output_bytes = auto_output.encode('utf-8') + b'\x00'
    cdef char* auto_output_string = auto_output_bytes

    cross_output_bytes = cross_output.encode('utf-8') + b'\x00'
    cdef char* cross_output_string = cross_output_bytes


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
    cdef int argc = 4
    cdef char* argv[4]
    argv[0] = arg0_str
    argv[1] = conf_string
    argv[2] = auto_output_string
    argv[3] = cross_output_string

    # Create CATA structure (involves data copying)
    cdef int ndata = 2;
    cdef int i = 0
    cdef CATA* cat = cata_init(ndata)
    
    cdef double* sumw2_dat = numpy_to_cata(np.c_[data_1_x, data_1_y, data_1_z, data_1_wcomp, data_1_wfkp, data_1_nz],
                    cat,
                    i,
                    0,
                    False)

    cdef double* sumw2_ran = numpy_to_cata(np.c_[rand_1_x, rand_1_y, rand_1_z, rand_1_wcomp, rand_1_wfkp, rand_1_nz],
                    cat,
                    i,
                    1,
                    False)
    
    cat.alpha[i] = cat.wdata[i] / cat.wrand[i]
    # Shot noise from both data and random 
    cat.shot[i] = sumw2_dat[0] + cat.alpha[i] * cat.alpha[i] * sumw2_ran[0]
    # Normalise using the random. 
    if (sumw2_dat[1] == 0): cat.norm[i] = cat.alpha[i] * sumw2_ran[1]
    #Normalise using the data. 
    elif (sumw2_ran[1] == 0): cat.norm[i] = sumw2_dat[1]
    # Check consistency (TODO) and Normalise using the random. 
    else: cat.norm[i] = cat.alpha[i] * sumw2_ran[1]
    print(sumw2_dat[0], sumw2_dat[1])
    print(sumw2_ran[0], sumw2_ran[1])
    print(cat.shot[i])

    # Repeat with second set of catalogs

    i = 1
    sumw2_dat[0] = 0
    sumw2_dat[1] = 0
    sumw2_ran[0] = 0
    sumw2_ran[1] = 0

    sumw2_dat = numpy_to_cata(np.c_[data_2_x, data_2_y, data_2_z, data_2_wcomp, data_2_wfkp, data_2_nz],
                    cat,
                    i,
                    0,
                    False)

    sumw2_ran = numpy_to_cata(np.c_[rand_2_x, rand_2_y, rand_2_z, rand_2_wcomp, rand_2_wfkp, rand_2_nz],
                    cat,
                    i,
                    1,
                    False)
    
    cat.alpha[i] = cat.wdata[i] / cat.wrand[i]
    # Shot noise from both data and random 
    cat.shot[i] = sumw2_dat[0] + cat.alpha[i] * cat.alpha[i] * sumw2_ran[0]
    # Normalise using the random. 
    if (sumw2_dat[1] == 0): cat.norm[i] = cat.alpha[i] * sumw2_ran[1]
    #Normalise using the data. 
    elif (sumw2_ran[1] == 0): cat.norm[i] = sumw2_dat[1]
    # Check consistency (TODO) and Normalise using the random. 
    else: cat.norm[i] = cat.alpha[i] * sumw2_ran[1]
    print(sumw2_dat[0], sumw2_dat[1])
    print(sumw2_ran[0], sumw2_ran[1])
    print(cat.shot[i])
    


    cdef PK* pk = compute_pk(cat, <bint> save_auto, True, argc, argv)
    pk_result = {}
    pk_result['n_data_objects'] = [cat.ndata[i] for i in range(cat.num)]
    pk_result['w_data_objects'] = [cat.wdata[i] for i in range(cat.num)]
    pk_result['shot_noise'] = [cat.shot[i] for i in range(cat.num)]
    pk_result['normalisation'] = [cat.norm[i] for i in range(cat.num)]
    pk_result['k'] = np.empty(pk.nbin, dtype = np.double) 
    pk_result['kmin'] = np.empty(pk.nbin, dtype = np.double) 
    pk_result['kmax'] = np.empty(pk.nbin, dtype = np.double) 
    pk_result['kavg'] = np.empty(pk.nbin, dtype = np.double) 
    pk_result['nmodes'] = np.empty(pk.nbin, dtype = np.double) 
    pk_result['auto_multipoles'] = np.empty((ndata, pk.nbin, pk.nl), dtype = np.double) 
    pk_result['cross_multipoles'] = np.empty((pk.nbin, pk.nl), dtype = np.double) 
    

   
    for i in range(pk.nbin):
        pk_result['k'][i] = pk.k[i]
        pk_result['kavg'][i] = pk.km[i]
        pk_result['kmin'][i] = pk.kedge[i]
        pk_result['kmax'][i] = pk.kedge[i+1]
        pk_result['nmodes'][i] = pk.cnt[i]
        for j in range(pk.nl):
            pk_result['cross_multipoles'][i,j] = pk.xpl[j][i]
            for k in range(ndata):
                pk_result['auto_multipoles'][k,i,j] = pk.pl[k][j][i]
        
    cata_destroy(cat)
    powspec_destroy(pk)
    free(sumw2_dat)
    free(sumw2_ran)
    
    return pk_result

def compute_auto_box_rand(double[:] data_x, #Assumes double precision input/FFTW!
                          double[:] data_y, 
                          double[:] data_z, 
                          double[:] data_w,
                          double[:] rand_x, #Assumes double precision input/FFTW!
                          double[:] rand_y, 
                          double[:] rand_z, 
                          double[:] rand_w,
                          powspec_conf_file,
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
    cdef int argc = 3
    cdef char* argv[3]
    argv[0] = arg0_str
    argv[1] = conf_string
    argv[2] = test_output_string

    # Create CATA structure (involves data copying)
    cdef int ndata = 1;
    cdef int i = 0;
    cdef CATA* cat = cata_init(ndata)
    
    
    cdef double* sumw2_dat = numpy_to_cata(np.c_[data_x, data_y, data_z, data_w],
                    cat,
                    i,
                    0,
                    True)

    cdef double* sumw2_ran = numpy_to_cata(np.c_[rand_x, rand_y, rand_z, rand_w],
                    cat,
                    i,
                    1,
                    True)
    
    cat.alpha[i] = cat.wdata[i] / cat.wrand[i]
    # Shot noise from both data and random 
    cat.shot[i] = sumw2_dat[0] + cat.alpha[i] * cat.alpha[i] * sumw2_ran[0]
    # Normalise using the random. 
    if (sumw2_dat[1] == 0): cat.norm[i] = cat.alpha[i] * sumw2_ran[1]
    #Normalise using the data. 
    elif (sumw2_ran[1] == 0): cat.norm[i] = sumw2_dat[1]
    # Check consistency (TODO) and Normalise using the random. 
    else: cat.norm[i] = cat.alpha[i] * sumw2_ran[1]

    cdef PK* pk = compute_pk(cat, <bint> save_out, True, argc, argv)
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
    free(sumw2_dat)
    free(sumw2_ran)
    
    return pk_result

def compute_cross_box_rand(double[:] data_1_x, #Assumes double precision input/FFTW!
                           double[:] data_1_y, 
                           double[:] data_1_z, 
                           double[:] data_1_w,
                           double[:] rand_1_x, #Assumes double precision input/FFTW!
                           double[:] rand_1_y, 
                           double[:] rand_1_z, 
                           double[:] rand_1_w,
                           double[:] data_2_x, #Assumes double precision input/FFTW!
                           double[:] data_2_y, 
                           double[:] data_2_z, 
                           double[:] data_2_w,
                           double[:] rand_2_x, #Assumes double precision input/FFTW!
                           double[:] rand_2_y, 
                           double[:] rand_2_z, 
                           double[:] rand_2_w,
                           powspec_conf_file,
                           output_auto = None,
                           output_cross = None):

    
    save_auto = output_auto is not None
    if not save_auto:
        # Define dummy names for IO so conf does not crash
        auto_output = "--auto=[test/auto_1.out,test/auto_2.out]"
        cross_output = "--cross=test/cross.out"
    else:
        auto_output = f"--auto=[{','.join(output_auto)}]"
        cross_output = f"--cross={output_cross}"
    auto_output_bytes = auto_output.encode('utf-8') + b'\x00'
    cdef char* auto_output_string = auto_output_bytes

    cross_output_bytes = cross_output.encode('utf-8') + b'\x00'
    cdef char* cross_output_string = cross_output_bytes

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
    cdef int argc = 4
    cdef char* argv[4]
    argv[0] = arg0_str
    argv[1] = conf_string
    argv[2] = auto_output_string
    argv[3] = cross_output_string

    # Create CATA structure (involves data copying)
    cdef int ndata = 2;
    cdef int i = 0;
    cdef CATA* cat = cata_init(ndata)
    
    
    cdef double* sumw2_dat = numpy_to_cata(np.c_[data_1_x, data_1_y, data_1_z, data_1_w],
                    cat,
                    i,
                    0,
                    True)

    cdef double* sumw2_ran = numpy_to_cata(np.c_[rand_1_x, rand_1_y, rand_1_z, rand_1_w],
                    cat,
                    i,
                    1,
                    True)
    
    cat.alpha[i] = cat.wdata[i] / cat.wrand[i]
    # Shot noise from both data and random 
    cat.shot[i] = sumw2_dat[0] + cat.alpha[i] * cat.alpha[i] * sumw2_ran[0]
    # Normalise using the random. 
    if (sumw2_dat[1] == 0): cat.norm[i] = cat.alpha[i] * sumw2_ran[1]
    #Normalise using the data. 
    elif (sumw2_ran[1] == 0): cat.norm[i] = sumw2_dat[1]
    # Check consistency (TODO) and Normalise using the random. 
    else: cat.norm[i] = cat.alpha[i] * sumw2_ran[1]

    i = 1;
    sumw2_dat[0] = 0
    sumw2_dat[1] = 0
    sumw2_ran[0] = 0
    sumw2_ran[1] = 0

    sumw2_dat = numpy_to_cata(np.c_[data_2_x, data_2_y, data_2_z, data_2_w],
                    cat,
                    i,
                    0,
                    True)

    sumw2_ran = numpy_to_cata(np.c_[rand_2_x, rand_2_y, rand_2_z, rand_2_w],
                    cat,
                    i,
                    1,
                    True)


    cat.alpha[i] = cat.wdata[i] / cat.wrand[i]
    # Shot noise from both data and random 
    cat.shot[i] = sumw2_dat[0] + cat.alpha[i] * cat.alpha[i] * sumw2_ran[0]
    # Normalise using the random. 
    if (sumw2_dat[1] == 0): cat.norm[i] = cat.alpha[i] * sumw2_ran[1]
    #Normalise using the data. 
    elif (sumw2_ran[1] == 0): cat.norm[i] = sumw2_dat[1]
    # Check consistency (TODO) and Normalise using the random. 
    else: cat.norm[i] = cat.alpha[i] * sumw2_ran[1]

    cdef PK* pk = compute_pk(cat, <bint> save_auto, True, argc, argv)
    pk_result = {}
    pk_result['n_data_objects'] = [cat.ndata[i] for i in range(cat.num)]
    pk_result['w_data_objects'] = [cat.wdata[i] for i in range(cat.num)]
    pk_result['shot_noise'] = [cat.shot[i] for i in range(cat.num)]
    pk_result['normalisation'] = [cat.norm[i] for i in range(cat.num)]
    pk_result['k'] = np.empty(pk.nbin, dtype = np.double) 
    pk_result['kmin'] = np.empty(pk.nbin, dtype = np.double) 
    pk_result['kmax'] = np.empty(pk.nbin, dtype = np.double) 
    pk_result['kavg'] = np.empty(pk.nbin, dtype = np.double) 
    pk_result['nmodes'] = np.empty(pk.nbin, dtype = np.double) 
    pk_result['auto_multipoles'] = np.empty((ndata, pk.nbin, pk.nl), dtype = np.double) 
    pk_result['cross_multipoles'] = np.empty((pk.nbin, pk.nl), dtype = np.double) 
    

   
    for i in range(pk.nbin):
        pk_result['k'][i] = pk.k[i]
        pk_result['kavg'][i] = pk.km[i]
        pk_result['kmin'][i] = pk.kedge[i]
        pk_result['kmax'][i] = pk.kedge[i+1]
        pk_result['nmodes'][i] = pk.cnt[i]
        for j in range(pk.nl):
            pk_result['cross_multipoles'][i,j] = pk.xpl[j][i]
            for k in range(ndata):
                pk_result['auto_multipoles'][k,i,j] = pk.pl[k][j][i]
        
        #printf("%lf %lf %lf %lf\n", pk.k[i], pk.pl[0][0][i], pk.pl[0][1][i], pk.pl[0][2][i])

    cata_destroy(cat)
    powspec_destroy(pk)
    free(sumw2_dat)
    free(sumw2_ran)
    
    return pk_result