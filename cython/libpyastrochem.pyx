from libc.stdlib cimport malloc, free
from libc.string cimport strcmp
from cpython.string cimport PyString_AsString

cdef extern from "../src/libastrochem.h":
    cdef double CHI_DEFAULT
    cdef double COSMIC_DEFAULT
    cdef double GRAIN_SIZE_DEFAULT
    cdef double ABS_ERR_DEFAULT
    cdef double REL_ERR_DEFAULT

    ctypedef struct params_t:
        double nh
        double tgas
        double tdust

    ctypedef struct phys_t:
        double chi
        double cosmic
        double grain_size
        double grain_abundance

    ctypedef struct net_t:
        int n_species
        char **species_names

    ctypedef struct cell_t:
        double av
        double nh
        double tgas
        double tdust

    ctypedef struct astrochem_mem_t:
        params_t params

    void read_network ( const char *chem_file, net_t* network, const int verbose)
    void free_network ( net_t* network)
    int alloc_abundances( const net_t* network, double** abundances )
    void free_abundances( double* abundances )
    int set_initial_abundances( char** species, int n_initialized_abundances, const double* initial_abundances, const net_t* network, double* abundances )
    int solver_init( const cell_t* cell, const net_t* network, const phys_t* phys, const double* abundances , double density, double abs_err, double rel_err, astrochem_mem_t* astrochem_mem )
    int solve( const astrochem_mem_t* astrochem_mem, const net_t* network, double* abundances, double time , const cell_t* new_cell, int verbose )
    void solver_close( astrochem_mem_t* astrochem_mem )

_ABS_ERR_DEFAULT = ABS_ERR_DEFAULT
_REL_ERR_DEFAULT = REL_ERR_DEFAULT

cdef class Network:
    """
    Chemical Network class

    Network( const char* chem_file, int verbose )
    Create a chemical network

    Arguments:
    chem_file       -- File containing a network to load
    verbose         -- Verbose if 1, Quiet if 0

    """
    cdef net_t thisstruct
    def __cinit__( self, const char* chem_file, int verbose ):
        read_network( chem_file, &self.thisstruct, verbose )

    def __dealloc__(self):
        free_network( &self.thisstruct )

cdef class Phys:
    """
    Physical parameters to use in chemical reaction solver
    """
    cdef public phys_t thisstruct

    def __cinit__( self ):
        self.thisstruct.chi = CHI_DEFAULT
        self.thisstruct.cosmic = COSMIC_DEFAULT
        self.thisstruct.grain_size = GRAIN_SIZE_DEFAULT
        self.thisstruct.grain_abundance = 0

    property chi:
        """
        Chi physical property
        """
        def __get__(self):
            return self.thisstruct.chi
        def __set__(self, double chi):
            self.thisstruct.chi = chi

    property cosmic:
        """
        Cosmic physical property
        """
        def __get__(self):
            return self.thisstruct.cosmic
        def __set__(self, double cosmic):
            self.thisstruct.cosmic = cosmic

    property grain_size:
        """
        Grain Size physical property
        """
        def __get__(self):
            return self.thisstruct.grain_size
        def __set__(self, double grain_size):
            self.thisstruct.grain_size = grain_size

    property grain_abundance:
        """
        Grain Abundance physical property
        """
        def __get__(self):
            return self.thisstruct.grain_abundance
        def __set__(self, double grain_abundance):
            self.thisstruct.grain_abundance = grain_abundance

cdef class Cell:
    """
    Chemical cell class

    Cell( double av, double nh, double tgas, double tdust )
    Create a chemical cell class

    Arguments:
    av    -- av to use
    nh    -- nh to use
    tgas  -- tgas to use
    tdust -- tdust to use
    """
    cdef public cell_t thisstruct

    def __cinit__( self, double av, double nh, double tgas, double tdust ):
        self.thisstruct.av = av
        self.thisstruct.nh = nh
        self.thisstruct.tdust = tgas
        self.thisstruct.tgas = tdust

    property av:
        """
        av cell property
        """
        def __get__(self):
            return self.thisstruct.av
        def __set__(self, double av):
            self.thisstruct.av = av

    property nh:
        """
        nh cell property
        """
        def __get__(self):
            return self.thisstruct.nh
        def __set__(self, double nh):
            self.thisstruct.nh = nh

    property tgas:
        """
        tgas cell property
        """
        def __get__(self):
            return self.thisstruct.tgas
        def __set__(self, double tgas):
            self.thisstruct.tgas = tgas

    property tdust:
        """
        tdust cell property
        """
        def __get__(self):
            return self.thisstruct.tdust
        def __set__(self, double tdust):
            self.thisstruct.tdust = tdust

cdef class Solver:
    """
    Chemical reaction solver

    Solver( cell , const char* chem_file, phys, abs_err, rel_err, initial_abundances , double density, int verbose )
    Create a chemical reaction solver

    Arguments:
    cell                 -- Chemical cell class to use in solver
    chem_file            -- Chemical network file string to load network from and use in solver
    phys                 -- Physical properties class to use in solver
    abs_err              -- Absolute acceptable error to use in solver
    rel_err              -- Relative acceptable error to use in solver
    initial_abundances   -- Initial abundances dictionnary with format {Species:Value}
    density              -- Density to use in solver
    verbose              -- verbose if 1, quiet if 0
    """
    cdef astrochem_mem_t astrochemstruct
    cdef double* abundances
    cdef Network network
    cdef int verbose

    def __cinit__( self , cell , const char* chem_file, phys, abs_err, rel_err, initial_abundances , double density, int verbose ):
        self.verbose = verbose
        self.network = Network( chem_file, verbose )
        cdef net_t c_net = self.network.thisstruct
        cdef phys_t c_phys = phys.thisstruct
        cdef cell_t c_cell = cell.thisstruct

        if alloc_abundances( &c_net , &self.abundances ) != 0 :
            raise MemoryError
        cdef char **initial_abundances_str = <char **>malloc(len(initial_abundances) * sizeof(char *))
        cdef double* initial_abundances_val = <double*>malloc(len(initial_abundances) * sizeof(double))
        cdef int j = 0
        for i in initial_abundances:
            initial_abundances_str[j] = PyString_AsString(i)
            initial_abundances_val[j] = initial_abundances[i]
            j+=1
        set_initial_abundances( < const char** >initial_abundances_str, len(initial_abundances), initial_abundances_val,  &c_net, self.abundances )
        free( initial_abundances_str )
        free( initial_abundances_val )

        solver_init( &c_cell, &c_net, &c_phys , self.abundances, density, abs_err, rel_err, &self.astrochemstruct )

    def __dealloc__(self):
        free_abundances( self.abundances )
        solver_close(  &self.astrochemstruct )

    def solve(self, time , new_cell=0 ):
        """
        solve( time , new_cell )
        Solve chemical reaction for a certain time

        Arguments:
        time     -- time to solve the system at
        new_cell -- cell class to use in solver, optionnal
        """
        cdef net_t c_net = self.network.thisstruct
        cdef cell_t c_new_cell
        if new_cell== 0 :
            if solve( &self.astrochemstruct, &c_net , self.abundances, time , NULL, self.verbose ) != 0:
                raise ArithmeticError("Solve Error")
        else:
            c_new_cell = new_cell.thisstruct
            if solve( &self.astrochemstruct, &c_net , self.abundances, time , &c_new_cell, self.verbose ) != 0:
                raise ArithmeticError("Solve Error")
        cdef int i
        cdef bytes py_string
        ret = {}
        for i in range( c_net.n_species ):
            py_string =  c_net.species_names[i]
            ret[py_string] = self.abundances[i]
        return ret
