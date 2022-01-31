#!/usr/bin/python
import TTCF
import numpy as np


class MDInterface:
    """ This is a template class to interface with whatever MD backend code we decide to use.

        It should only include the broadest outline of functionality, which can be
        overridden to support, e.g. LAMMPS or Debra's MD code."""

    def __init__(self):
        self._tr     = 1.0
        self._drw    = 1.0
        self._drf    = 1.0
        self._delta  = 5.0E-003
        self._latt   = 1
        self._npart  = 500

        self._fe0    = 1.0
        self._rcut   = 2.5
        self._kh     = 1
        self._nprint = 100

        self._mix    = 0.2
        self._eps1   = 1
        self._eps2   = 1
        self._qvol   = 100
        
        self._eps_LJG = 1.5
        self._r_LJG = 1.47

        self._kf     = 1.0
        self._r0     = 1.0
        self._limol  = 1
        self._yzdivx = 1.0

        self._dxxdiv = 1.0

        self._ntype  = 1
        self._non    = 0
        self._ngaus  = 4
        self._e0     = 1.0

        self._nplot  = 1
        self._maxtau = 1000
        self._eqtim  = 1000
        self._ncyc   = 1

        # Lastly, set the internal flag determining whether to turn on the NEMD field
        self._do_nemd = False
        self._iflag = 0
        self._box_bounds = (None, None, None)

        self._x = None
        self._y = None
        self._z = None

        self._rdf_arrays = None

        # Finally, pass these values into the simulation
        self.setup()

    # Getters and setters (properties)
    #
    # Box coords/initial conditions. Need to reset the simulation if these change
    @property
    def yzdivx(self):
        return(self._yzdivx)
    @yzdivx.setter
    def yzdivx(self, value):
        self._yzdivx = value
        self.reset_and_update_parameters()

    @property
    def reduced_density(self):
        return(self._drf)
    @reduced_density.setter
    def reduced_density(self, value):
        self._drf=value
        self.reset_and_update_parameters()

    @property
    def npart(self):
        return(self._npart)
    @npart.setter
    def npart(self, value):
        self._npart=value
        self.reset_and_update_parameters()

    # WCA parameters. These do not require a restart when they change
    @property
    def eps(self):
        return(self._kf)
    @eps.setter
    def eps(self, value):
        self._kf=value
        self.update_parameters()
    
    # LJG parameters. These do not require a restart when they change (?)
    @property
    def epsilon_LJG(self):
    	return(self._eps_LJG)
    @epsilon_LJG.setter
    def epsilon_LJG(self, value):
    	self._eps_LJG=value
    	self.update_parameters()
    	
    @property
    def r_LJG(self):
    	return(self._r_LJG)
    @r_LJG.setter
    def r_LJG(self, value):
    	self._r_LJG=value
    	self.update_parameters()

    # Misc simulation parameters. Also don't require a restart
    @property
    def fieldstrength(self):
        return(self._fe0)
    @fieldstrength.setter
    def fieldstrength(self, value):
        self._fe0 = value
        self.update_parameters()

    @property
    def temp(self):
        return(TTCF.averg.temp)
    @temp.setter
    def temp(self, value):
        self._tr = value
        self.update_parameters()

    # Get the status of whether we're doing nonequilubrium calculations.
    # No setter, since this is handled by toggle_nemd()
    @property
    def do_nemd(self):
        return(self._do_nemd)

    # Atomic coordinates. These should be read-only to external classes
    @property
    def x(self):
        return(TTCF.coord.x[:self._npart])

    @property
    def y(self):
        return(TTCF.coord.y[:self._npart])

    @property
    def z(self):
        return(TTCF.coord.z[:self._npart])

    # Momenta. Also read-only
    @property
    def px(self):
        return(TTCF.moment.px[:self._npart])

    @property
    def py(self):
        return(TTCF.moment.py[:self._npart])

    @property
    def pz(self):
        return(TTCF.moment.pz[:self._npart])
    
    @property
    def box_bounds(self):
        return(self._box_bounds)

    @property
    def vol(self):
        return(TTCF.nopart.npart / TTCF.parm.drf)

    @property
    def timestep(self):
        return(TTCF.simul.delta)

    # Callable methods
    def setup(self):
        # This function gets called at the start of the program's run, as well as
        # whenever we restart the simulation.
        self.update_parameters()
        TTCF.setup()
        TTCF.md(1, 0)

        # Zero array for g(2) RDF
        self.delta_r = np.zeros([self.npart, self.npart])

        # Now update all of our thermodynamic parameters with values from LAMMPS
        self.get_params_from_TTCF()

    def get_params_from_TTCF(self):
        # Update this class's parameters with the values from LAMMPS. This doesn't
        # change the underlying simulation
        self._temp = TTCF.averg.temp
        self._vol = TTCF.nopart.npart / TTCF.parm.drf
        # And particle positions
        self._x = TTCF.coord.x[:self._npart]
        self._y = TTCF.coord.y[:self._npart]
        self._z = TTCF.coord.z[:self._npart]
        self._box_bounds = (TTCF.parm.cubex, TTCF.parm.cubey, TTCF.parm.cubez)

    # Radial distribution function. This returns into two arrays for r and g(2)(r), stored as a tuple
    def compute_rdf(self):
        # RDF already computed by Fortran backend during Force calculation, just need to normalise it
        rdf = self.vol / (2*np.pi*self.npart**2) * TTCF.averg.rij_hist[np.nonzero(TTCF.averg.rij_hist)]
        # Now calculate the bin coordinates
        r = np.array([i/TTCF.averg.rbin_inv for i in range(len(rdf))])
        # Finally, divide through by the width of each successive shell to normalise the probability
        for i in range(1,len(r)):
            rn = (i - 0.5)*(r[i] - r[i-1])
            rdf[i] = rdf[i]/rn
        return({'r': r, 'rdf': rdf})

    # Mean-square displacement
    def compute_msd(self):
        msdx = TTCF.count.msdx
        msdy = TTCF.count.msdy
        msdz = TTCF.count.msdz
        return((float(msdx), float(msdy), float(msdz)))

    # The backend actually calculates and stores the pressure tensor, so we need to do some processing
    # to get it into the right form for output
    def compute_pressure(self):
        p_tensor = TTCF.flux.pt
        PXY = 0.5 * (p_tensor[0,1] + p_tensor[1,0])
        PXZ = 0.5 * (p_tensor[0,1] + p_tensor[2,0])
        PYZ = 0.5 * (p_tensor[1,2] + p_tensor[2,1])
        PP = 0.5 * (p_tensor[0,0] + p_tensor[1,1] + p_tensor[2,2]) / 3
        pressure = {"PXY": PXY, "PXZ": PXZ, "PYZ": PYZ, "PP": PP} # Dict to hold all the different pressure components
        return(pressure)

    def reset_and_update_parameters(self):
        # Reset the simulation and initialise the parameters
        TTCF.inener.tr      = self._tr
        TTCF.parm.drw       = self._drw
        TTCF.parm.drf       = self._drf
        TTCF.simul.delta    = self._delta
        TTCF.iparm.latt     = self._latt
        TTCF.nopart.npart   = self._npart

        TTCF.parm.fe0       = self._fe0
        TTCF.parm.rcut      = self._rcut
        TTCF.parm.kh        = self._kh
        TTCF.iparm.nprint   = self._nprint

        TTCF.parm.mix       = self._mix
        TTCF.parm.eps1      = self._eps1
        TTCF.parm.eps2      = self._eps2
        TTCF.coord.qvol     = self._qvol
        
        TTCF.parm.eps_LJG		= self._eps_LJG
        TTCF.parm.r_LJG		= self._r_LJG

        TTCF.parm.kf        = self._kf
        TTCF.parm.r0        = self._r0
        TTCF.nopart.limol   = self._limol
        TTCF.parm.yzdivx    = self._yzdivx

        TTCF.parm.dxxdiv    = self._dxxdiv
        TTCF.iparm.ntype    = self._ntype

        TTCF.iparm.non      = self._non
        TTCF.iparm.ngaus    = self._ngaus
        TTCF.inener.e0      = self._e0

        TTCF.iparm.nplot    = self._nplot
        TTCF.iparm.maxtau   = self._maxtau
        TTCF.iparm.eqtim    = self._eqtim
        TTCF.iparm.ncyc     = self._ncyc
        TTCF.setup()

    def update_parameters(self):
        # Only update the parameters which don't require a reset
        TTCF.inener.tr      = self._tr
        TTCF.parm.drw       = self._drw
        TTCF.parm.drf       = self._drf
        TTCF.simul.delta    = self._delta
        TTCF.iparm.latt     = self._latt
        TTCF.nopart.npart   = self._npart

        TTCF.parm.fe0       = self._fe0
        TTCF.parm.rcut      = self._rcut
        TTCF.parm.kh        = self._kh
        TTCF.iparm.nprint   = self._nprint

        TTCF.parm.mix       = self._mix
        TTCF.parm.eps1      = self._eps1
        TTCF.parm.eps2      = self._eps2
        TTCF.coord.qvol     = self._qvol
        
        TTCF.parm.eps_LJG		= self._eps_LJG
        TTCF.parm.r_LJG		= self._r_LJG

        TTCF.parm.kf        = self._kf
        TTCF.parm.r0        = self._r0
        TTCF.nopart.limol   = self._limol
        TTCF.parm.yzdivx    = self._yzdivx

        TTCF.parm.dxxdiv    = self._dxxdiv
        TTCF.iparm.ntype    = self._ntype

        TTCF.iparm.non      = self._non
        TTCF.iparm.ngaus    = self._ngaus
        TTCF.inener.e0      = self._e0

        TTCF.iparm.nplot    = self._nplot
        TTCF.iparm.maxtau   = self._maxtau
        TTCF.iparm.eqtim    = self._eqtim
        TTCF.iparm.ncyc     = self._ncyc

    def toggle_nemd(self):
        if(self._do_nemd):
            self._do_nemd = False
            self._iflag = 0
            self.update_parameters()
        else:
            self._do_nemd = True
            self._iflag = 1
            self.update_parameters()

    def format_params(self):
        """ Make a format string for the input which can be either written to file or displayed in Qt."""

        # Write the input parameters, making sure to preserve whitespace and newlines
        # Write the input parameters, making sure to preserve whitespace and newlines
        param_str  =   f"""{self._tr:.3f}  {self._drw:.3f}  {self._drf:.3f}   {self._delta:.3f} {self._latt} {self._npart} 25"""
        param_str += f"""\n{self._fe0:.3f} {self._rcut:.3f} 1 {self._kh:.3f}  {self._nprint} """
        param_str += f"""\n{self._mix:.3f} {self._eps1:.3f} {self._eps2:.3f}  {self._qvol:.3f}"""
        param_str += f"""\n{self._kf:.3f}  {self._r0:.3f}   {self._limol:.3f} {self._yzdivx:.3f}"""
        param_str += f"""\n{self._dxxdiv:.3f}"""
        param_str += f"""\n{self._ntype} {self._non}    {self._ngaus} {self._e0:.3f}"""
        param_str += f"""\n{self._nplot} {self._maxtau} {self._eqtim} {self._ncyc}"""

        # Now write the header
        param_str += """\nTR,DRW,DRF,DELTA,LATT,NPART,NLP
FE0,RCUT,NLAYER,KH,NPRINT
MIX, EPS1, EPS2, QVOL
KF,R0,LIMOL,YZDIVX
DXXDIV
NTYPE,NON,NGAUS,E0
NPLOT,MAXTAU,EQTIM,NCYC"""
        return(param_str)
    def read_from_file(self, input_file):
        # Read the input file line by line and assign to internal variables
        with open(input_file, "r") as ifp:
            # Manually iterate through the data, since each line needs to go into a different set of
            # variables
            data = next(ifp).split()
            tr, drw, drf, delta, latt, npart, nlp = data

            data = next(ifp).split()
            fe0, rcut, nlayer, kh, nprint = data

            data = next(ifp).split()
            mix, eps1, eps2, qvol = data
            
            #data = next(ifp).split() # can't have these until the input file is found??? template_SSGK
            #eps_LJG, r_LJG = data

            data = next(ifp).split()
            kf, r0, limol, xydivz = data

            dxxdiv = float(next(ifp))

            data = next(ifp).split()
            ntype, non, ngaus, e0 = data

            data = next(ifp).split()
            nplot, maxtau, eqtim, ncyc = data

        # now cast the data we read from the file into the proper numpy types
        self._tr, self._drw, self._drf, self._delta = np.array([tr, drw, drf, delta]).astype(np.double)
        self._latt, self._npart = np.array([latt, npart]).astype(np.int)
        self._fe0, self._rcut, self._kh = np.array([fe0, rcut, kh]).astype(np.double)
        self._nprint = np.int(nprint)
        self._mix, self._eps1, self._eps2, self._qvol = np.array([mix, eps1, eps2, qvol]).astype(np.double)
        #self._eps_LJG, self._r_LJG = np.array([eps_0, r_0]).astype(np.double) # same as above

        self._kf, self._r0, self._limol, self._xydivz = np.array([kf,r0,limol,xydivz]).astype(np.double)
        self._dxxdiv = np.double(dxxdiv)
        self._ntype, self._non, self._ngaus = np.array([ntype,non,ngaus]).astype(np.int)
        self._e0 = np.double(e0)
        self._nplot, self._maxtau, self._eqtim, self._ncyc = np.array([nplot,maxtau,eqtim,ncyc]).astype(np.int)
        self.reset_and_update_parameters()

    def run(self, nsteps):
        # First, advance the simulation
        TTCF.md(nsteps, self._iflag)

        # Now update all of our thermodynamic parameters with values from LAMMPS
        self.get_params_from_TTCF()
