# Input specification for Walls_SSGK - non-equilibrium molecular dynamics code
## General notes
- The input file is plain text, and is read one line at a time via Fortran's 
  formatted input routines. The template file is given by `template_SSGK.in`.
- The executable takes one argument from the command line, which is the
  name and path of the input file to use for this calculation.
    * For example: `./bin/Walls_SSGK template_SSGK.in`
- The data fields are either ints or floats. This is not documented in the
  input file itself (currently hard-coded into the source code), so it's safest
  to make sure that if you're editing a field which is a float in the template
  file that it stays a float in the new file.

- Fields are separated by whitespace (doesn't matter whether it's tabs or 
  spaces).
- The files have 2N lines: the first N are the actual input values, while the
  last N are text descriptions of each line, in corresponding order. For
  example, the first data line has seven elements, which the first line of text
  (line 8, in this case) indicates correspond to the variables 
  `TR,DRW,DRF,DELTA,LATT,NPART,NLP`.

## Equilibrium vs non-equilibrium calculations
- The usual process for running simulations is to start with an equilibrium
  simulation, a process known as "equilibrating the system". This can sometimes
  require very long run-times and will produce several output files in the
  directory where the calc is run.
- After equilibration, we then want to run non-equilibrium dynamics, using the
  output from the equilibrium system as a starting point. To do this, we re-run
  the calculation in the same directory as the original, equilibrium
  calculation, but with different a input file.
- The equilibrium vs non-equilibrium part of the calculation is controlled via
  the variable "NTYPE" in the input file. Currently, this is the first variable
  in line 6 of the input file.
    * NTYPE=1 denotes an equilibrium MD calculation,
    * NTYPE=2 denotes non-equilibrium MD
- I've handled this so far by making two input files, one for NTYPE=1 and one
  for NTYPE=2, and then running the code twice in the same directory. This may
  change though, as it is sometimes more efficient to do one equilibrium
  calculation and use that as the starting point for multiple non-equilibrium
  calculations.

## List of input parameters
There are multiple versions of the code floating around with slightly different
input options. Until I finish unifying the different versions, treat this
section as a work-in-progress.

- Line one:
    * TR (float): Parameterises the velocity distribution
    * DRW (float): Wall density (relative to the number of particles)
    * DRF (float): Fluid density (relative to the number of particles)
    * DELTA (float): Timestep used when integrating equations of motion.
    * LATT (int): Determines whether to arrange wall particles in FCC (LATT=1)
      or BCC (LATT=2) lattice.
    * NPART (int): Total number of particles. Includes both wall and fluid 
      particles. Hardcoded to be <2500 because of a few static arrays (this is
      on the todo list to fix).
      NWALL is calculated using NLAYER and NLP:  NFLUID = NPART - NWALL
    * NLP (int): Number of particles per layer of wall. Used to calculate the 
      number of "wall" particles by the formula: NWALL = NLP*NLAYER
- Line two: 
    * FE0 (float): Initial value of the NE colour-field
    * RCUT (float): Cutoff radius to use when calculating intermolecular forces
    * NLAYER (int): Number of layers in the walls. Used to calculate the number
      of "wall" particles by the formula: NWALL = NLP*NLAYER
    * KH (float): Restoring force for particles in the wall lattice.
    * NPRINT (int): Time increment at which to write particle positions to the
      `.xyz` files (e.g. NPRINT == 100 will write position values every 100 
      time-steps).
- Line three:
    * MIX, EPS1, EPS2 (float): Control fluid mixing during force calculations.
- Line three: 
    * KF (float): Parameter which determines the strength of the intermolecular
      forces between liquid particles.
    * R0 (float): Radius over which the WCA potential acts.
    * LIMOL (int): number of particles per fluid molecule. Every fluid particle
      must be part of a molecule, so nfluid%limol must = 0
    * XYDIVZ (float): Determines the length of the box used in periodic 
      boundary conditions. For each box, x_length = y_length, 
      z_length = x_length/XYDIVZ (note: internally, the box dimensions are 
      stored as CUBEX, CUBEY, CUBEZ).
- Line four: 
    DXXDIV (float): Appears to be currently unused.

- Line five:
    * NTYPE (int): type of calculation to do. 1 = equilibrium, 
      2 = non-equilibrium, 3 = restart from previous run.
    * NON (int): if NON == 0, then we re-scale kinetic energies/temperatures, 
      rather than using a thermostat.
    * NGAUS (int): Control variable to decide which form of Gaussian thermostat
      to use. TODO: Document what each value does.
    * E0 (float): Initial average kinetic energy of particles. Each particle's
      initial kinetic energy is generated from a Gaussian distribution centred
      on E0. This also means that, initially, the sum of squares of momenta 
      equals 2*E0*NPART.
- Final line:
    * NPLOT (int): time increment at which to print instantaneous values during
      MD-loop (e.g. NPLOT == 100 will write values every 100 time-steps)
    * MAXTAU (int): Maximum length of time over which to calculate correlation
      functions. Essentially, this sets the bounds on the discrete sum over
      instantaneous quantities.
    * EQTIM (int): Number of time steps to spend equilibrating the system. This
      "burn-in" is carried out before any statistical properties are
      calculated.
    * NCYC (int): number of statistical trials to carry out.

## Output
- The code produces multiple streams of output, both to stdout and to some
  files in the calculation directory. Both forms out output are necessary:
    * stdout has some consistency checks and instantaneous values of
      thermodynamic quantities of interest.
    * The output files contain information on molecular trajectories, as well
      as ensemble averages of thermodynamics quantities as a function of time.
      They are a mix of plaint text and binary and are useful for analysis and
      visualisation.
- The output files contain useful scientific information, and can also be
  used as a checkpoint to re-start a calculation.
- We want to average the results from multiple runs to get better statistical
  uncertainties. At the moment, it will just be enough to ensure that the
  output files from each run are available *somewhere* on the cluster(s). I
  have some scripts to aggregate the results, but the code is still sort of in
  flux, so this could change in a few weeks.

## Parameter sweeps
Here are the parameters we intend to sweep over for our preliminary
experiments. They may change later, but once we're set up with Nimrod we'll be
able to use preliminary results to guide future parameter sweeps.
- NPART
- RCUT
- FE0
- NON: will sometimes want to do multiple sweeps, one for NON=0 and NON=1
- EQTIM
- MAXTAU
