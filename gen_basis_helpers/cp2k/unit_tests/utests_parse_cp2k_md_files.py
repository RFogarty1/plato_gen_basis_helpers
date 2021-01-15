
import itertools as it
import unittest 
import unittest.mock as mock

import plato_pylib.shared.ucell_class as uCellHelp
import plato_pylib.shared.unit_convs as uConvHelp

import gen_basis_helpers.analyse_md.thermo_data as thermoDataHelp
import gen_basis_helpers.analyse_md.traj_core as trajHelp

import gen_basis_helpers.cp2k.parse_md_files as tCode

class TestParseInfoFromCpout(unittest.TestCase):

	def setUp(self):
		self.createTestObjs()

	def createTestObjs(self):
		self.cpoutFileAsListA = _getCpoutMdNPTFileStrA().split("\n")

	@mock.patch("gen_basis_helpers.cp2k.parse_md_files.parseCP2KHelp._getFileAsListFromInpFile")
	def testExpectedThermoInfoParsed(self, mockedReadFileIntoList):
		mockedReadFileIntoList.side_effect = lambda *args,**kwargs: self.cpoutFileAsListA
		expObj = self._loadExpectedThermalDataA()
		actObj = tCode.parseCpoutForMDJob("fake_path")["thermo_data"]
		mockedReadFileIntoList.assert_called_with("fake_path")
		self.assertEqual(expObj, actObj)

	@mock.patch("gen_basis_helpers.cp2k.parse_md_files.parseCP2KHelp._getFileAsListFromInpFile")
	def testExpectedTrajectoryParsedFileA(self, mockedReadFileIntoList):
		mockedReadFileIntoList.side_effect = lambda *args,**kwargs : self.cpoutFileAsListA
		expObj = self._loadExpectedTrajectoryObjA()
		actObj = tCode.parseCpoutForMDJob("fake_path")["trajectory"]
		mockedReadFileIntoList.assert_called_with("fake_path")
		self.assertEqual(expObj,actObj)

	@unittest.skip("May or may not bother....")
	def testExpectedTrajectoryParsedFile_NVT_a(self):
		self.assertTrue(False)


	def _loadExpectedThermalDataA(self):
		steps = [1,2,3]
		times = [0.5,1.0,1.5]
		temps = [300.206, 300.323, 300.349]
		eKinetic = [uConvHelp.RYD_TO_EV*2*x for x in [0.142604641825E-02, 0.142659934591E-02, 0.142672695109E-02]]
		ePot = [uConvHelp.RYD_TO_EV*2*x for x in [-0.176102057981E+01, -0.176101613342E+01, -0.176101115757E+01]]
		pressure = [0.806499998740E+04, 0.829617125241E+04, 0.852673138917E+04]

		outKwargs = {"step":steps, "time":times, "temp":temps,
		             "eKinetic":eKinetic, "ePot":ePot, "pressure":pressure}
		outObj = thermoDataHelp.ThermoDataStandard.fromStdKwargs(**outKwargs)
		return outObj


	def _loadExpectedTrajectoryObjA(self):
		cellA = uCellHelp.UnitCell(lattParams=[6.038771, 6.038771, 9.768012],lattAngles=[90,90,120])
		cellB = uCellHelp.UnitCell(lattParams=[6.037550, 6.037550, 9.766037],lattAngles=[90,90,120])
		cellC = uCellHelp.UnitCell(lattParams=[6.036338, 6.036338, 9.764076],lattAngles=[90,90,120])
		steps = [1,2,3]
		times = [0.5,1.0,1.5]
		cells = [cellA, cellB, cellC]
		trajSteps = [trajHelp.TrajStepBase(unitCell=cell, step=step, time=time) for cell,step,time in it.zip_longest(cells, steps, times)]
		outObj = trajHelp.TrajectoryInMemory(trajSteps)
		return outObj

class TestParseInitMDThermoSection(unittest.TestCase):

	def setUp(self):
		self.createTestObjs()

	def createTestObjs(self):
		self.fileAsListA = self._loadFileStrA().split("\n")
		self.startIdxA = 1

	def _loadFileStrA(self):
		return """
 ******************************** GO CP2K GO! **********************************
 INITIAL POTENTIAL ENERGY[hartree]     =                     -0.176102449402E+01
 INITIAL KINETIC ENERGY[hartree]       =                      0.142506690448E-02
 INITIAL TEMPERATURE[K]                =                                 300.000
 INITIAL BAROSTAT TEMP[K]              =                      0.300000000000E+03
 INITIAL PRESSURE[bar]                 =                      0.783325915339E+04
 INITIAL VOLUME[bohr^3]                =                      0.308673305529E+03
 INITIAL CELL LNTHS[bohr]   =      0.6040000E+01   0.6040000E+01   0.9770000E+01
 INITIAL CELL ANGLS[deg]    =      0.9000000E+02   0.9000000E+02   0.1200000E+03
 ******************************** GO CP2K GO! **********************************
"""

	def testExpectedSectionA(self):
		expEndIdx = 11
		expDict = {"lattAngles":[90,90,120], "lattParams":[x*uConvHelp.BOHR_TO_ANG for x in [6.04,6.04,9.77]]}
		actDict, actEndIdx = tCode._parseMDInitSection(self.fileAsListA, self.startIdxA)
		expCell = uCellHelp.UnitCell(lattParams=expDict["lattParams"], lattAngles=expDict["lattAngles"])
		actCell = uCellHelp.UnitCell(lattParams=actDict["lattParams"], lattAngles=actDict["lattAngles"])
		self.assertEqual(expEndIdx, actEndIdx)
		self.assertEqual(expCell,actCell)


class TestParseMdStepGenInfo(unittest.TestCase):

	def setUp(self):
		self.createTestObjs()

	def createTestObjs(self):
		self.fileAsListNPT_a = self._loadFileStrNPT_a().split("\n")
		self.fileAsListNVT_a = self._loadFileStrNVT_a().split("\n")
		self.startIdxNPT_a = 2
		self.startIdxNVT_a = 2

	def _loadFileStrNPT_a(self):
		return """ 
 *******************************************************************************
 ENSEMBLE TYPE                =                                            NPT_I
 STEP NUMBER                  =                                                2
 TIME [fs]                    =                                         1.000000
 CONSERVED QUANTITY [hartree] =                              -0.175817330894E+01

                                              INSTANTANEOUS             AVERAGES
 CPU TIME [s]                 =                       28.43                42.06
 ENERGY DRIFT PER ATOM [K]    =          0.333491073427E-03   0.166745536714E-03
 POTENTIAL ENERGY[hartree]    =         -0.176101613342E+01  -0.176101835662E+01
 KINETIC ENERGY [hartree]     =          0.142659934591E-02   0.142632288208E-02
 TEMPERATURE [K]              =                     300.323              300.264
 PRESSURE [bar]               =          0.829617125241E+04   0.818058561990E+04
 BAROSTAT TEMP[K]             =          0.292352496174E+03   0.294280744292E+03
 VOLUME[bohr^3]               =          0.308297884123E+03   0.308391407993E+03
 CELL LNTHS[bohr]             =    0.6037550E+01   0.6037550E+01   0.9766037E+01
 AVE. CELL LNTHS[bohr]        =    0.6038161E+01   0.6038161E+01   0.9767025E+01
 *******************************************************************************
"""

	def _loadFileStrNVT_a(self):
		return """
 *******************************************************************************
 ENSEMBLE TYPE                =                                              NVT
 STEP NUMBER                  =                                                5
 TIME [fs]                    =                                         2.500000
 CONSERVED QUANTITY [hartree] =                              -0.175912673291E+01

                                              INSTANTANEOUS             AVERAGES
 CPU TIME [s]                 =                       10.18                33.08
 ENERGY DRIFT PER ATOM [K]    =         -0.367576766846E+00  -0.146931593261E+00
 POTENTIAL ENERGY[hartree]    =         -0.176101943050E+01  -0.176102275265E+01
 KINETIC ENERGY [hartree]     =          0.141612684746E-02   0.142099809148E-02
 TEMPERATURE [K]              =                     298.118              299.143
 PRESSURE [bar]               =          0.783741673349E+04   0.783387164598E+04
 *******************************************************************************
"""

	def testExpectedDictNPT_a(self):
		expEndIdx = 19
		haToEv = uConvHelp.RYD_TO_EV*2
		expDict = {"step": 2, "time":1.0, "temp": 300.323,
		           "eKinetic": 0.142659934591E-02*haToEv,
		           "ePot": -0.176101613342E+01*haToEv,
		           "pressure": 0.829617125241E+04,
		           "lattParams": [0.6037550E+01,0.6037550E+01,0.9766037E+01] }
		actDict, actEndIdx = tCode._parseMdStepInfo(self.fileAsListNPT_a, self.startIdxNPT_a)
		self.assertEqual(expEndIdx, actEndIdx)
		for key in expDict.keys():
			if key!="lattParams":
				self.assertAlmostEqual( expDict[key], actDict[key] )
			else:
				[self.assertAlmostEqual(exp,act) for exp,act in it.zip_longest(expDict["lattParams"],actDict["lattParams"])]

	def testExpectedDictNVT_a(self):
		expEndIdx = 15
		haToEv = uConvHelp.RYD_TO_EV*2
		expDict = {"step": 5, "time": 2.5, "ePot": -0.176101943050E+01*haToEv,
		           "eKinetic": 0.141612684746E-02*haToEv, "temp": 298.118, "pressure": 0.783741673349E+04}
		actDict, actEndIdx = tCode._parseMdStepInfo(self.fileAsListNVT_a, self.startIdxNVT_a)
		self.assertEqual(expEndIdx, actEndIdx)
		for key in expDict.keys():
			self.assertAlmostEqual( expDict[key], actDict[key] )


def _getXyzMdNPTFileStrA():
	return """
       2
 i =        0, time =        0.000, E =        -1.7610244940
 Mg         0.0000000000        0.0000000000        0.0000000000
 Mg        -0.0000000160        1.8453444568        2.5850306640
       2
 i =        1, time =        0.500, E =        -1.7610205798
 Mg         0.0008159084       -0.0017784971        0.0001434100
 Mg        -0.0008159242        1.8467474920        2.5843612925
       2
 i =        2, time =        1.000, E =        -1.7610161334
 Mg         0.0016322252       -0.0035569587        0.0002867308
 Mg        -0.0016322408        1.8481529833        2.5836955004
       2
 i =        3, time =        1.500, E =        -1.7610111576
 Mg         0.0024488239       -0.0053351324        0.0004298829
 Mg        -0.0024488394        1.8495607336        2.5830334449
"""


def _getCpoutMdNPTFileStrA():
	return """
 DBCSR| Multiplication driver                                               BLAS
 DBCSR| Multrec recursion limit                                              512
 DBCSR| Multiplication stack size                                           1000
 DBCSR| Maximum elements for images                                    UNLIMITED
 DBCSR| Multiplicative factor virtual images                                   1
 DBCSR| Multiplication size stacks                                             3


  **** **** ******  **  PROGRAM STARTED AT               2020-12-17 16:44:29.908
 ***** ** ***  *** **   PROGRAM STARTED ON                             archlinux
 **    ****   ******    PROGRAM STARTED BY                               richard
 ***** **    ** ** **   PROGRAM PROCESS ID                                 23177
  **** **  *******  **  PROGRAM STARTED IN /home/richard/work/Random_Stuff/temp/
                                           temp2/sept_2019_temp_work/cp2k_test_m
                                           d/temp/npt/print_cell

 CP2K| version string:                                          CP2K version 6.1
 CP2K| source code revision number:                                    svn:18464
 CP2K| cp2kflags: fftw3 libxc max_contr=4                                       
 CP2K| is freely available from                            https://www.cp2k.org/
 CP2K| Program compiled at                          Tue Nov  3 16:31:14 GMT 2020
 CP2K| Program compiled on                                             archlinux
 CP2K| Program compiled for                                Linux-x86-64-gfortran
 CP2K| Data directory path    /home/richard/work/Random_Stuff/cp2k/cp2k-6.1.0/da
 CP2K| Input file name                                              geom_opt.inp

 GLOBAL| Force Environment number                                              1
 GLOBAL| Basis set file name                                         PLATO_BASIS
 GLOBAL| Potential file name                                      GTH_POTENTIALS
 GLOBAL| MM Potential file name                                     MM_POTENTIAL
 GLOBAL| Coordinate file name                                      __STD_INPUT__
 GLOBAL| Method name                                                        CP2K
 GLOBAL| Project name                                                   geom_opt
 GLOBAL| Preferred FFT library                                             FFTW3
 GLOBAL| Preferred diagonalization lib.                                       SL
 GLOBAL| Run type                                                             MD
 GLOBAL| All-to-all communication in single precision                          F
 GLOBAL| FFTs using library dependent lengths                                  F
 GLOBAL| Global print level                                                  LOW
 GLOBAL| Total number of message passing processes                             1
 GLOBAL| Number of threads for this process                                    1
 GLOBAL| This output is from process                                           0
 GLOBAL| CPU model name :  Intel(R) Core(TM) i7-7700HQ CPU @ 2.80GHz

 MEMORY| system memory details [Kb]
 MEMORY|                        rank 0           min           max       average
 MEMORY| MemTotal             16141344      16141344      16141344      16141344
 MEMORY| MemFree               8845756       8845756       8845756       8845756
 MEMORY| Buffers                249696        249696        249696        249696
 MEMORY| Cached                3072828       3072828       3072828       3072828
 MEMORY| Slab                   271824        271824        271824        271824
 MEMORY| SReclaimable           138808        138808        138808        138808
 MEMORY| MemLikelyFree        12307088      12307088      12307088      12307088


 GENERATE|  Preliminary Number of Bonds generated:                             0
 GENERATE|  Achieved consistency in connectivity generation.

 *** WARNING in cryssym.F:163 :: Symmetry library SPGLIB not available ***


 *******************************************************************************
 *******************************************************************************
 **                                                                           **
 **     #####                         ##              ##                      **
 **    ##   ##            ##          ##              ##                      **
 **   ##     ##                       ##            ######                    **
 **   ##     ##  ##   ##  ##   #####  ##  ##   ####   ##    #####    #####    **
 **   ##     ##  ##   ##  ##  ##      ## ##   ##      ##   ##   ##  ##   ##   **
 **   ##  ## ##  ##   ##  ##  ##      ####     ###    ##   ######   ######    **
 **    ##  ###   ##   ##  ##  ##      ## ##      ##   ##   ##       ##        **
 **     #######   #####   ##   #####  ##  ##  ####    ##    #####   ##        **
 **           ##                                                    ##        **
 **                                                                           **
 **                                                ... make the atoms dance   **
 **                                                                           **
 **            Copyright (C) by CP2K developers group (2000 - 2018)           **
 **                                                                           **
 *******************************************************************************


 TOTAL NUMBERS AND MAXIMUM NUMBERS

  Total number of            - Atomic kinds:                                   1
                             - Atoms:                                          2
                             - Shell sets:                                     6
                             - Shells:                                        10
                             - Primitive Cartesian functions:                 24
                             - Cartesian basis functions:                     24
                             - Spherical basis functions:                     22

  Maximum angular momentum of- Orbital basis functions:                        2
                             - Local part of the GTH pseudopotential:          0
                             - Non-local part of the GTH pseudopotential:      2


 SCF PARAMETERS         Density guess:                                    ATOMIC
                        --------------------------------------------------------
                        max_scf:                                             300
                        max_scf_history:                                       0
                        max_diis:                                              4
                        --------------------------------------------------------
                        eps_scf:                                        1.00E-06
                        eps_scf_history:                                0.00E+00
                        eps_diis:                                       1.00E-01
                        eps_eigval:                                     1.00E-05
                        --------------------------------------------------------
                        level_shift [a.u.]:                                 0.00
                        added MOs                                        10    0
                        --------------------------------------------------------
                        Mixing method:                            BROYDEN_MIXING
                                                charge density mixing in g-space
                        --------------------------------------------------------
                        Smear method:                                FERMI_DIRAC
                        Electronic temperature [K]:                        157.9
                        Electronic temperature [a.u.]:                  5.00E-04
                        Accuracy threshold:                             1.00E-10
                        --------------------------------------------------------
                        No outer SCF

 MD| Molecular Dynamics Protocol 
 MD| Ensemble Type                                                         NPT_I
 MD| Number of Time Steps                                                      3
 MD| Time Step [fs]                                                         0.50
 MD| Temperature [K]                                                      300.00
 MD| Temperature tolerance [K]                                              0.00
 MD| Pressure [Bar]                                                         1.00
 MD| Barostat time constant [  fs]                                       1000.00
 MD| Print MD information every                                        1 step(s)
 MD| File type     Print frequency[steps]                             File names
 MD| Coordinates            1                                 geom_opt-pos-1.xyz
 MD| Simulation Cel         1                                    geom_opt-1.cell
 MD| Velocities             1                                 geom_opt-vel-1.xyz
 MD| Energies               1                                    geom_opt-1.ener
 MD| Dump                  20                                 geom_opt-1.restart

 ROT| Rotational Analysis Info 
 ROT| Principal axes and moments of inertia in atomic units:
 ROT|                                1                 2                 3
 ROT| EIGENVALUES            0.291038305E-10   0.798021796E+06   0.798021796E+06
 ROT|      X                     0.000000005      -0.000000007       1.000000000
 ROT|      Y                    -0.581007586       0.813898142       0.000000009
 ROT|      Z                    -0.813898142      -0.581007586       0.000000000
 ROT| Numer of Rotovibrational vectors:     5
 ROT| Linear Molecule detected..

 Calculation of degrees of freedom
                                                      Number of atoms:         2
                                 Number of Intramolecular constraints:         0
                                 Number of Intermolecular constraints:         0
                                  Invariants(translation + rotations):         3
                                                   Degrees of freedom:         3


 Restraints Information
                                  Number of Intramolecular restraints:         0
                                  Number of Intermolecular restraints:         0

 THERMOSTAT| Thermostat Info for PARTICLES
 THERMOSTAT| Type of thermostat                               Nose-Hoover-Chains
 THERMOSTAT| Nose-Hoover-Chain length                                          3
 THERMOSTAT| Nose-Hoover-Chain time constant [  fs]                      1000.00
 THERMOSTAT| Order of Yoshida integrator                                       3
 THERMOSTAT| Number of multiple time steps                                     2
 THERMOSTAT| Initial Potential Energy                                   0.000000
 THERMOSTAT| Initial Kinetic Energy                                     0.000475
 THERMOSTAT| End of Thermostat Info for PARTICLES


 THERMOSTAT| Thermostat Info for BAROSTAT
 THERMOSTAT| Type of thermostat                               Nose-Hoover-Chains
 THERMOSTAT| Nose-Hoover-Chain length                                          3
 THERMOSTAT| Nose-Hoover-Chain time constant [  fs]                      1000.00
 THERMOSTAT| Order of Yoshida integrator                                       3
 THERMOSTAT| Number of multiple time steps                                     2
 THERMOSTAT| Initial Potential Energy                                   0.000000
 THERMOSTAT| Initial Kinetic Energy                                     0.000475
 THERMOSTAT| End of Thermostat Info for BAROSTAT

 ************************** Velocities initialization **************************
 Initial Temperature                                                    300.00 K
 COM velocity:            0.000000000000     -0.000000000000      0.000000000000
 *******************************************************************************


 Number of electrons:                                                          4
 Number of occupied orbitals:                                                  2
 Number of molecular orbitals:                                                12

 Number of orbital functions:                                                 22
 Number of independent orbital functions:                                     22

 Extrapolation method: initial_guess


 SCF WAVEFUNCTION OPTIMIZATION

  Step     Update method      Time    Convergence         Total energy    Change
  ------------------------------------------------------------------------------
     1 NoMix/Diag. 0.40E+00    1.5     0.44539833        -1.8870307133 -1.89E+00
     2 Broy./Diag. 0.40E+00    2.9     0.00894568        -1.8551873484  3.18E-02
     3 Broy./Diag. 0.40E+00    2.8     0.02938061        -1.7621492041  9.30E-02
     4 Broy./Diag. 0.40E+00    2.8     0.00031024        -1.7621726963 -2.35E-05
     5 Broy./Diag. 0.40E+00    2.8     0.00005930        -1.7610216019  1.15E-03
     6 Broy./Diag. 0.40E+00    2.9     0.00000379        -1.7610258574 -4.26E-06
     7 Broy./Diag. 0.40E+00    3.0     0.00001301        -1.7610242101  1.65E-06
     8 Broy./Diag. 0.40E+00    3.0     0.00000029        -1.7610244775 -2.67E-07

  *** SCF run converged in     8 steps ***


  Electronic density on regular grids:         -4.0000000000       -0.0000000000
  Core density on regular grids:                4.0000000000       -0.0000000000
  Total charge density on r-space grids:       -0.0000000001
  Total charge density g-space grids:          -0.0000000001

  Overlap energy of the core charge distribution:               0.00000000000129
  Self energy of the core charge distribution:                 -3.91146295972394
  Core Hamiltonian energy:                                      0.96207479510633
  Hartree energy:                                               2.04889947620939
  Exchange-correlation energy:                                 -0.86049207451313
  Electronic entropic energy:                                  -0.00004371461430
  Fermi energy:                                                 0.10433249107928

  Total energy:                                                -1.76102447753185

 ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):               -1.761024494023267


 ATOMIC FORCES in [a.u.]

 # Atom   Kind   Element          X              Y              Z
      1      1      Mg          0.00011597    -0.00006695     0.00000000
      2      1      Mg         -0.00011591     0.00006692     0.00000000
 SUM OF ATOMIC FORCES           0.00000006    -0.00000004     0.00000000     0.00000007

 MD_ENERGIES| Initialization proceeding


 ******************************** GO CP2K GO! **********************************
 INITIAL POTENTIAL ENERGY[hartree]     =                     -0.176102449402E+01
 INITIAL KINETIC ENERGY[hartree]       =                      0.142506690448E-02
 INITIAL TEMPERATURE[K]                =                                 300.000
 INITIAL BAROSTAT TEMP[K]              =                      0.300000000000E+03
 INITIAL PRESSURE[bar]                 =                      0.783325915339E+04
 INITIAL VOLUME[bohr^3]                =                      0.308673305529E+03
 INITIAL CELL LNTHS[bohr]   =      0.6040000E+01   0.6040000E+01   0.9770000E+01
 INITIAL CELL ANGLS[deg]    =      0.9000000E+02   0.9000000E+02   0.1200000E+03
 ******************************** GO CP2K GO! **********************************

 Number of electrons:                                                          4
 Number of occupied orbitals:                                                  2
 Number of molecular orbitals:                                                12

 Number of orbital functions:                                                 22
 Number of independent orbital functions:                                     22

 Extrapolation method: initial_guess


 SCF WAVEFUNCTION OPTIMIZATION

  Step     Update method      Time    Convergence         Total energy    Change
  ------------------------------------------------------------------------------
     1 Broy./Diag. 0.40E+00    1.5     0.44541236        -1.8872542317 -1.89E+00
     2 Broy./Diag. 0.40E+00    2.9     0.00895579        -1.8551447233  3.21E-02
     3 Broy./Diag. 0.40E+00    2.9     0.02942643        -1.7621441458  9.30E-02
     4 Broy./Diag. 0.40E+00    2.8     0.00031119        -1.7621675328 -2.34E-05
     5 Broy./Diag. 0.40E+00    2.8     0.00005787        -1.7610177918  1.15E-03
     6 Broy./Diag. 0.40E+00    2.9     0.00000313        -1.7610220057 -4.21E-06
     7 Broy./Diag. 0.40E+00    2.9     0.00001196        -1.7610203929  1.61E-06
     8 Broy./Diag. 0.40E+00    3.1     0.00000135        -1.7610204753 -8.23E-08
     9 Broy./Diag. 0.40E+00    3.0     0.00000002        -1.7610205858 -1.11E-07

  *** SCF run converged in     9 steps ***


  Electronic density on regular grids:         -4.0000000000       -0.0000000000
  Core density on regular grids:                4.0000000000       -0.0000000000
  Total charge density on r-space grids:       -0.0000000001
  Total charge density g-space grids:          -0.0000000001

  Overlap energy of the core charge distribution:               0.00000000000130
  Self energy of the core charge distribution:                 -3.91146295972394
  Core Hamiltonian energy:                                      0.96245657246016
  Hartree energy:                                               2.04868706862421
  Exchange-correlation energy:                                 -0.86065757912155
  Electronic entropic energy:                                  -0.00004368804810
  Fermi energy:                                                 0.10445299115967

  Total energy:                                                -1.76102058580888

 ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):               -1.761020579812398


 ATOMIC FORCES in [a.u.]

 # Atom   Kind   Element          X              Y              Z
      1      1      Mg          0.00009101    -0.00001697    -0.00001555
      2      1      Mg         -0.00009101     0.00001697     0.00001555
 SUM OF ATOMIC FORCES           0.00000000    -0.00000000    -0.00000000     0.00000000

 *******************************************************************************
 ENSEMBLE TYPE                =                                            NPT_I
 STEP NUMBER                  =                                                1
 TIME [fs]                    =                                         0.500000
 CONSERVED QUANTITY [hartree] =                              -0.175817330912E+01

                                              INSTANTANEOUS             AVERAGES
 CPU TIME [s]                 =                       55.69                55.69
 ENERGY DRIFT PER ATOM [K]    =          0.305765887366E-03   0.000000000000E+00
 POTENTIAL ENERGY[hartree]    =         -0.176102057981E+01  -0.176102057981E+01
 KINETIC ENERGY [hartree]     =          0.142604641825E-02   0.142604641825E-02
 TEMPERATURE [K]              =                     300.206              300.206
 PRESSURE [bar]               =          0.806499998740E+04   0.806499998740E+04
 BAROSTAT TEMP[K]             =          0.296208992411E+03   0.296208992411E+03
 VOLUME[bohr^3]               =          0.308484931864E+03   0.308484931864E+03
 CELL LNTHS[bohr]             =    0.6038771E+01   0.6038771E+01   0.9768012E+01
 AVE. CELL LNTHS[bohr]        =    0.6038771E+01   0.6038771E+01   0.9768012E+01
 *******************************************************************************


 Number of electrons:                                                          4
 Number of occupied orbitals:                                                  2
 Number of molecular orbitals:                                                12

 Number of orbital functions:                                                 22
 Number of independent orbital functions:                                     22

 Extrapolation method: initial_guess


 SCF WAVEFUNCTION OPTIMIZATION

  Step     Update method      Time    Convergence         Total energy    Change
  ------------------------------------------------------------------------------
     1 Broy./Diag. 0.40E+00    1.4     0.44542373        -1.8874326600 -1.89E+00
     2 Broy./Diag. 0.40E+00    2.8     0.00896587        -1.8551018552  3.23E-02
     3 Broy./Diag. 0.40E+00    2.8     0.02947208        -1.7621386542  9.30E-02
     4 Broy./Diag. 0.40E+00    2.8     0.00031219        -1.7621618041 -2.31E-05
     5 Broy./Diag. 0.40E+00    2.8     0.00005764        -1.7610137268  1.15E-03
     6 Broy./Diag. 0.40E+00    2.8     0.00000315        -1.7610177379 -4.01E-06
     7 Broy./Diag. 0.40E+00    2.8     0.00001014        -1.7610161102  1.63E-06
     8 Broy./Diag. 0.40E+00    2.8     0.00000313        -1.7610159554  1.55E-07
     9 Broy./Diag. 0.40E+00    2.8     0.00000001        -1.7610161384 -1.83E-07

  *** SCF run converged in     9 steps ***


  Electronic density on regular grids:         -4.0000000000       -0.0000000000
  Core density on regular grids:                4.0000000000       -0.0000000000
  Total charge density on r-space grids:       -0.0000000001
  Total charge density g-space grids:          -0.0000000001

  Overlap energy of the core charge distribution:               0.00000000000132
  Self energy of the core charge distribution:                 -3.91146295972394
  Core Hamiltonian energy:                                      0.96283686745205
  Hartree energy:                                               2.04847568180867
  Exchange-correlation energy:                                 -0.86082207040541
  Electronic entropic energy:                                  -0.00004365749066
  Fermi energy:                                                 0.10457291151314

  Total energy:                                                -1.76101613835888

 ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):               -1.761016133419581


 ATOMIC FORCES in [a.u.]

 # Atom   Kind   Element          X              Y              Z
      1      1      Mg          0.00006622     0.00003251    -0.00003115
      2      1      Mg         -0.00006622    -0.00003251     0.00003115
 SUM OF ATOMIC FORCES           0.00000000    -0.00000000    -0.00000000     0.00000000

 *******************************************************************************
 ENSEMBLE TYPE                =                                            NPT_I
 STEP NUMBER                  =                                                2
 TIME [fs]                    =                                         1.000000
 CONSERVED QUANTITY [hartree] =                              -0.175817330894E+01

                                              INSTANTANEOUS             AVERAGES
 CPU TIME [s]                 =                       28.43                42.06
 ENERGY DRIFT PER ATOM [K]    =          0.333491073427E-03   0.166745536714E-03
 POTENTIAL ENERGY[hartree]    =         -0.176101613342E+01  -0.176101835662E+01
 KINETIC ENERGY [hartree]     =          0.142659934591E-02   0.142632288208E-02
 TEMPERATURE [K]              =                     300.323              300.264
 PRESSURE [bar]               =          0.829617125241E+04   0.818058561990E+04
 BAROSTAT TEMP[K]             =          0.292352496174E+03   0.294280744292E+03
 VOLUME[bohr^3]               =          0.308297884123E+03   0.308391407993E+03
 CELL LNTHS[bohr]             =    0.6037550E+01   0.6037550E+01   0.9766037E+01
 AVE. CELL LNTHS[bohr]        =    0.6038161E+01   0.6038161E+01   0.9767025E+01
 *******************************************************************************


 Number of electrons:                                                          4
 Number of occupied orbitals:                                                  2
 Number of molecular orbitals:                                                12

 Number of orbital functions:                                                 22
 Number of independent orbital functions:                                     22

 Extrapolation method: initial_guess


 SCF WAVEFUNCTION OPTIMIZATION

  Step     Update method      Time    Convergence         Total energy    Change
  ------------------------------------------------------------------------------
     1 Broy./Diag. 0.40E+00    1.5     0.44543251        -1.8876097064 -1.89E+00
     2 Broy./Diag. 0.40E+00    2.8     0.00897591        -1.8550587538  3.26E-02
     3 Broy./Diag. 0.40E+00    2.8     0.02951758        -1.7621327317  9.29E-02
     4 Broy./Diag. 0.40E+00    2.8     0.00031324        -1.7621555126 -2.28E-05
     5 Broy./Diag. 0.40E+00    2.8     0.00005739        -1.7610094109  1.15E-03
     6 Broy./Diag. 0.40E+00    2.8     0.00000320        -1.7610130397 -3.63E-06
     7 Broy./Diag. 0.40E+00    2.8     0.00000891        -1.7610112434  1.80E-06
     8 Broy./Diag. 0.40E+00    2.8     0.00000432        -1.7610109342  3.09E-07
     9 Broy./Diag. 0.40E+00    2.8     6.8258E-09        -1.7610111619 -2.28E-07

  *** SCF run converged in     9 steps ***


  Electronic density on regular grids:         -4.0000000000       -0.0000000000
  Core density on regular grids:                4.0000000000       -0.0000000000
  Total charge density on r-space grids:       -0.0000000001
  Total charge density g-space grids:          -0.0000000001

  Overlap energy of the core charge distribution:               0.00000000000134
  Self energy of the core charge distribution:                 -3.91146295972394
  Core Hamiltonian energy:                                      0.96321560435026
  Hartree energy:                                               2.04826534423863
  Exchange-correlation energy:                                 -0.86098552791690
  Electronic entropic energy:                                  -0.00004362286393
  Fermi energy:                                                 0.10469222959823

  Total energy:                                                -1.76101116191538

 ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):               -1.761011157570904


 ATOMIC FORCES in [a.u.]

 # Atom   Kind   Element          X              Y              Z
      1      1      Mg          0.00004099     0.00008182    -0.00004680
      2      1      Mg         -0.00004099    -0.00008182     0.00004680
 SUM OF ATOMIC FORCES           0.00000000    -0.00000000    -0.00000000     0.00000000

 *******************************************************************************
 ENSEMBLE TYPE                =                                            NPT_I
 STEP NUMBER                  =                                                3
 TIME [fs]                    =                                         1.500000
 CONSERVED QUANTITY [hartree] =                              -0.175817330875E+01

                                              INSTANTANEOUS             AVERAGES
 CPU TIME [s]                 =                       28.38                37.50
 ENERGY DRIFT PER ATOM [K]    =          0.363596349041E-03   0.232362474156E-03
 POTENTIAL ENERGY[hartree]    =         -0.176101115757E+01  -0.176101595693E+01
 KINETIC ENERGY [hartree]     =          0.142672695109E-02   0.142645757175E-02
 TEMPERATURE [K]              =                     300.349              300.293
 PRESSURE [bar]               =          0.852673138917E+04   0.829596754299E+04
 BAROSTAT TEMP[K]             =          0.288432751623E+03   0.292331413403E+03
 VOLUME[bohr^3]               =          0.308112188496E+03   0.308298334828E+03
 CELL LNTHS[bohr]             =    0.6036338E+01   0.6036338E+01   0.9764076E+01
 AVE. CELL LNTHS[bohr]        =    0.6037553E+01   0.6037553E+01   0.9766042E+01
 *******************************************************************************


 -------------------------------------------------------------------------------
 -                                                                             -
 -                                DBCSR STATISTICS                             -
 -                                                                             -
 -------------------------------------------------------------------------------
 COUNTER                                    TOTAL       BLAS       SMM       ACC
 flops inhomo. stacks                           0       0.0%      0.0%      0.0%
 flops total                         0.000000E+00       0.0%      0.0%      0.0%
 flops max/rank                      0.000000E+00       0.0%      0.0%      0.0%
 matmuls inhomo. stacks                         0       0.0%      0.0%      0.0%
 matmuls total                                  0       0.0%      0.0%      0.0%
 number of processed stacks                     0       0.0%      0.0%      0.0%
 average stack size                                     0.0       0.0       0.0
 marketing flops                     0.000000E+00
 -------------------------------------------------------------------------------

 MEMORY| Estimated peak process memory [MiB]                                  79

 -------------------------------------------------------------------------------
 -                                                                             -
 -                           R E F E R E N C E S                               -
 -                                                                             -
 -------------------------------------------------------------------------------
 
 CP2K version 6.1, the CP2K developers group (2018).
 CP2K is freely available from https://www.cp2k.org/ .

 Schuett, Ole; Messmer, Peter; Hutter, Juerg; VandeVondele, Joost. 
 Electronic Structure Calculations on Graphics Processing Units, John Wiley & Sons, Ltd, 173-190 (2016). 
 GPU-Accelerated Sparse Matrix-Matrix Multiplication for Linear Scaling Density Functional Theory.
 http://dx.doi.org/10.1002/9781118670712.ch8


 Borstnik, U; VandeVondele, J; Weber, V; Hutter, J. 
 PARALLEL COMPUTING, 40 (5-6), 47-58 (2014). 
 Sparse matrix multiplication: The distributed block-compressed sparse
 row library.
 http://dx.doi.org/10.1016/j.parco.2014.03.012


 Hutter, J; Iannuzzi, M; Schiffmann, F; VandeVondele, J. 
 WILEY INTERDISCIPLINARY REVIEWS-COMPUTATIONAL MOLECULAR SCIENCE, 4 (1), 15-25 (2014). 
 CP2K: atomistic simulations of condensed matter systems.
 http://dx.doi.org/10.1002/wcms.1159


 Krack, M. 
 THEORETICAL CHEMISTRY ACCOUNTS, 114 (1-3), 145-152 (2005). 
 Pseudopotentials for H to Kr optimized for gradient-corrected
 exchange-correlation functionals.
 http://dx.doi.org/10.1007/s00214-005-0655-y


 VandeVondele, J; Krack, M; Mohamed, F; Parrinello, M; Chassaing, T;
 Hutter, J. COMPUTER PHYSICS COMMUNICATIONS, 167 (2), 103-128 (2005). 
 QUICKSTEP: Fast and accurate density functional calculations using a
 mixed Gaussian and plane waves approach.
 http://dx.doi.org/10.1016/j.cpc.2004.12.014


 Frigo, M; Johnson, SG. 
 PROCEEDINGS OF THE IEEE, 93 (2), 216-231 (2005). 
 The design and implementation of FFTW3.
 http://dx.doi.org/10.1109/JPROC.2004.840301


 Hartwigsen, C; Goedecker, S; Hutter, J. 
 PHYSICAL REVIEW B, 58 (7), 3641-3662 (1998). 
 Relativistic separable dual-space Gaussian pseudopotentials from H to Rn.
 http://dx.doi.org/10.1103/PhysRevB.58.3641


 Lippert, G; Hutter, J; Parrinello, M. 
 MOLECULAR PHYSICS, 92 (3), 477-487 (1997). 
 A hybrid Gaussian and plane wave density functional scheme.
 http://dx.doi.org/10.1080/002689797170220


 Perdew, JP; Burke, K; Ernzerhof, M. 
 PHYSICAL REVIEW LETTERS, 77 (18), 3865-3868 (1996). 
 Generalized gradient approximation made simple.
 http://dx.doi.org/10.1103/PhysRevLett.77.3865


 Goedecker, S; Teter, M; Hutter, J. 
 PHYSICAL REVIEW B, 54 (3), 1703-1710 (1996). 
 Separable dual-space Gaussian pseudopotentials.
 http://dx.doi.org/10.1103/PhysRevB.54.1703


 NOSE, S. JOURNAL OF CHEMICAL PHYSICS, 81 (1), 511-519 (1984). 
 A UNIFIED FORMULATION OF THE CONSTANT TEMPERATURE MOLECULAR-DYNAMICS
 METHODS.
 http://dx.doi.org/10.1063/1.447334


 NOSE, S. MOLECULAR PHYSICS, 52 (2), 255-268 (1984). 
 A MOLECULAR-DYNAMICS METHOD FOR SIMULATIONS IN THE CANONICAL ENSEMBLE.
 http://dx.doi.org/10.1080/00268978400101201


 -------------------------------------------------------------------------------
 -                                                                             -
 -                                T I M I N G                                  -
 -                                                                             -
 -------------------------------------------------------------------------------
 SUBROUTINE                       CALLS  ASD         SELF TIME        TOTAL TIME
                                MAXIMUM       AVERAGE  MAXIMUM  AVERAGE  MAXIMUM
 CP2K                                 1  1.0    0.017    0.017  112.554  112.554
 qs_mol_dyn_low                       1  2.0    0.001    0.001  112.510  112.510
 qs_forces                            4  3.8    0.001    0.001  112.503  112.503
 qs_energies                          4  4.8    0.006    0.006  102.522  102.522
 scf_env_do_scf                       4  5.8    0.000    0.000  100.476  100.476
 scf_env_do_scf_inner_loop           35  6.8    0.010    0.010  100.476  100.476
 velocity_verlet                      3  3.0    0.000    0.000   86.316   86.316
 rebuild_ks_matrix                   39  8.5    0.000    0.000   51.413   51.413
 qs_ks_build_kohn_sham_matrix        39  9.5    0.019    0.019   51.413   51.413
 sum_up_and_integrate                39 10.5    0.000    0.000   51.255   51.255
 integrate_v_rspace                  39 11.5   50.787   50.787   51.254   51.254
 qs_rho_update_rho                   39  7.8    0.000    0.000   48.955   48.955
 calculate_rho_elec                  39  8.8   48.941   48.941   48.955   48.955
 qs_ks_update_qs_env                 35  7.8    0.000    0.000   43.544   43.544
 qs_scf_new_mos_kp                   35  7.8    0.000    0.000    7.968    7.968
 do_general_diag_kp                  35  8.8    0.034    0.034    7.968    7.968
 qs_ks_update_qs_env_forces           4  4.8    0.000    0.000    7.938    7.938
 rskp_transform                    5040  9.8    4.340    4.340    4.340    4.340
 kpoint_density_transform            39  9.5    0.019    0.019    2.391    2.391
 -------------------------------------------------------------------------------

 The number of warnings for this run is : 1
 
 -------------------------------------------------------------------------------
  **** **** ******  **  PROGRAM ENDED AT                 2020-12-17 16:46:22.496
 ***** ** ***  *** **   PROGRAM RAN ON                                 archlinux
 **    ****   ******    PROGRAM RAN BY                                   richard
 ***** **    ** ** **   PROGRAM PROCESS ID                                 23177
  **** **  *******  **  PROGRAM STOPPED IN /home/richard/work/Random_Stuff/temp/
                                           temp2/sept_2019_temp_work/cp2k_test_m
                                           d/temp/npt/print_cell
"""
