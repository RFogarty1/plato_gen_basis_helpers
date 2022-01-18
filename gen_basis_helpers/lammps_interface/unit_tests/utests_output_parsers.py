

import unittest
import unittest.mock as mock


import gen_basis_helpers.analyse_md.thermo_data as thermoDataObjs
import gen_basis_helpers.lammps_interface.lammps_parsers as tCode

class TestParseLogFileData(unittest.TestCase):

	def setUp(self):
		self.mockFilePath = "fake_path"
		self.createTestObjs()

	def createTestObjs(self):
		self.thermoSectAsListA = _loadThermoSectStrA().split("\n")
		self.fullFileAsListA = _loadSmallFullFileStrA().split("\n")

	@mock.patch("gen_basis_helpers.lammps_interface.lammps_parsers._getFileAsListFromInpPath")
	def testParseThermoSection(self, mockedGetFileAsList):
		mockedGetFileAsList.side_effect = lambda *args: self.fullFileAsListA
		expObj = self._loadExpThermoObjA()
		actObj = tCode.parseLammpsLogFile(self.mockFilePath)["thermo_data"]
		self.assertEqual(expObj,actObj)

	def _loadExpThermoObjA(self):
		expSteps = [0, 50, 100, 150, 200]
		expTemp = [300, 5142.1978, 1916.2449, 981.71584, 698.04681]
		expETot = [12712.121, 4186.2035, 2012.0162, 1631.1898, 993.1143]
		expPressure = [2186416, 121501.67, 49433.114, 54753.846, 40100.158]
		expTime = [2*x for x in expSteps]
		currKwargs = {"step":expSteps, "time":expTime, "temp":expTemp,
		              "eTotal":expETot, "pressure":expPressure}
		return thermoDataObjs.ThermoDataStandard.fromStdKwargs(**currKwargs)

	def testParseThermoSectionWithPotAndKinEnergy(self):
		thermoSection = _loadThermoSectStrB().split("\n")
		inpLineIdx = 0
		expObj = self._loadExpThermoObjWithPotAndKinEnergy()
		unused, actObj = tCode._parseThermoSection(thermoSection, inpLineIdx)
		self.assertEqual(expObj, actObj)

	def _loadExpThermoObjWithPotAndKinEnergy(self):
		expSteps = [0, 50, 100, 150, 200]
		expTemp = [300, 299.96484, 384.28311, 413.1955, 452.95121]
		expETot = [-438.07601, -438.26344, -438.33272, -438.49712, -438.99766]
		expKinE = [334.44689, 334.40769, 428.40764, 460.63983, 504.9604]
		expPotE = [-772.5229, -772.67113, -866.74036, -899.13695, -943.95807]
		currKwargs = {"step":expSteps, "temp":expTemp, "eTotal":expETot,
		              "eKinetic":expKinE, "ePot":expPotE}
		return thermoDataObjs.ThermoDataStandard.fromStdKwargs(**currKwargs)


def _loadThermoSectStrA():
	outStr="""Step Temp E_pair E_mol TotEng Press 
       0          300     12557.08    59.357396    12712.121      2186416 
      50    5142.1978    1580.3068    965.80983    4186.2035    121501.67 
     100    1916.2449    1102.9612     297.8751    2012.0162    49433.114 
     150    981.71584    1128.1042    189.97061    1631.1898    54753.846 
     200    698.04681     674.3612    96.113385     993.1143    40100.158 
Loop time of 23.25 on 1 procs for 50000 steps with 108 atoms"""
	return outStr

def _loadThermoSectStrB():
	outStr ="""Step Temp PotEng KinEng TotEng 
       0          300    -772.5229    334.44689   -438.07601 
      50    299.96484   -772.67113    334.40769   -438.26344 
     100    384.28311   -866.74036    428.40764   -438.33272 
     150     413.1955   -899.13695    460.63983   -438.49712 
     200    452.95121   -943.95807     504.9604   -438.99766 
Loop time of 28.906 on 1 procs for 10000 steps with 375 atoms"""
	return outStr




def _loadSmallFullFileStrA():
	outStr = """LAMMPS (29 Oct 2020)
OMP_NUM_THREADS environment is not set. Defaulting to 1 thread. (src/mylammps/src/comm.cpp:94)
  using 1 OpenMP thread(s) per MPI task
units real
atom_style full
boundary p p f
read_data inp_file_1.data
Reading data file ...
  triclinic box = (0.0000000 0.0000000 0.0000000) to (9.6300000 8.3398200 11.477460) with tilt (-4.8101800 0.0000000 0.0000000)
  1 by 1 by 1 MPI processor grid
  reading atoms ...
  108 atoms
  scanning bonds ...
  2 = max bonds/atom
  scanning angles ...
  1 = max angles/atom
  reading bonds ...
  72 bonds
  reading angles ...
  36 angles
Finding 1-2 1-3 1-4 neighbors ...
  special bond factors lj:    0.0      0.0      0.0     
  special bond factors coul:  0.0      0.0      0.0     
     2 = max # of 1-2 neighbors
     1 = max # of 1-3 neighbors
     1 = max # of 1-4 neighbors
     2 = max # of special neighbors
  special bonds CPU = 0.000 seconds
  read_data CPU = 0.018 seconds
pair_style lj/cut/coul/cut 10 10
bond_style harmonic
angle_style harmonic
pair_coeff 1 1 0.1521 3.1507
pair_coeff 1 2 0.0000 0.0000
pair_coeff 2 2 0.0000 0.0000
bond_coeff 1 450.0000 0.9572
angle_coeff 1 55.0000 104.5200
velocity all create 300.0 200 dist uniform
timestep 2.0
thermo_style one
thermo 50
fix 1 all nvt temp 300.0 300.0 200.0
fix 2 all wall/reflect zlo EDGE
fix 3 all wall/reflect zhi EDGE
dump myDump all atom 50 dump.lammpstrj
dump_modify myDump scale no
run 50000
Neighbor list info ...
  update every 1 steps, delay 10 steps, check yes
  max neighbors/atom: 2000, page size: 100000
  master list distance cutoff = 12
  ghost atom cutoff = 12
  binsize = 6, bins = 3 2 2
  1 neighbor lists, perpetual/occasional/extra = 1 0 0
  (1) pair lj/cut/coul/cut, perpetual
      attributes: half, newton on
      pair build: half/bin/newton/tri
      stencil: half/bin/3d/newton/tri
      bin: standard
Per MPI rank memory allocation (min/avg/max) = 7.653 | 7.653 | 7.653 Mbytes
Step Temp E_pair E_mol TotEng Press 
       0          300     12557.08    59.357396    12712.121      2186416 
      50    5142.1978    1580.3068    965.80983    4186.2035    121501.67 
     100    1916.2449    1102.9612     297.8751    2012.0162    49433.114 
     150    981.71584    1128.1042    189.97061    1631.1898    54753.846 
     200    698.04681     674.3612    96.113385     993.1143    40100.158 
Loop time of 23.25 on 1 procs for 50000 steps with 108 atoms

Performance: 371.612 ns/day, 0.065 hours/ns, 2150.533 timesteps/s
99.0% CPU use with 1 MPI tasks x 1 OpenMP threads

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Pair    | 19.386     | 19.386     | 19.386     |   0.0 | 83.38
Bond    | 0.17117    | 0.17117    | 0.17117    |   0.0 |  0.74
Neigh   | 2.8825     | 2.8825     | 2.8825     |   0.0 | 12.40
Comm    | 0.35889    | 0.35889    | 0.35889    |   0.0 |  1.54
Output  | 0.25069    | 0.25069    | 0.25069    |   0.0 |  1.08
Modify  | 0.12637    | 0.12637    | 0.12637    |   0.0 |  0.54
Other   |            | 0.07417    |            |       |  0.32

Nlocal:        108.000 ave         108 max         108 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Nghost:        1530.00 ave        1530 max        1530 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Neighs:        28629.0 ave       28629 max       28629 min
Histogram: 1 0 0 0 0 0 0 0 0 0

Total # of neighbors = 28629
Ave neighs/atom = 265.08333
Ave special neighs/atom = 2.0000000
Neighbor list builds = 2987
Dangerous builds = 113
Total wall time: 0:00:23"""
	return outStr
