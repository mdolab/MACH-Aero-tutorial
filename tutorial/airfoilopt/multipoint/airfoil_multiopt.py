# ======================================================================
#         Import modules
# ======================================================================
#rst Imports (beg)
import numpy
from mpi4py import MPI
from baseclasses import *
from adflow import ADFLOW
from pygeo import *
from pyspline import *
from repostate import *
from pyoptsparse import Optimization, OPT
from pywarp import *
from multipoint import *
#rst Imports (end)
# ======================================================================
#         Specify parameters for optimization
# ======================================================================
#rst parameters (beg)
# specify flight conditions and constraints
mach = [0.4, 0.5]
alt = [10000, 10000]
alpha = [1, 1]
mycl = [0.5, 0.5]
# number of points in multipoint optimization
nFlowCases = len(mach)
# assign number of processors
nGroup = 1
nProcPerGroup = 4
#rst parameters (end)
# ======================================================================
#         Create multipoint communication object
# ======================================================================
#rst multipoint (beg)
MP = multiPointSparse(MPI.COMM_WORLD)
MP.addProcessorSet('cruise', nMembers=nGroup, memberSizes=nProcPerGroup)
comm, setComm, setFlags, groupFlags, ptID = MP.createCommunicators()
#rst multipoint (end)
# ======================================================================
#         ADflow Set-up
# ======================================================================
#rst adflow (beg)
aeroOptions = {
    # Common Parameters
    'gridFile':'n0012.cgns',
    'outputDirectory':'output',

    # Physics Parameters
    'equationType':'RANS',
    #'smoother':'dadi',

    'smoother':'runge kutta',
    'rkreset':True,
    'nrkreset':200,
    'CFL':0.8,
    'CFLCoarse':0.4,
    'MGCycle':'sg',
    'MGStartLevel':-1,
    'nCyclesCoarse':2500,
    'nCycles':20000,
    'monitorvariables':['resrho','cl','cd','cmz','yplus'],
    'useNKSolver':True,
    'useanksolver' : True,
    'nsubiterturb' : 10,
    'liftIndex':2,
    # Convergence Parameters
    'L2Convergence':1e-15,
    'L2ConvergenceCoarse':1e-4,
    # Adjoint Parameters
    'adjointSolver':'gmres', #gmres,tfqmr,rechardson,bcgs,ibcgs
    'adjointL2Convergence':1e-12,
    'ADPC':True,
    #'ADPC':False, #hxl
    'adjointMaxIter': 5000,
    'adjointSubspaceSize':400,
    'ILUFill':3,
    #'ILUFill':2, #hxl
    'ASMOverlap':3,
    'outerPreconIts':3,
    #'innerPreconIts':2, #hxl
    'NKSubSpaceSize':400,
    'NKASMOverlap':4,
    'NKPCILUFill':4,
    'NKJacobianLag':5,
    'nkswitchtol':1e-6, #2e-4,
    'nkouterpreconits':3,
    'NKInnerPreConIts':3,
    'writeSurfaceSolution':False,
    'writeVolumeSolution':False,
    'frozenTurbulence':False,
    'restartADjoint':False,
    }
# Create solver
CFDSolver = ADFLOW(options=aeroOptions, comm=comm)
CFDSolver.addLiftDistribution(200, 'z')
#rst adflow (end)
# ======================================================================
#         Set up flow conditions with AeroProblem
# ======================================================================
#rst aeroproblem (beg)
aeroProblems = []
for i in range(nFlowCases):
    ap = AeroProblem(name='fc%d'%i, alpha=alpha[i], mach=mach[i], altitude=alt[i],
                 areaRef=1.0, chordRef=1.0, evalFuncs=['cl','cd'])
    # Add angle of attack variable
    ap.addDV('alpha', value=alpha[i], lower=0, upper=10.0, scale=1.0)
    aeroProblems.append(ap)
#rst aeroproblem (end)
# ======================================================================
#         Geometric Design Variable Set-up
# ======================================================================
#rst dvgeo (beg)
# Create DVGeometry object
FFDFile = 'ffd.xyz'

DVGeo = DVGeometry(FFDFile)#DVGeo = DVGeometry_FFD(FFDFile)
DVGeo.addGeoDVLocal('shape', lower=-0.05, upper=0.05, axis='y', scale=1.0)

span = 1.0
pos = numpy.array([0.5])*span
CFDSolver.addSlices('z',pos,sliceType='absolute')

# Add DVGeo object to CFD solver
CFDSolver.setDVGeo(DVGeo)
#rst dvgeo (end)
# ======================================================================
#         DVConstraint Setup
# ======================================================================
#rst dvcon (beg)
DVCon = DVConstraints()#DVCon = DVConstraints_FFD_data()
DVCon.setDVGeo(DVGeo)

# Only ADflow has the getTriangulatedSurface Function
DVCon.setSurface(CFDSolver.getTriangulatedMeshSurface())

# Le/Te constraints
lIndex = DVGeo.getLocalIndex(0)
indSetA = []; indSetB = [];
#print('lIndex.shape[0]',lIndex.shape[0])
for k in range(0,1):
    indSetA.append(lIndex[0, 0, k]) #all DV for upper and lower should be same but different sign
    indSetB.append(lIndex[0, 1, k])
for k in range(0,1):
    indSetA.append(lIndex[-1, 0, k])
    indSetB.append(lIndex[-1, 1, k])
#print(indSetA)
DVCon.addLeTeConstraints(0, indSetA=indSetA, indSetB=indSetB)

# DV should be same along spanwise
lIndex = DVGeo.getLocalIndex(0)
indSetA = []; indSetB = [];
for i in range(lIndex.shape[0]):
    indSetA.append(lIndex[i, 0, 0])
    indSetB.append(lIndex[i, 0, 1])
for i in range(lIndex.shape[0]):
    indSetA.append(lIndex[i, 1, 0])
    indSetB.append(lIndex[i, 1, 1])
DVCon.addLinearConstraintsShape(indSetA, indSetB,factorA=1.0, factorB=-1.0,lower=0, upper=0)


le=0.0001
leList = [[le    , 0, le], [le    , 0, 1.0-le]]
teList = [[1.0-le, 0, le], [1.0-le, 0, 1.0-le]]


DVCon.addVolumeConstraint(leList, teList, 2,100, lower=0.064837137176294343, upper=0.07,scaled=False)
DVCon.addThicknessConstraints2D(leList, teList, 2,100, lower=0.1, upper=3.0)#lower=0.01, upper=3.0)


if comm.rank == 0:
    fileName = 'output/constraints.dat'
    DVCon.writeTecplot(fileName)

#rst dvcon (end)
# ======================================================================
#         Mesh Warping Set-up
# ======================================================================
#rst warp (beg)
meshOptions = {'gridFile':'n0012.cgns', 'warpType':'algebraic',}

mesh = MBMesh(options=meshOptions, comm=comm)
CFDSolver.setMesh(mesh)
#rst warp (end)
# ======================================================================
#         Functions:
# ======================================================================
#rst funcs (beg)
def cruiseFuncs(x):
    if MPI.COMM_WORLD.rank == 0:
        print x
    # Set design vars
    DVGeo.setDesignVars(x)
    # Evaluate functions
    funcs = {}
    DVCon.evalFunctions(funcs)

    for i in range(nFlowCases):
        if i%nGroup == ptID:
            aeroProblems[i].setDesignVars(x)
            CFDSolver(aeroProblems[i])
            CFDSolver.evalFunctions(aeroProblems[i], funcs)
            CFDSolver.checkSolutionFailure(aeroProblems[i], funcs)
    if MPI.COMM_WORLD.rank == 0:
        print funcs
    return funcs

def cruiseFuncsSens(x, funcs):
    funcsSens = {}
    DVCon.evalFunctionsSens(funcsSens)
    for i in range(nFlowCases):
        if i%nGroup == ptID:
            CFDSolver.evalFunctionsSens(aeroProblems[i], funcsSens)
    if MPI.COMM_WORLD.rank == 0:
        print funcsSens
    return funcsSens

def objCon(funcs, printOK):
    # Assemble the objective and any additional constraints:
    funcs['obj'] = 0.0
    for i in range(nFlowCases):
        ap = aeroProblems[i]
        funcs['obj'] += funcs[ap['cd']] / nFlowCases
        funcs['cl_con_'+ap.name] = funcs[ap['cl']] - mycl[i]
    if printOK:
       print 'funcs in obj:', funcs
    return funcs
#rst funcs (end)
# ======================================================================
#         Optimization Problem Set-up
# ======================================================================
#rst optprob (beg)
# Create optimization problem
optProb = Optimization('opt', MP.obj, comm=MPI.COMM_WORLD)

# Add objective
optProb.addObj('obj', scale=1e4)

# Add variables from the AeroProblem
ap.addVariablesPyOpt(optProb)

# Add DVGeo variables
DVGeo.addVariablesPyOpt(optProb)

# Add constraints
DVCon.addConstraintsPyOpt(optProb)
for i in range(nFlowCases):
    ap = aeroProblems[i]
    optProb.addCon('cl_con_'+ap.name, lower=0.0, upper=0.0, scale=1.0)

# The MP object needs the 'obj' and 'sens' function for each proc set,
# the optimization problem and what the objcon function is:
MP.setProcSetObjFunc('cruise', cruiseFuncs)
MP.setProcSetSensFunc('cruise', cruiseFuncsSens)
MP.setObjCon(objCon)
MP.setOptProb(optProb)
optProb.printSparsity()
#rst optprob (end)
#rst optimizer
# Set up optimizer
optOptions = {
    'Major iterations limit':200,
    'Major step limit':2.0,
    'Major feasibility tolerance':1.0e-6,
    'Major optimality tolerance':1.0e-6,
}
opt = OPT('snopt', options=optOptions)

# Run Optimization
sol = opt(optProb, MP.sens, storeHistory='opt.hst')
if MPI.COMM_WORLD.rank == 0:
   print sol
