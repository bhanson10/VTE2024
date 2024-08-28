import Monte # type: ignore
import poincare # type: ignore
import gbeespy_m as gbees
import math 

DIM_f = 4 # State dimension
DIM_h = 3 # Measurement dimension

km = Monte.UnitDbl(1, "km")
sec = Monte.UnitDbl(1, "sec")
begin_time = Monte.Epoch('1-jan-2000 00:00:00et')
boa = Monte.BoaLoad()
Monte.DefaultGm.addAll(boa)
boa = Monte.BoaLoad( ["./ephem/sat375l.boa", "./ephem/de430.boa"] )
sa_en_sys = poincare.makeCR3BP( boa, 'saturn', 'enceladus', epoch = begin_time)

# This function defines the dynamics model - required
def PCR3BP(x, dx, coef):
    x_nd = poincare.NONDIMENSIONAL_ST(x[0], x[1], 0, x[2], x[3], 0)
    x_d = poincare.makeDimensional(x_nd, sa_en_sys, 'DPO', epoch = begin_time)
    poincare.propagate(x_d, 0.01 * sec)
    query_DPO = Monte.TrajQuery(sa_en_sys.boa, 'DPO', 'Saturn Barycenter', 'EME2000')
    x_d = query_DPO.pva(begin_time, sa_en_sys.rotatingFrame)
    v3 = x_d.acc()[0] * km/sec**2 * sa_en_sys.Tstar**2/sa_en_sys.Lstar
    v4 = x_d.acc()[1] * km/sec**2 * sa_en_sys.Tstar**2/sa_en_sys.Lstar
    return [x[2], x[3], v3, v4]

# This function defines the measurement model - required
def rtrr(x, dx, coef):
    v1 = ((x[0] - (1- coef[0]))**2 + (x[1])**2)**0.5 
    v2 = math.atan2(x[1],  x[0] - (1 - coef[0]))
    v3 = ((x[0] - (1 - coef[0]))*x[2] + x[1]*x[3])/v1
    return [v1, v2, v3]

# This function defines the initial grid boundaries - optional
def PCR3BP_J(x, coef):
    r1 = ((x[0]+coef[0])**2+(x[1])**2)**0.5
    r2 = ((x[0]-1+coef[0])**2+(x[1])**2)**0.5
    J = (x[0])**2.0 + (x[1])**2.0 + (2*(1-coef[0])/r1) + (2*coef[0]/r2) + coef[0]*(1 - coef[0]) - ((x[2])**2.0 + (x[3])**2.0)
    return J

#==================================== Read in initial discrete measurement ==================================#
print("Reading in initial discrete measurement...\n")

P_DIR = "./results/monte"    # Saved PDFs path
M_DIR = "./measurements"     # Measurement path
M_FILE = "measurement0.txt"; # Measurement file
M = gbees.Meas_create(DIM_f, M_DIR, M_FILE)
#============================================================================================================#

#=========================================== Read in user inputs ============================================#
print("Reading in user inputs...\n")

dx = [None] * DIM_f                             # Grid width, default is half of the std. dev. from the initial measurement 
for i in range(DIM_f):
    dx[i] = (M.cov[i][i]**(0.5))/2
G = gbees.Grid_create(DIM_f, 5E-8, M.mean, dx); # Inputs: (dimension, probability threshold, center, grid width)    
 
coef = [1.901109735892602E-07]                # PCR3BP trajectory attributes (mu)
T = gbees.Traj_create(len(coef), coef);       # Inputs: (# of coefficients, coefficients)

NUM_DIST = 8;                                 # Number of distributions recorded per measurement
NUM_MEAS = 4;                                 # Number of measurements
DEL_STEP = 20;                                # Number of steps per deletion procedure
OUTPUT_FREQ = 20;                             # Number of steps per output to terminal
OUTPUT = False;                               # Write info to terminal
RECORD = True;                                # Write PDFs to .txt file
MEASURE = True;                               # Take discrete measurement updates
BOUNDS = True;                                # Add inadmissible regions to grid
#============================================================================================================#

#================================================== GBEES ===================================================#
gbees.run_gbees(PCR3BP, rtrr, PCR3BP_J, G, M, T, P_DIR, M_DIR, NUM_DIST, NUM_MEAS, DEL_STEP, OUTPUT_FREQ, DIM_h, OUTPUT, RECORD, MEASURE, BOUNDS)