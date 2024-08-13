import Monte
import poincare
import gbeespy as gbees

DIM = 4
km = Monte.UnitDbl(1, "km")
sec = Monte.UnitDbl(1, "sec")
begin_time = Monte.Epoch('1-jan-2000 00:00:00et')
boa = Monte.BoaLoad()
Monte.DefaultGm.addAll(boa)
boa = Monte.BoaLoad( ["./ephem/jup310.boa", "./ephem/de433.boa"] )
jupiter_europa_system = poincare.makeCR3BP( boa, 'jupiter', 'europa', epoch = begin_time)

# This function defines the dynamics model - required
def PCR3BP(x, dx, coef):
    x_nd = poincare.NONDIMENSIONAL_ST(x[0], x[1], 0, x[2], x[3], 0)
    x_d = poincare.makeDimensional(x_nd, jupiter_europa_system, 'LPO', epoch = begin_time)
    poincare.propagate(x_d, 0.01 * sec)
    query_LPO = Monte.TrajQuery(jupiter_europa_system.boa, 'LPO', 'Jupiter Barycenter', 'EME2000')
    x_d = query_LPO.pva(begin_time, jupiter_europa_system.rotatingFrame)
    v1 = x_d.vel()[0] * km/sec * jupiter_europa_system.Tstar/jupiter_europa_system.Lstar
    v2 = x_d.vel()[1] * km/sec * jupiter_europa_system.Tstar/jupiter_europa_system.Lstar
    v3 = x_d.acc()[0] * km/sec**2 * jupiter_europa_system.Tstar**2/jupiter_europa_system.Lstar
    v4 = x_d.acc()[1] * km/sec**2 * jupiter_europa_system.Tstar**2/jupiter_europa_system.Lstar
    return [v1, v2, v3, v4]

# This function defines the initial grid boundaries - optional
def PCR3BP_J(x, coef):
    r1 = ((x[0]+coef[0])**2+((x[1])**2))**0.5
    r2 = ((x[0]-1+coef[0])**2+(x[1])**2)**0.5
    J = (x[0])**2.0 + (x[1])**2.0 + (2*(1-coef[0])/r1) + (2*coef[0]/r2) + coef[0]*(1 - coef[0]) - ((x[2])**2.0 + (x[3])**2.0)
    return J

#==================================== Read in initial discrete measurement ==================================#
print("\nReading in initial discrete measurement...\n\n")

P_DIR = "./results"           # Saved PDFs path
M_DIR = "."                   # Measurement path
M_FILE = "/measurement0.txt"; # Measurement file
M = gbees.Meas_create(DIM, M_DIR, M_FILE)
#============================================================================================================#

#=========================================== Read in user inputs ============================================#
print("Reading in user inputs...\n\n")

dx = [None] * DIM                             # Grid width, default is half of the std. dev. from the initial measurement 
for i in range(DIM):
    dx[i] = (M.cov[i][i]**(0.5))/2
G = gbees.Grid_create(DIM, 5E-8, M.mean, dx); # Inputs: (dimension, probability threshold, center, grid width)    
 
coef = [2.528017528540000E-5]                 # PCR3BP trajectory attributes (mu)
T = gbees.Traj_create(len(coef), coef);       # Inputs: (# of coefficients, coefficients)

OUTPUT_FREQ = 20;                             # Number of steps per output to terminal
DEL_STEP = 20;                                # Number of steps per deletion procedure
NUM_DIST = 17;                                # Number of distributions recorded per measurement
NUM_MEAS = 1;                                 # Number of measurements
OUTPUT = False;                               # Write info to terminal
RECORD = True;                                # Write PDFs to .txt file
MEASURE = True;                               # Take discrete measurement updates
BOUNDS = False;                                # Add inadmissible regions to grid
#============================================================================================================#

#================================================== GBEES ===================================================#
gbees.run_gbees(PCR3BP, PCR3BP_J, G, M, T, P_DIR, M_DIR, NUM_DIST, NUM_MEAS, DEL_STEP, OUTPUT_FREQ, DIM, OUTPUT, RECORD, MEASURE, BOUNDS)