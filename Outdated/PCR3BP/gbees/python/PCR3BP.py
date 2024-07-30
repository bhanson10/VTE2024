import gbeespy as gbees

def PCR3BP(x, dx, T):
    r1 = ((x[0]+T.coef[0])**2+x[1]**2)**1.5
    r2 = ((x[0]-1+T.coef[0])**2+x[1]**2)**1.5

    v1 = x[2]
    v2 = x[3]
    v3 = 2*x[3]+x[0]-(T.coef[0]*(x[0]-1+T.coef[0])/r2)-((1-T.coef[0])*(x[0]+T.coef[0])/r1)
    v4 = -2*x[2]+x[1]-(T.coef[0]*x[1]/r2)-((1-T.coef[0])*x[1]/r1)

    return [v1, v2, v3, v4]
 
FILE_PATH = "./Data"                          # Define directory where Measurements and PDFs are stored
G = gbees.Grid(dim = 4, thresh = 1E-7)        # Define grid probability threshold
T = gbees.Traj(coef = [2.528017528540000E-5]) # Define trajectory coefficients

NM          = 1     # Number of discrete measurements
OUTPUT      = True  # Write info to terminal
RECORD      = True  # Write PDFs to .txt file
MEASURE     = True  # Take discrete measurement updates
OUTPUT_FREQ = 20    # Number of steps per output to terminal
DEL_STEP    = 10    # Number of steps per deletion procedure
NUM_DIST    = 17    # Number of distributions recorded per measurement

gbees.propagate_uncertainty(G, T, PCR3BP, NM, FILE_PATH, NUM_DIST, RECORD, DEL_STEP, OUTPUT, OUTPUT_FREQ, MEASURE)