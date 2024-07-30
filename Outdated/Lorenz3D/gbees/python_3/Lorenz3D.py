import gbeespy as gbees

def Lorenz3D(x, dx, T):
    v1 = T.coef[0]*(x[1] - (x[0] + (dx[0]/2.0)))
    v2 = -(x[1] + (dx[1]/2.0)) - x[0]*x[2]
    v3 = -T.coef[1]*(x[2] + (dx[2]/2.0)) + x[0]*x[1] - T.coef[1]*T.coef[2]
    return [v1, v2, v3]
 
FILE_PATH = "./Data"                    # Define directory where Measurements and PDFs are stored
G = gbees.Grid(dim = 3, thresh = 4E-5)  # Define grid probability threshold
T = gbees.Traj(coef = [4.0, 1.0, 48.0]) # Define trajectory coefficients

NM          = 1     # Number of discrete measurements
OUTPUT      = True  # Write info to terminal
RECORD      = True  # Write PDFs to .txt file
MEASURE     = True  # Take discrete measurement updates
OUTPUT_FREQ = 20    # Number of steps per output to terminal
DEL_STEP    = 25    # Number of steps per deletion procedure
NUM_DIST    = 6     # Number of distributions recorded per measurement

gbees.propagate_uncertainty(G, T, Lorenz3D, NM, FILE_PATH, NUM_DIST, RECORD, DEL_STEP, OUTPUT, OUTPUT_FREQ, MEASURE)