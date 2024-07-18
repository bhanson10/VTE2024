import ctypes as ct
import time

DIM = 4 # Dimensionality

class TreeNode(ct.Structure):
    pass

class Cell(ct.Structure):
    _fields_ = [
        ("prob", ct.c_double),
        ("v", ct.c_double * DIM),
        ("u", ct.c_double * DIM),
        ("w", ct.c_double * DIM),
        ("ctu", ct.c_double * DIM),
        ("state", ct.c_int * DIM),
        ("i_nodes", ct.POINTER(TreeNode) * DIM),
        ("k_nodes", ct.POINTER(TreeNode) * DIM),
        ("dcu", ct.c_double),
        ("cfl_dt", ct.c_double),
        ("new_f", ct.c_int),
        ("ik_f", ct.c_int),
        ("del_f", ct.c_int)
    ]

TreeNode._fields_ = [
        ("key", ct.c_ulong),
        ("cell", Cell),
        ("left", ct.POINTER(TreeNode)),
        ("right", ct.POINTER(TreeNode))
    ]

class Grid(ct.Structure):
    _fields_ = [
        ("thresh", ct.c_double),
        ("dt", ct.c_double),
        ("center", ct.c_double * DIM),
        ("dx", ct.c_double * DIM)
    ]

class Traj(ct.Structure):
    _fields_ = [
        ("mu", ct.c_double)
    ]

class BST(ct.Structure):
    _fields_ = [
        ("dead", ct.POINTER(TreeNode)),
        ("root", ct.POINTER(TreeNode)),
        ("a_count", ct.c_int), 
        ("tot_count", ct.c_int), 
        ("max_key", ct.c_ulong), 
        ("cfl_min_dt", ct.c_float)
    ]

def initialize_vuw(r, G, T):
    if not r:
        return 
    
    initialize_vuw(r.contents.left, G, T)
    initialize_vuw(r.contents.right, G, T)
    
    if r.contents.cell.new_f == 0:
        x = [0, 0, 0, 0]
        for i in range(DIM):
            x[i] = G.dx[i]*r.contents.cell.state[i] + G.center[i]

        r1 = ((x[0]+T.mu)**2+x[1]**2)**1.5
        r2 = ((x[0]-1+T.mu)**2+x[1]**2)**1.5

        v = [x[2], x[3], 2*x[3]+x[0]-(T.mu*(x[0]-1+T.mu)/r2)-((1-T.mu)*(x[0]+T.mu)/r1), -2*x[2]+x[1]-(T.mu*x[1]/r2)-((1-T.mu)*x[1]/r1)]

        for i in range(DIM):
            r.contents.cell.v[i] = v[i]
            r.contents.cell.u[i] = min(v[i], 0.0)
            r.contents.cell.w[i] = max(v[i], 0.0)
        r.contents.cell.new_f = 1

        sum = 0
        for q in range(DIM):
            sum += abs(r.contents.cell.v[q])/G.dx[q]
        r.contents.cell.cfl_dt = 1/sum

#==================================== ctypes initialization ================================================#
gbees = ct.CDLL("gbees.so")
gbees.initialize_gbees.argtypes = [ct.c_char_p, ct.POINTER(Grid), ct.POINTER(ct.c_double), Traj]
gbees.initialize_gbees.restype = BST
gbees.record_data.argtypes = [ct.POINTER(TreeNode), ct.c_char_p, Grid, ct.c_double]
gbees.record_data.restype = None
gbees.grow_tree.argtypes = [ct.POINTER(BST), Grid, Traj]
gbees.grow_tree.restype = None
gbees.update_prob.argtypes = [ct.POINTER(BST), ct.POINTER(Grid), ct.c_double]
gbees.update_prob.restype = None
gbees.prune_tree.argtypes = [ct.POINTER(BST), ct.POINTER(Grid)]
gbees.prune_tree.restype = None
gbees.get_tree_info.argtypes = [ct.POINTER(TreeNode), ct.POINTER(BST), ct.POINTER(Grid)]
gbees.get_tree_info.restype = None
gbees.measurement_update.argtypes = [ct.c_char_p, ct.POINTER(BST), ct.POINTER(ct.c_double)]
gbees.measurement_update.restype = None
#===========================================================================================================#

#======================================= Read in user inputs ===============================================#
FILE_PATH = "./Data"
FILE_NAME = FILE_PATH + "/Measurements/measurement0.txt"
G = Grid(thresh = 1E-7)
T = Traj(mu = 2.528017528540000E-5)

NM          = 1     # Number of discrete measurements
OUTPUT      = True  # Write info to terminal
RECORD      = True  # Write PDFs to .txt file
MEASURE     = True  # Take discrete measurement updates
OUTPUT_FREQ = 20    # Number of steps per output to terminal
DEL_STEP    = 10    # Number of steps per deletion procedure
NUM_DIST    = 17    # Number of distributions recorded per measurement
#===========================================================================================================#

#=============================================== GBEES =====================================================#
measure_time = ct.c_double(0); 
P = gbees.initialize_gbees(FILE_NAME.encode('utf-8'), ct.pointer(G), ct.pointer(measure_time), T)
initialize_vuw(P.root, G, T)
record_time = (measure_time.value)/(NUM_DIST-1.0)

print("Entering time marching...\n")

start = time.time()
tt = 0.0
for nm in range(NM):
    gbees.get_tree_info(P.root, ct.pointer(P),  ct.pointer(G))
    finish = time.time()
    print("Timestep: " + str(nm) + "-0, Program time: " + str(finish - start) + " s, Sim. time: " + str(tt) + " TU, Active/Total Cells: " + str(P.a_count) + "/" + str(P.tot_count) + ", Max key %: " + str(P.max_key/((2**64)-1)*100) + "%")
    if RECORD:
        FILE_NAME = FILE_PATH + "/PDFs/P" + str(nm) + "/pdf_0.txt"
        gbees.record_data(P.root, FILE_NAME.encode('utf-8'), G, ct.c_double(0.0))

    mt = 0.0
    record_count = 1
    step_count = 1
    while(mt < measure_time.value):

        rt = 0.0
        while(rt < record_time):
            
            gbees.grow_tree(ct.pointer(P), G, T)
            initialize_vuw(P.root, G, T)
            gbees.update_prob(ct.pointer(P), ct.pointer(G), ct.c_double(record_time - rt))

            if (step_count % DEL_STEP == 0):
                gbees.prune_tree(ct.pointer(P), ct.pointer(G))

            if OUTPUT and step_count % OUTPUT_FREQ == 0:
                gbees.get_tree_info(P.root, ct.pointer(P),  ct.pointer(G))
                finish = time.time()
                print("Timestep: " + str(nm) + "-" + str(step_count) + ", Program time: " + str(finish - start) + " s, Sim. time: " + str(tt + mt + rt) + " TU, Active/Total Cells: " + str(P.a_count) + "/" + str(P.tot_count) + ", Max key %: " + str(P.max_key/((2**64)-1)*100) + "%")
            step_count += 1

            rt += G.dt

        if OUTPUT and step_count % OUTPUT_FREQ != 0:
            gbees.get_tree_info(P.root, ct.pointer(P),  ct.pointer(G))
            finish = time.time()
            print("Timestep: " + str(nm) + "-" + str(step_count) + ", Program time: " + str(finish - start) + " s, Sim. time: " + str(tt + mt + rt) + " TU, Active/Total Cells: " + str(P.a_count) + "/" + str(P.tot_count) + ", Max key %: " + str(P.max_key/((2**64)-1)*100) + "%")

        if RECORD:
            print("\nRECORDING PDF AT: " + str(tt + mt + rt) + " TU...\n")
            FILE_NAME = FILE_PATH + "/PDFs/P" + str(nm) + "/pdf_" + str(record_count) + ".txt"
            gbees.record_data(P.root, FILE_NAME.encode('utf-8'), G, ct.c_double(tt + mt + rt))
            record_count += 1
    
        mt += rt

    tt += mt
    if MEASURE and nm < NM-1:
        print("\nPERFORMING BAYESIAN UPDATE AT: " + str(tt) + " TU...\n")
        nm += 1

        FILE_NAME = FILE_PATH + "/Measurements/measurement" + str(nm) + ".txt"

        gbees.measurement_update(FILE_NAME.encode('utf-8'), ct.pointer(P), G, ct.pointer(measure_time));