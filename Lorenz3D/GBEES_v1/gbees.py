import ctypes as ct
import time
import cProfile
import pstats

DIM = 3 # Dimensionality

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
        ("sigma", ct.c_double),
        ("b", ct.c_double),
        ("r", ct.c_double)
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
        x = [0, 0, 0]
        for i in range(DIM):
            x[i] = G.dx[i]*r.contents.cell.state[i] + G.center[i]
        v1 = T.sigma*(x[1] - (x[0] + (G.dx[0]/2.0)))
        v2 = -(x[1] + (G.dx[1]/2.0)) - x[0]*x[2]
        v3 = -T.b*(x[2] + (G.dx[2]/2.0)) + x[0]*x[1] - T.b * T.r 
        r.contents.cell.v = (ct.c_double * DIM)(*[v1, v2, v3])
        r.contents.cell.u = (ct.c_double * DIM)(*[min(v1, 0.0), min(v2, 0.0), min(v3, 0.0)])
        r.contents.cell.w = (ct.c_double * DIM)(*[max(v1, 0.0), max(v2, 0.0), max(v3, 0.0)])
        r.contents.cell.new_f = 1

        sum = 0
        for q in range(DIM):
            sum += abs(r.contents.cell.v[q])/G.dx[q]
        r.contents.cell.cfl_dt = 1/sum

if __name__ == "__main__":
    profiler = cProfile.Profile()
    profiler.enable()
    #==================================== ctypes initialization ================================================#
    gbees = ct.CDLL("gbees.so")
    gbees.initialize_gbees.argtypes = [ct.c_char_p, ct.POINTER(Grid), ct.POINTER(ct.c_double)]
    gbees.initialize_gbees.restype = BST
    gbees.record_data.argtypes = [ct.POINTER(TreeNode), ct.c_char_p, Grid, ct.c_double]
    gbees.record_data.restype = None
    gbees.grow_tree.argtypes = [ct.POINTER(BST), Grid]
    gbees.grow_tree.restype = None
    gbees.initialize_vuw.argtypes = [ct.POINTER(TreeNode), Grid, Traj]
    gbees.initialize_vuw.restype = None
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
    FILE_NAME = FILE_PATH + "/measurements0.txt"
    G = Grid(thresh = 2E-5)
    T = Traj(sigma = 4, b = 1, r = 48)

    NM          = 1     # Number of discrete measurements
    OUTPUT      = True  # Write info to terminal
    RECORD      = True  # Write PDFs to .txt file
    MEASURE     = True  # Take discrete measurement updates
    OUTPUT_FREQ = 10    # Number of steps per output to terminal
    DEL_STEP    = 25    # Number of steps per deletion procedure
    NUM_DIST    = 6     # Number of distributions recorded per measurement
    #===========================================================================================================#

    #=============================================== GBEES =====================================================#
    measure_time = ct.c_double(0); 
    P = gbees.initialize_gbees(FILE_NAME.encode('utf-8'), ct.POINTER(Grid)(G), ct.POINTER(ct.c_double)(measure_time))
    gbees.initialize_vuw(P.root, G, T)
    # initialize_vuw(P.root, G, T)
    record_time = (measure_time.value)/(NUM_DIST-1.0)

    print("Entering time marching...\n")

    start = time.time()
    tt = 0.0
    for nm in range(NM):
        gbees.get_tree_info(P.root, ct.POINTER(BST)(P),  ct.POINTER(Grid)(G))
        finish = time.time()
        print("Timestep: " + str(nm) + "-0, Program time: " + str(finish - start) + " s, Sim. time: " + str(tt) + " TU, Active/Total Cells: " + str(P.a_count) + "/" + str(P.tot_count) + ", Max key %: " + str(P.max_key/((2**64)-1)*100) + "%")
        if RECORD:
            FILE_NAME = FILE_PATH + "/M" + str(nm) + "/pdf_0.txt"
            gbees.record_data(P.root, FILE_NAME.encode('utf-8'), G, ct.c_double(0.0))

        mt = 0.0
        record_count = 1
        step_count = 1
        while(mt < measure_time.value):

            rt = 0.0
            while(rt < record_time):
                
                gbees.grow_tree(ct.POINTER(BST)(P), G)
                gbees.initialize_vuw(P.root, G, T)
                # initialize_vuw(P.root, G, T)
                gbees.update_prob(ct.POINTER(BST)(P), ct.POINTER(Grid)(G), ct.c_double(record_time - rt))

                if (step_count % DEL_STEP == 0):
                    gbees.prune_tree(ct.POINTER(BST)(P), ct.POINTER(Grid)(G))

                if OUTPUT and step_count % OUTPUT_FREQ == 0:
                    gbees.get_tree_info(P.root, ct.POINTER(BST)(P),  ct.POINTER(Grid)(G))
                    finish = time.time()
                    print("Timestep: " + str(nm) + "-" + str(step_count) + ", Program time: " + str(finish - start) + " s, Sim. time: " + str(tt + mt + rt) + " TU, Active/Total Cells: " + str(P.a_count) + "/" + str(P.tot_count) + ", Max key %: " + str(P.max_key/((2**64)-1)*100) + "%")
                step_count += 1

                rt += G.dt

            if OUTPUT and step_count % OUTPUT_FREQ != 0:
                gbees.get_tree_info(P.root, ct.POINTER(BST)(P),  ct.POINTER(Grid)(G))
                finish = time.time()
                print("Timestep: " + str(nm) + "-" + str(step_count) + ", Program time: " + str(finish - start) + " s, Sim. time: " + str(tt + mt + rt) + " TU, Active/Total Cells: " + str(P.a_count) + "/" + str(P.tot_count) + ", Max key %: " + str(P.max_key/((2**64)-1)*100) + "%")

            if RECORD:
                print("\nRECORDING PDF AT: " + str(tt + mt + rt) + " TU...\n")
                FILE_NAME = FILE_PATH + "/M" + str(nm) + "/pdf_" + str(record_count) + ".txt"
                gbees.record_data(P.root, FILE_NAME.encode('utf-8'), G, ct.c_double(tt + mt + rt))
                record_count += 1
        
            mt += rt

        tt += mt
        if MEASURE and nm < NM-1:
            print("\nPERFORMING BAYESIAN UPDATE AT: " + str(tt) + " TU...\n")
            nm += 1

            FILE_NAME = FILE_PATH + "/measurements" + str(nm) + ".txt"

            gbees.measurement_update(FILE_NAME.encode('utf-8'), ct.POINTER(BST)(P), G, ct.POINTER(ct.c_double)(measure_time));

    profiler.disable()
    results = pstats.Stats(profiler)
    results.dump_stats("./Data/results_wrapper.prof")