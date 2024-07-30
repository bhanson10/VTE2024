import ctypes as ct
import time 

DIM = 4    # Dimensionality
TOL = 1E-8 # Tolerance 

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
        ("dim", ct.c_int),
        ("thresh", ct.c_double),
        ("dt", ct.c_double),
        ("center", ct.c_double * DIM),
        ("dx", ct.c_double * DIM)
    ]

class Traj(ct.Structure):
    _fields_ = [("coef", ct.POINTER(ct.c_double))]
    max_length = 50

    def __init__(self, coef):
        if len(coef) > self.max_length:
            raise ValueError(f"The list of coefficients cannot exceed {self.max_length} elements.")
        self._coef = (ct.c_double * len(coef))(*coef)
        self.coef = ct.cast(self._coef, ct.POINTER(ct.c_double))
        self.length = len(coef)

class BST(ct.Structure):
    _fields_ = [
        ("dead", ct.POINTER(TreeNode)),
        ("root", ct.POINTER(TreeNode)),
        ("a_count", ct.c_int),
        ("tot_count", ct.c_int),
        ("max_key", ct.c_ulong),
        ("cfl_min_dt", ct.c_float)
    ]

def initialize_vuw(r, G, T, adv):
    if not r:
        return

    initialize_vuw(r.contents.left, G, T, adv)
    initialize_vuw(r.contents.right, G, T, adv)

    if r.contents.cell.new_f == 0:
        x = [None] * G.dim
        for i in range(G.dim):
            x[i] = G.dx[i] * r.contents.cell.state[i] + G.center[i]

        v = adv(x, G.dx, T)

        for i in range(G.dim):
            r.contents.cell.v[i] = v[i]
            r.contents.cell.u[i] = min(v[i], 0.0)
            r.contents.cell.w[i] = max(v[i], 0.0)
        r.contents.cell.new_f = 1

        sum = 0
        for q in range(G.dim):
            sum += abs(r.contents.cell.v[q]) / G.dx[q]
        r.contents.cell.cfl_dt = 1 / sum

def propagate_uncertainty(G, T, advection, NM, FILE_PATH, NUM_DIST=6, RECORD=True, DEL_STEP=25, OUTPUT=True, OUTPUT_FREQ=20, MEASURE=True):
    lib = ct.CDLL("gbees.so")
    lib.initialize_gbees.argtypes = [ct.c_char_p, ct.POINTER(Grid), ct.POINTER(ct.c_double)]
    lib.initialize_gbees.restype = BST
    lib.record_data.argtypes = [ct.POINTER(TreeNode), ct.c_char_p, Grid, ct.c_double]
    lib.record_data.restype = None
    lib.grow_tree.argtypes = [ct.POINTER(BST), Grid]
    lib.grow_tree.restype = None
    lib.update_prob.argtypes = [ct.POINTER(BST), ct.POINTER(Grid), ct.c_double]
    lib.update_prob.restype = None
    lib.prune_tree.argtypes = [ct.POINTER(BST), ct.POINTER(Grid)]
    lib.prune_tree.restype = None
    lib.get_tree_info.argtypes = [ct.POINTER(TreeNode), ct.POINTER(BST), ct.POINTER(Grid)]
    lib.get_tree_info.restype = None
    lib.measurement_update.argtypes = [ct.c_char_p, ct.POINTER(BST), ct.POINTER(ct.c_double)]
    lib.measurement_update.restype = None

    measure_time = ct.c_double(0); 
    FILE_NAME = FILE_PATH + "/Measurements/measurement0.txt"
    P = lib.initialize_gbees(FILE_NAME.encode('utf-8'), ct.pointer(G), ct.pointer(measure_time))
    initialize_vuw(P.root, G, T, advection)
    record_time = (measure_time.value)/(NUM_DIST-1.0)
    
    print("Entering time marching...\n")

    start = time.time()
    tt = 0.0
    for nm in range(NM):
        lib.get_tree_info(P.root, ct.pointer(P),  ct.pointer(G))
        finish = time.time()
        if OUTPUT:
            print("Timestep: " + str(nm) + "-0, Program time: " + str(finish - start) + " s, Sim. time: " + str(tt) + " TU, Active/Total Cells: " + str(P.a_count) + "/" + str(P.tot_count) + ", Max key %: " + str(P.max_key/((2**64)-1)*100) + "%")
        if RECORD:
            FILE_NAME = FILE_PATH + "/PDFs/P" + str(nm) + "/pdf_0.txt"
            lib.record_data(P.root, FILE_NAME.encode('utf-8'), G, ct.c_double(0.0))

        mt = 0.0
        record_count = 1
        step_count = 1
        while(abs(mt - measure_time.value) > TOL):

            rt = 0.0
            while(rt < record_time):
                
                lib.grow_tree(ct.pointer(P), G)
                initialize_vuw(P.root, G, T, advection)
                lib.update_prob(ct.pointer(P), ct.pointer(G), ct.c_double(record_time - rt))

                if (step_count % DEL_STEP == 0):
                    lib.prune_tree(ct.pointer(P), ct.pointer(G))

                if OUTPUT and step_count % OUTPUT_FREQ == 0:
                    lib.get_tree_info(P.root, ct.pointer(P),  ct.pointer(G))
                    finish = time.time()
                    print("Timestep: " + str(nm) + "-" + str(step_count) + ", Program time: " + str(finish - start) + " s, Sim. time: " + str(tt + mt + rt) + " TU, Active/Total Cells: " + str(P.a_count) + "/" + str(P.tot_count) + ", Max key %: " + str(P.max_key/((2**64)-1)*100) + "%")
                step_count += 1

                rt += G.dt

            if OUTPUT and step_count % OUTPUT_FREQ != 0:
                lib.get_tree_info(P.root, ct.pointer(P),  ct.pointer(G))
                finish = time.time()
                print("Timestep: " + str(nm) + "-" + str(step_count) + ", Program time: " + str(finish - start) + " s, Sim. time: " + str(tt + mt + rt) + " TU, Active/Total Cells: " + str(P.a_count) + "/" + str(P.tot_count) + ", Max key %: " + str(P.max_key/((2**64)-1)*100) + "%")

            if RECORD:
                if OUTPUT:
                    print("\nRECORDING PDF AT: " + str(tt + mt + rt) + " TU...\n")
                FILE_NAME = FILE_PATH + "/PDFs/P" + str(nm) + "/pdf_" + str(record_count) + ".txt"
                lib.record_data(P.root, FILE_NAME.encode('utf-8'), G, ct.c_double(tt + mt + rt))
                record_count += 1
        
            mt += rt

        tt += mt
        if MEASURE and nm < NM-1:
            if OUTPUT:
                print("\nPERFORMING BAYESIAN UPDATE AT: " + str(tt) + " TU...\n")
            nm += 1

            FILE_NAME = FILE_PATH + "/Measurements/measurement" + str(nm) + ".txt"

            lib.measurement_update(FILE_NAME.encode('utf-8'), ct.pointer(P), G, ct.pointer(measure_time))

    print("Time marching complete.\n")