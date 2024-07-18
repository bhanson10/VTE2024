import numpy as np
import ctypes as ct
import sys
import time

DIM = 3  # Dimensionality

def rosenberg_pair(state, d, m):
    if d==1:
        return state[0]
    
    new_state = np.empty(d-1)
    for i in range(d-1):
        new_state[i] = state[i]

    new_m = np.max(new_state)
    return rosenberg_pair(new_state, d-1, new_m) + m**d + (m - state[d-1])*((m + 1)**(d - 1) - m**(d - 1))

def state_conversion(state):
    shift_state = np.empty(DIM)
    for i in range(DIM):
        if(state[i] < 0):
            shift_state[i] = -2*state[i] - 1
        else:
            shift_state[i] = 2*state[i]

    m = np.max(shift_state)
    return int(rosenberg_pair(shift_state, DIM, m))

def mc(th):
    return max(0, min((1 + th)/2, 2.0, 2*th))

def get_size(root):
    if root == None:
        return 0
    
    return 1 + get_size(root.left) + get_size(root.right)

class Cell(ct.Structure):
    pass

Cell._fields_ = [
        ("key", ct.c_int),
        ("prob", ct.c_float),
        ("v", ct.c_float * DIM),
        ("u", ct.c_float * DIM),
        ("w", ct.c_float * DIM),
        ("ctu", ct.c_float * DIM),
        ("state", ct.c_int * DIM),
        ("i_nodes", ct.POINTER(Cell) * DIM),
        ("k_nodes", ct.POINTER(Cell) * DIM),
        ("dcu", ct.c_float),
        ("cfl_dt", ct.c_float),
        ("vuw_f", ct.c_int),
        ("ik_f", ct.c_int),
        ("del_f", ct.c_int)
    ]

class TreeNode:
    def __init__(self, CELL):
        self.cell = CELL
        self.cell.i_nodes = (ct.POINTER(Cell) * DIM)(*[ct.POINTER(Cell)(), ct.POINTER(Cell)(), ct.POINTER(Cell)()])
        self.cell.k_nodes = (ct.POINTER(Cell) * DIM)(*[ct.POINTER(Cell)(), ct.POINTER(Cell)(), ct.POINTER(Cell)()])

        self.left = None
        self.right = None
        self.height = 1

class BST:
    def __init__(self, NODE=None):
        self.root = NODE
        self.dead = None
        self.a_count = 0 
        self.tot_count = 0 
        self.prob_sum = 0
        self.max_key = -1
        self.cfl_min_dt = sys.float_info.max

    #========================================== BST Functions ===============================================#
    def insert(self, new_node):
        self.root = self._insert_recursive(self.root, new_node)

    def _insert_recursive(self, root, new_node):
        if not root:
            return new_node
        
        if new_node.cell.key < root.cell.key:
            root.left = self._insert_recursive(root.left, new_node)
        elif new_node.cell.key > root.cell.key:
            root.right = self._insert_recursive(root.right, new_node)
        else:
            return root

        return root
    
    def _min_value_node(self, root):
        current = root
        while current.left:
            current = current.left
        return current

    def delete(self, key):
        self.root = self._delete_recursive(self.root, key)

    def _delete_recursive(self, root, key):
        if not root:
            return root
        elif key < root.cell.key:
            root.left = self._delete_recursive(root.left, key)
        elif key > root.cell.key:
            root.right = self._delete_recursive(root.right, key)
        else:
            if not root.left:
                temp = root.right
                root = None
                return temp
            elif not root.right:
                temp = root.left
                root = None
                return temp
            else:
                temp = self._min_value_node(root.right)
                root.cell = temp.cell
                root.right = self._delete_recursive(root.right, temp.cell.key)

        return root

    def _get_height(self, root):
        if not root:
            return 0
        
        return 1 + max(self._get_height(root.left), self._get_height(root.right))
    
    def _get_difference(self, root):
        l_height = self._get_height(root.left)
        r_height = self._get_height(root.right)
        b_factor = l_height - r_height
        return b_factor
    
    def _rr_rotate(self, root):
        t = root.right
        root.right = t.left
        t.left = root
        return t
    
    def _ll_rotate(self, root):
        t = root.left
        root.left = t.right
        t.right = root
        return t
    
    def _lr_rotate(self, root):
        t = root.left
        root.left = self._rr_rotate(t)
        return self._ll_rotate(root)
    
    def _rl_rotate(self, root):
        t = root.right
        root.right = self._ll_rotate(t)
        return self._rr_rotate(root)
    
    def balance(self, root):
        bal_factor = self._get_difference(root)
        while(abs(bal_factor) >  1):
            if bal_factor > 1:
                if self._get_difference(root.left) > 0:
                    root = self._ll_rotate(root)
                else:
                    root = self._lr_rotate(root)
            elif bal_factor < -1:
                if self._get_difference(root.right) > 0:
                    root = self._rl_rotate(root)
                else:
                    root = self._rr_rotate(root)
            bal_factor = self._get_difference(root)

        return root

    def search(self, key):
        return self._search_recursive(self.root, key)
    
    def _search_recursive(self, root, key):
        if not root:
            return ct.POINTER(Cell)()
        if root.cell.key == key:
            return ct.pointer(root.cell)
        if root.cell.key < key:
            return self._search_recursive(root.right, key)
        return self._search_recursive(root.left, key)

    def _write_file(self, file, G):
        return self._write_file_recursive(self.root, file, G)

    def _write_file_recursive(self, node, file, G):
        if not node:
            return
        else:
            self._write_file_recursive(node.left, file, G)
            self._write_file_recursive(node.right, file, G)

            if node.cell.prob >= G.thresh:
                file.write(str(node.cell.prob) + " ")
                for i in range(DIM):
                    file.write(str(G.dx[i]*node.cell.state[i] + G.center[i]) + " ")
                file.write("\n")
                                
    def print_tree(self):
        return self._print_tree_recursive(self.root)
    
    def _print_tree_recursive(self, root, indent=None, last=None):
        if root:
            if not indent:
                indent = ""
                
            print(indent, end="")
            if indent == "":
                print("ROOT----", end="")
                indent += "        "
            else:
                if last:
                    print("R----", end="")
                    indent += "     "
                else:
                    print("L----", end="")
                    indent += "|    "
            print(root.cell.key)
            self._print_tree_recursive(root.left, indent, False)
            self._print_tree_recursive(root.right, indent, True)
    #========================================================================================================#

    #========================================= GBEES Functions ==============================================#
    def initialize_grid(self, G, traj, m):
        blank_c = Cell(key = -1, prob = 0, v = (ct.c_float * DIM)(*[0.0,0.0,0.0]), u = (ct.c_float * DIM)(*[0.0,0.0,0.0]), w = (ct.c_float * DIM)(*[0.0,0.0,0.0]), ctu = (ct.c_float * DIM)(*[0.0,0.0,0.0]), dcu = 0)
        self.dead = TreeNode(blank_c)
        self.dead.cell.i_nodes = (ct.POINTER(Cell) * DIM)(*[ct.pointer(self.dead.cell), ct.pointer(self.dead.cell), ct.pointer(self.dead.cell)])
        self.dead.cell.k_nodes = (ct.POINTER(Cell) * DIM)(*[ct.pointer(self.dead.cell), ct.pointer(self.dead.cell), ct.pointer(self.dead.cell)])

        for i in range(round(-3*(m.cov[0][0]**0.5)/G.dx[0]), round(3*(m.cov[0][0]**0.5)/G.dx[0]) + 1):
            for j in range(round(-3*(m.cov[1][1]**0.5)/G.dx[1]), round(3*(m.cov[1][1]**0.5)/G.dx[1]) + 1):
                for k in range(round(-3*(m.cov[2][2]**0.5)/G.dx[2]), round(3*(m.cov[2][2]**0.5)/G.dx[2]) + 1):
                    current_state = np.array([i,j,k])
                    key = state_conversion(current_state)
                    current_state_vec = current_state*G.dx
                    x = current_state_vec @ (np.linalg.inv(m.cov)) @ current_state_vec
                    c = Cell(key = key, prob = np.exp(-x/2), state = (ct.c_int * DIM)(*current_state), vuw_f = 0, ik_f = 0, del_f = 0)
                    new_node = TreeNode(c)
                    self.insert(new_node)

        self.initialize_vuw(G, traj, self.root)
        self.initialize_ik_nodes(self.root)

    def initialize_vuw(self, G, traj, root):
        if not root:
            return 
        
        self.initialize_vuw(G, traj, root.left)
        self.initialize_vuw(G, traj, root.right)

        if root.cell.vuw_f == 0:
            x = G.dx*np.array(root.cell.state) + G.center

            v1 = traj.sigma*(x[1] - (x[0] + (G.dx[0]/2.0)))
            v2 = -(x[1] + (G.dx[1]/2.0)) - x[0]*x[2]
            v3 = -traj.b*(x[2] + (G.dx[2]/2.0)) + x[0]*x[1] - traj.b * traj.r 
            root.cell.v = (ct.c_float * DIM)(*[v1, v2, v3])
            root.cell.u = (ct.c_float * DIM)(*[min(v1, 0.0), min(v2, 0.0), min(v3, 0.0)])
            root.cell.w = (ct.c_float * DIM)(*[max(v1, 0.0), max(v2, 0.0), max(v3, 0.0)])
            root.cell.vuw_f = 1

            sum = 0
            for q in range(DIM):
                sum += abs(root.cell.v[q])/G.dx[q]
            
            root.cell.cfl_dt = 1/sum

    def initialize_ik_nodes(self, root):
        if not root: 
            return 
        
        self.initialize_ik_nodes(root.left)
        self.initialize_ik_nodes(root.right)

        if root.cell.ik_f == 0: 
            for q in range(DIM):
                # Initializing I Node
                i_state = np.array(root.cell.state) 
                i_state[q] -= 1
                i_key = state_conversion(i_state)
                i_node = self.search(i_key)
                if not i_node:
                    root.cell.i_nodes[q] = ct.pointer(self.dead.cell)
                else:
                    root.cell.i_nodes[q] = i_node
                    i_node.contents.k_nodes[q] = ct.pointer(root.cell)

                # Initializing K Node
                k_state = np.array(root.cell.state) 
                k_state[q] += 1
                k_key = state_conversion(k_state)
                k_node = self.search(k_key)

                if not k_node:
                    root.cell.k_nodes[q] = ct.pointer(self.dead.cell)
                else:
                    root.cell.k_nodes[q] = k_node
                    k_node.contents.i_nodes[q] = ct.pointer(root.cell)

            root.ik_f = 1
    
    def grow_tree(self, G, Lor):
        self.create_neighbors(G, self.root)
        self.initialize_ik_nodes(self.root)
        self.initialize_vuw(G, Lor, self.root)
        self.root = self.balance(self.root)
    
    def create_neighbors(self, G, root):
        if not root:
            return
        
        self.create_neighbors(G, root.left)
        self.create_neighbors(G, root.right)
        
        if root.cell.prob >= G.thresh:
            current_v = np.array(root.cell.v)
            for q in range(DIM):
                new_state = np.array(root.cell.state)

                # Checking Forward Faces
                if current_v[q] > 0:
                    if root.cell.k_nodes[q].contents.key == -1:
                        new_state[q] += 1
                        new_key = state_conversion(new_state)
                        c = Cell(key = new_key, prob = 0, state = (ct.c_int * DIM)(*new_state), vuw_f = 0, ik_f = 0, del_f = 0)
                        new_node = TreeNode(c)
                        self.insert(new_node)
                        
                        # Checking Edges
                        for e in range(DIM):
                            new_state = np.array(root.cell.state)
                            new_state[q] += 1
                            if e!=q:
                                if current_v[e] > 0:
                                    new_state[e] += 1
                                    new_key = state_conversion(new_state)
                                    c = Cell(key = new_key, prob = 0, state = (ct.c_int * DIM)(*new_state), vuw_f = 0, ik_f = 0, del_f = 0)
                                    new_node = TreeNode(c)
                                    self.insert(new_node)

                                elif current_v[e] < 0:
                                    new_state[e] -= 1
                                    new_key = state_conversion(new_state)
                                    c = Cell(key = new_key, prob = 0, state = (ct.c_int * DIM)(*new_state), vuw_f = 0, ik_f = 0, del_f = 0)
                                    new_node = TreeNode(c)
                                    self.insert(new_node)
                    else:
                        # Checking Edges
                        for e in range(DIM):
                            new_state = np.array(root.cell.state)
                            new_state[q] += 1
                            if e!=q:
                                if current_v[e] > 0:
                                    if root.cell.k_nodes[q].contents.k_nodes[e].contents.key == -1:
                                        new_state[e] += 1
                                        new_key = state_conversion(new_state)
                                        c = Cell(key = new_key, prob = 0, state = (ct.c_int * DIM)(*new_state), vuw_f = 0, ik_f = 0, del_f = 0)
                                        new_node = TreeNode(c)
                                        self.insert(new_node)

                                elif current_v[e] < 0:
                                    if root.cell.k_nodes[q].contents.i_nodes[e].contents.key == -1:
                                        new_state[e] -= 1
                                        new_key = state_conversion(new_state)
                                        c = Cell(key = new_key, prob = 0, state = (ct.c_int * DIM)(*new_state), vuw_f = 0, ik_f = 0, del_f = 0)
                                        new_node = TreeNode(c)
                                        self.insert(new_node)

                # Checking Backward Faces
                elif current_v[q] < 0:
                    if root.cell.i_nodes[q].contents.key == -1:
                        new_state[q] -= 1
                        new_key = state_conversion(new_state)
                        c = Cell(key = new_key, prob = 0, state = (ct.c_int * DIM)(*new_state), vuw_f = 0, ik_f = 0, del_f = 0)
                        new_node = TreeNode(c)
                        self.insert(new_node)
                        
                        # Checking Edges
                        for e in range(DIM):
                            new_state = np.array(root.cell.state)
                            new_state[q] -= 1
                            if e!=q:
                                if current_v[e] > 0:
                                    new_state[e] += 1
                                    new_key = state_conversion(new_state)
                                    c = Cell(key = new_key, prob = 0, state = (ct.c_int * DIM)(*new_state), vuw_f = 0, ik_f = 0, del_f = 0)
                                    new_node = TreeNode(c)
                                    self.insert(new_node)

                                elif current_v[e] < 0:
                                    new_state[e] -= 1
                                    new_key = state_conversion(new_state)
                                    c = Cell(key = new_key, prob = 0, state = (ct.c_int * DIM)(*new_state), vuw_f = 0, ik_f = 0, del_f = 0)
                                    new_node = TreeNode(c)
                                    self.insert(new_node)
                    else:
                        # Checking Edges
                        for e in range(DIM):
                            new_state = np.array(root.cell.state)
                            new_state[q] -= 1
                            if e!=q:
                                if current_v[e] > 0:
                                    if root.cell.k_nodes[q].contents.k_nodes[e].contents.key == -1:
                                        new_state[e] += 1
                                        new_key = state_conversion(new_state)
                                        c = Cell(key = new_key, prob = 0, state = (ct.c_int * DIM)(*new_state), vuw_f = 0, ik_f = 0, del_f = 0)
                                        new_node = TreeNode(c)
                                        self.insert(new_node)

                                elif current_v[e] < 0:
                                    if root.cell.k_nodes[q].contents.i_nodes[e].contents.key == -1:
                                        new_state[e] -= 1
                                        new_key = state_conversion(new_state)
                                        c = Cell(key = new_key, prob = 0, state = (ct.c_int * DIM)(*new_state), vuw_f = 0, ik_f = 0, del_f = 0)
                                        new_node = TreeNode(c)
                                        self.insert(new_node)

    def prune_tree(self, G):
        self._mark_cells(G, self.root)
        self._delete_cells(G, self.root)
        self.normalize_tree(G)
        self.initialize_ik_nodes(self.root)

    def _mark_cells(self, G, root):
        if not root:
            return
        
        self._mark_cells(G, root.left)
        self._mark_cells(G, root.right)

        root.cell.ik_f = 0
        DELETE = True

        if root.cell.prob < G.thresh:
            for q in range(DIM):
                # Looking at Backwards Node
                if root.cell.i_nodes[q].contents.key != -1:
                    if (root.cell.i_nodes[q].contents.v[q] > 0) and (root.cell.i_nodes[q].contents.prob >= G.thresh):
                        DELETE = False
                        break
                    else:
                        for e in range(DIM):
                            if e!=q:
                                if root.cell.i_nodes[q].contents.i_nodes[e].contents.v[q] > 0 and root.cell.i_nodes[q].contents.i_nodes[e].contents.v[e] > 0 and root.cell.i_nodes[q].contents.i_nodes[e].contents.prob >= G.thresh:
                                    DELETE = False
                                    break

                                if root.cell.i_nodes[q].contents.k_nodes[e].contents.v[q] > 0 and root.cell.i_nodes[q].contents.k_nodes[e].contents.v[e] < 0 and root.cell.i_nodes[q].contents.k_nodes[e].contents.prob >= G.thresh:
                                    DELETE = False
                                    break

                # Looking at Forwards Node
                if root.cell.k_nodes[q].contents.key != -1:
                    if (root.cell.k_nodes[q].contents.v[q] < 0) and (root.cell.k_nodes[q].contents.prob >= G.thresh):
                        DELETE = False
                        break
                    else:
                        for e in range(DIM):
                            if e!=q:
                                if root.cell.k_nodes[q].contents.i_nodes[e].contents.v[q] < 0 and root.cell.k_nodes[q].contents.i_nodes[e].contents.v[e] > 0 and root.cell.k_nodes[q].contents.i_nodes[e].contents.prob >= G.thresh:
                                    DELETE = False
                                    break

                                if root.cell.k_nodes[q].contents.k_nodes[e].contents.v[q] < 0 and root.cell.k_nodes[q].contents.k_nodes[e].contents.v[e] < 0 and root.cell.k_nodes[q].contents.k_nodes[e].contents.prob >= G.thresh:
                                    DELETE = False
                                    break
            
            if DELETE:
                root.cell.del_f = 1
    
    def _delete_cells(self, G, root):
        if not root:
            return

        self._delete_cells(G, root.left)
        self._delete_cells(G, root.right)

        if root.cell.del_f == 1:
            self.delete(root.cell.key)
                
    def godunov_method(self, G):
        self._get_dcu(G, self.root)
        self._update_ctu(G, self.root)

    def _get_dcu(self, G, root):
        if not root:
            return
        
        self._get_dcu(G, root.left)
        self._get_dcu(G, root.right)

        root.cell.dcu = 0
        root.cell.ctu = (ct.c_float * DIM)(*[0.0,0.0,0.0])

        for q in range(DIM):
            i_node = root.cell.i_nodes[q]

            dcu_p = root.cell.w[q] * root.cell.prob + root.cell.u[q] * root.cell.k_nodes[q].contents.prob
            dcu_m = i_node.contents.w[q] * i_node.contents.prob + i_node.contents.u[q] * root.cell.prob

            root.cell.dcu -= (G.dt/G.dx[q])*(dcu_p - dcu_m)
    
    def _update_ctu(self, G, root):
        if not root:
            return 
        
        self._update_ctu(G, root.left)
        self._update_ctu(G, root.right)

        for q in range(DIM):
            i_node = root.cell.i_nodes[q]
            if i_node.contents.key != -1:
                F = G.dt*(root.cell.prob - i_node.contents.prob)/(2*G.dx[q])
                for e in range(DIM):
                    if e!=q:
                        j_node = root.cell.i_nodes[e]
                        p_node = i_node.contents.i_nodes[e]

                        root.cell.ctu[e]       -= i_node.contents.w[q]*root.cell.w[e]*F
                        j_node.contents.ctu[e] -= i_node.contents.w[q]*j_node.contents.u[e]*F
                        i_node.contents.ctu[e] -= i_node.contents.u[q]*i_node.contents.w[e]*F
                        p_node.contents.ctu[e] -= i_node.contents.u[q]*p_node.contents.u[e]*F

                # High-resolution correction terms
                if i_node.contents.v[q] > 0:
                    i_i_node = i_node.contents.i_nodes[q]
                    if root.cell.prob - i_node.contents.prob == 0:
                        th = float("nan")
                    else:
                        th = (i_node.contents.prob - i_i_node.contents.prob)/(root.cell.prob - i_node.contents.prob)
                else:
                    if root.cell.prob - i_node.contents.prob == 0:
                        th = float("nan")
                    else:
                        th = (root.cell.k_nodes[q].contents.prob - root.cell.prob)/(root.cell.prob - i_node.contents.prob)

                i_node.contents.ctu[q] += abs(i_node.contents.v[q])*(G.dx[q]/G.dt - abs(i_node.contents.v[q]))*F*mc(th)

    def update_prob(self, G, root):
        if not root:
            return
        
        self.update_prob(G, root.left)
        self.update_prob(G, root.right)

        root.cell.prob += root.cell.dcu
        for q in range(DIM):
            root.cell.prob -= (G.dt/G.dx[q])*(root.cell.ctu[q] - root.cell.i_nodes[q].contents.ctu[q])

    def normalize_tree(self, G):
        self.a_count = 0 
        self.tot_count = 0
        self.prob_sum = 0

        self._get_sum(G, self.root)
        self._divide_sum(G, self.root)

    def _get_sum(self, G, root):
        if not root:
            return
        
        self._get_sum(G, root.left)
        self._get_sum(G, root.right)

        self.prob_sum += root.cell.prob
    
    def _divide_sum(self, G, root):
        if not root:
            return

        self._divide_sum(G, root.left)
        self._divide_sum(G, root.right)

        root.cell.prob /= self.prob_sum

    def record_data(self, file_name, G, t):
        file = open(file_name, "w")
        file.write(str(t) + "\n")
        self._write_file(file, G)
        file.close()

    def get_max_key(self, G):
        self.max_key = 0
        self.a_count = 0
        self.tot_count = 0
        self._get_max_key_recursive(self.root, G)

    def _get_max_key_recursive(self, root, G):
        if not root:
            return
        
        self._get_max_key_recursive(root.left, G)
        self._get_max_key_recursive(root.right, G)

        self.max_key = max(self.max_key, root.cell.key)

        if root.cell.prob >= G.thresh:
            self.a_count += 1
        
        self.tot_count += 1

    def check_cfl_condition(self, root):
        if not root:
            return
        
        self.check_cfl_condition(root.left)
        self.check_cfl_condition(root.right)
        
        self.cfl_min_dt = min(self.cfl_min_dt, root.cell.cfl_dt)
    #========================================================================================================#

class Grid:
    def __init__(self):
        self.thresh = None
        self.dt     = None
        self.center = None
        self.dx     = None

class Traj: 
    sigma = None
    b     = None
    r     = None

class Measurement:
    def __init__(self):
        self.mean = None
        self.cov  = None
        self.T    = None

if __name__ == "__main__":
     #===================================== Read in measurement/trajectory info =================================#
    print("Reading in discrete measurements...\n")

    NM        = 1         # Number of measurements 
    FILE_PATH = "./Data"  # Measurement file path

    m = Measurement()
    measurement_file = open(FILE_PATH + "/Measurements/measurement0.txt", "r")
    measurement_file.readline(); # Skip label line
    m.mean = np.array([float(string) for string in measurement_file.readline().split(" ")]) # Mean
    measurement_file.readline(); # Skip blank space
    measurement_file.readline(); # Skip label line

    m.cov = np.empty((DIM,DIM))
    for i in range(DIM):
        m.cov[i] = [float(string) for string in measurement_file.readline().split(" ")] # Covariance
    measurement_file.readline(); # Skip blank space
    measurement_file.readline(); # Skip label line
    m.T = float(measurement_file.readline())
    measurement_file.close()
    #===========================================================================================================#
    
    #========================================== Begin User Input ===============================================#
    print( "Reading in user inputs...\n" )

    G           = Grid()  # Grid instance
    G.thresh    = 2E-5    # Probability threshold
    G.dx = np.array([x**0.5/2 for x in np.diag(m.cov)]) # Grid width
    G.center = m.mean

    Lor = Traj() # Trajectory instance
    Lor.sigma = 4
    Lor.b = 1
    Lor.r = 48

    OUTPUT      = True    # Write info to terminal
    RECORD      = True    # Write PDFs to .txt file
    MEASURE     = False   # Take discrete measurement updates
    OUTPUT_FREQ = 20      # Number of steps per output to terminal
    DEL_STEP    = 25      # Number of steps per deletion procedure
    NUM_DIST    = 6       # Number of distributions recorded per measurement
    record_time = m.T/(NUM_DIST-1)
    #===========================================================================================================#

    #=============================================== GBEES =====================================================#
    P = BST()

    print("Initializing distribution...\n")

    P.initialize_grid(G, Lor, m)
    P.normalize_tree(G)

    print("Entering time marching...\n")

    start = time.time()
    tt = 0
    for nm in range(NM):

        finish = time.time()
        P.get_max_key(G)
        print("Timestep: " + str(nm) + "-0, Program time: " + str(finish - start) + " s, Sim. time: " + str(tt) + " TU, Active/Total Cells: " + str(P.a_count) + "/" + str(P.tot_count) + ", Max key %: " + str(P.max_key/((2**64)-1)*100) + "%")
        if RECORD:
            file_name = FILE_PATH + "/PDFs/P" + str(nm) + "/pdf_0.txt"
            P.record_data(file_name, G, 0)

        mt = 0
        record_count = 1
        step_count = 1
        while(mt < m.T):

            rt = 0
            while(rt < record_time):
                
                P.grow_tree(G, Lor)
                P.check_cfl_condition(P.root)
                G.dt = min(P.cfl_min_dt, record_time - rt)
                rt += G.dt
                P.godunov_method(G)
                P.update_prob(G, P.root)
                P.normalize_tree(G)

                if (step_count % DEL_STEP == 0):
                    P.prune_tree(G)

                step_count += 1
                if OUTPUT and step_count % OUTPUT_FREQ == 0:
                    finish = time.time()
                    P.get_max_key(G)
                    print("Timestep: " + str(nm) + "-" + str(step_count) + ", Program time: " + str(finish - start) + " s, Sim. time: " + str(tt + mt + rt) + " TU, Active/Total Cells: " + str(P.a_count) + "/" + str(P.tot_count) + ", Max key %: " + str(P.max_key/((2**64)-1)*100) + "%")

                P.cfl_min_dt = sys.float_info.max

            if OUTPUT and step_count % OUTPUT_FREQ != 0:
                finish = time.time()
                P.get_max_key(G)
                print("Timestep: " + str(nm) + "-" + str(step_count) + ", Program time: " + str(finish - start) + " s, Sim. time: " + str(tt + mt + rt) + " TU, Active/Total Cells: " + str(P.a_count) + "/" + str(P.tot_count) + ", Max key %: " + str(P.max_key/((2**64)-1)*100) + "%")

            if RECORD:
                print("\nRECORDING PDF AT: " + str(tt + mt + rt) + " TU...\n")
                file_name = FILE_PATH + "/PDFs/P" + str(nm) + "/pdf_" + str(record_count) + ".txt"
                P.record_data(file_name, G, tt + mt + rt)
                record_count += 1
        
            mt += rt

        tt += mt
        if MEASURE and nm < NM-1:
            print("\nPERFORMING BAYESIAN UPDATE AT: " + str(tt) + " TU...\n")
            nm += 1

            m = Measurement()
            measurement_file = open(FILE_PATH + "/Measurements/measurement" + nm + ".txt", "r")
            measurement_file.readline(); # Skip label line
            m.mean = np.array([float(string) for string in measurement_file.readline().split(" ")]) # Mean
            measurement_file.readline(); # Skip blank space
            measurement_file.readline(); # Skip label line
            for i in range(DIM):
                m.cov[i] = [float(string) for string in measurement_file.readline().split(" ")] # Covariance
            measurement_file.readline(); # Skip blank space
            measurement_file.readline(); # Skip label line
            m.T = float(measurement_file.readline())
            measurement_file.close()