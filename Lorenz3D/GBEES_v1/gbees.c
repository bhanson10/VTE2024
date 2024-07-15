#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <stdio.h>
#include <lapacke.h>
#include <cblas.h>
#include <stdbool.h>
#include <string.h>

#define DIM 3

/*==============================================================================
Non-member Function DEFINITIONS
==============================================================================*/
uint64_t rosenberg_pair(const int* state, int d, int m) {
    if (d == 1) {
        return state[0];
    }

    int* new_state = (int*)malloc((d-1) * sizeof(int));
    for (int i = 0; i < d-1; i++) {
        new_state[i] = state[i];
    }

    int new_m = new_state[0];
    for (int i = 1; i < d-1; i++) {
        if (new_state[i] > new_m) {
            new_m = new_state[i];
        }
    }

    uint64_t result = rosenberg_pair(new_state, d-1, new_m) + (uint64_t)pow(m, d) + (m - state[d-1]) * ((uint64_t)pow(m+1, d-1) - (uint64_t)pow(m, d-1));
    free(new_state);
    return result;
}

uint64_t state_conversion(const int* state) {
    int shift_state[DIM];
    int m;
    uint64_t key;

    for (int i = 0; i < DIM; i++) {
        if (state[i] < 0) {
            shift_state[i] = -2 * state[i] - 1;
        } else {
            shift_state[i] = 2 * state[i];
        }
    }

    m = shift_state[0];
    for (int i = 1; i < DIM; i++) {
        if (shift_state[i] > m) {
            m = shift_state[i];
        }
    }

    key = rosenberg_pair(shift_state, DIM, m);
    return key;
}

double mc(double th){
    double min1 = fmin((1 + th)/2.0, 2.0);
    return fmax(0.0, fmin(min1, 2*th)); 
}

double gauss_probability(double* x, double* mat){
    lapack_int n = DIM;
    lapack_int ipiv[DIM];
    lapack_int info;
    info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, mat, n, ipiv);
    info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, mat, n, ipiv);
    double y[DIM];
    cblas_dgemv(CblasRowMajor, CblasNoTrans, DIM, DIM, 1.0, mat, DIM, x, 1, 0.0, y, 1);
    return exp(-cblas_ddot(DIM, x, 1, y, 1)/2);
}
/*==============================================================================
CLASS DEFINITIONS
==============================================================================*/
typedef struct TreeNode TreeNode;

typedef struct Cell { 
    double prob;
    double v[DIM];
    double u[DIM];
    double w[DIM];
    double ctu[DIM];
    int state[DIM];
    TreeNode* i_nodes[DIM];
    TreeNode* k_nodes[DIM];
    double dcu;
    double cfl_dt;
    int new_f;
    int ik_f;
    int del_f;
} Cell;

Cell Cell_create(double prob, double* v, double* u, double* w, double* ctu, int* state, double dcu, int new_f, int ik_f, int del_f) {
    Cell c;
    c.prob = prob;
    for(int i = 0; i < DIM; i++){
        c.v[i] = v[i]; c.u[i] = u[i]; c.w[i] = w[i]; c.ctu[i] = ctu[i]; c.state[i] = state[i]; 
    }
    c.dcu = dcu; 
    c.new_f = new_f; 
    c.ik_f = ik_f; 
    c.del_f = del_f; 

    return c;
}

typedef struct Grid {
    double thresh;
    double dt;
    double center[DIM];
    double dx[DIM];
} Grid;

typedef struct Traj {
    double sigma;
    double b;
    double r;
} Traj;

typedef struct TreeNode {
    uint64_t key;
    Cell cell;
    struct TreeNode* left;
    struct TreeNode* right;
} TreeNode;

TreeNode* TreeNode_create(uint64_t k, Cell c) {
    TreeNode* node = (TreeNode*)malloc(sizeof(TreeNode));
    node->key = k;
    node->cell = c;
    node->left = NULL;
    node->right = NULL;
    return node;
}

typedef struct Measurement{
    double mean[DIM];
    double cov[DIM][DIM];
    double T;
} Measurement;

Measurement Measurement_create() {
    Measurement M;
    for (int i = 0; i < DIM; i++) {
        M.mean[i] = 0.0; 
        for (int j = 0; j < DIM; j++) {
            M.cov[i][j] = 0.0;
        }
    }
    M.T = 0; 
    return M;
}

typedef struct BST{
    TreeNode* dead;
    TreeNode* root;
    int a_count;
    int tot_count;
    uint64_t max_key;
    double cfl_min_dt;
 } BST;

BST BST_create() {
    BST P;
    P.dead = NULL;
    P.root = NULL; 
    P.a_count = 0;
    P.tot_count = 0;
    P.max_key = (uint64_t) 0;
    P.cfl_min_dt = INT_MAX;

    return P; 
}

TreeNode* BST_insert_recursive(TreeNode* r, TreeNode* new_node){
    if(r == NULL){
        r = new_node; 
        return r; 
    }

    if (new_node->key < r->key){
        r->left = BST_insert_recursive(r->left, new_node);
    }else if (new_node->key > r->key){
        r->right = BST_insert_recursive(r->right, new_node);
    }else{
        return r;
    }
    return r;
}

TreeNode* BST_search_recursive(TreeNode* r, uint64_t k){ 
    if ((r == NULL)||(r->key == k)){
        return r;
    }else if (k < r->key){
        return BST_search_recursive(r->left, k);
    }else{
        return BST_search_recursive(r->right,k);
    }
}

TreeNode* min_value_node(TreeNode* node){ 
    TreeNode* current = node;
    while(current->left != NULL){
        current = current->left;
    }
    return current;
}

TreeNode* BST_delete_node(TreeNode* r, uint64_t k){ 
    if(r==NULL){
        return NULL;
    }else if(k < r->key){
        r->left = BST_delete_node(r->left, k);
    }else if(k > r->key){
        r->right = BST_delete_node(r->right,k);
    }else{
        if(r->left == NULL){
            TreeNode* temp = r->right;
            free(r);
            return temp;
        } else if (r->right == NULL){
            TreeNode* temp = r->left;
            free(r);
            return temp;
        }else{
            TreeNode* temp = min_value_node(r->right);
            r->key = temp->key; r->cell = temp->cell;
            r->right = BST_delete_node(r->right, temp->key);
        }
    }
    return r;
}

int BST_get_height(TreeNode* r){
    if (r == NULL){
        return 0;
    }

    return 1 + fmax(BST_get_height(r->left), BST_get_height(r->right));
}

int BST_get_difference(TreeNode* r){
    int l_height = BST_get_height(r->left);
    int r_height = BST_get_height(r->right);
    int b_factor = l_height - r_height;
    return b_factor;
}

TreeNode* BST_rr_rotate(TreeNode* parent){
    TreeNode* t;
    t = parent->right;
    parent->right = t->left;
    t->left = parent;
    return t;
}

TreeNode* BST_ll_rotate(TreeNode* parent){
    TreeNode* t;
    t = parent->left;
    parent->left = t->right;
    t->right = parent;
    return t;
}

TreeNode* BST_lr_rotate(TreeNode* parent){
    TreeNode* t;
    t = parent->left;
    parent->left = BST_rr_rotate(t);
    return BST_ll_rotate(parent);
}

TreeNode* BST_rl_rotate(TreeNode* parent){
    TreeNode* t;
    t = parent->right;
    parent->right = BST_ll_rotate(t);
    return BST_rr_rotate(parent);
}

TreeNode* BST_balance(TreeNode* r){
    int bal_factor = BST_get_difference(r);
    while(abs(bal_factor) > 1){
        if (bal_factor > 1){
            if (BST_get_difference(r->left) > 0){
                r = BST_ll_rotate(r);
            }else{
                r = BST_lr_rotate(r);
            }
        }else if (bal_factor < -1){
            if (BST_get_difference(r->right) > 0){
                r = BST_rl_rotate(r);
            }else{
                r = BST_rr_rotate(r);
            }
        }
        bal_factor = BST_get_difference(r);
    }
    return r;
}

void write_file(FILE* myfile, TreeNode* r, Grid G){
    if (r == NULL){
        return;
    }

    write_file(myfile, r->left, G);
    write_file(myfile, r->right, G);

    if (r->cell.prob >= G.thresh) {
        fprintf(myfile, "%lf", r->cell.prob);
        for (int i = 0; i < DIM; i++) {
            fprintf(myfile, " %lf", G.dx[i] * r->cell.state[i] + G.center[i]);
        }
        fprintf(myfile, "\n");
    }
}

void BST_get_sum(TreeNode* r, double* prob_sum){
    if (r == NULL){
        return;
    }
    BST_get_sum(r->left, prob_sum);
    BST_get_sum(r->right, prob_sum);

    *prob_sum += r->cell.prob;
}

void BST_divide_sum(TreeNode* r, BST* P, Grid* G, double prob_sum){
    if (r == NULL){
        return;
    }
    BST_divide_sum(r->left, P, G, prob_sum);
    BST_divide_sum(r->right, P, G, prob_sum);

    r->cell.prob /= prob_sum;
}

void BST_normalize_tree(TreeNode* r, BST* P, Grid* G){
    double prob_sum = 0;
    BST_get_sum(r, &prob_sum);
    BST_divide_sum(r, P, G, prob_sum);
}

void BST_initialize_ik_nodes(TreeNode* r, BST* P){
    if(r == NULL){
        return;
    }
    BST_initialize_ik_nodes(r->left, P); 
    BST_initialize_ik_nodes(r->right, P);

    if(r->cell.ik_f == 0){
        int l_state[DIM]; 
        memcpy(l_state, r->cell.state, DIM * sizeof(int)); 
        for(int q = 0; q < DIM; q++){
            // Initializing i, k nodes
            int i_state[DIM]; memcpy(i_state, l_state, DIM * sizeof(int)); i_state[q] = i_state[q] - 1; uint64_t i_key = state_conversion(i_state); 
            int k_state[DIM]; memcpy(k_state, l_state, DIM * sizeof(int)); k_state[q] = k_state[q] + 1; uint64_t k_key = state_conversion(k_state); 
            TreeNode* i_node = BST_search_recursive(P->root, i_key); TreeNode* k_node = BST_search_recursive(P->root, k_key); 

            if(i_node == NULL){
                i_node = P->dead; r->cell.i_nodes[q] = i_node; 
            }else{
                r->cell.i_nodes[q] = i_node; i_node->cell.k_nodes[q] = r; 
            }

            if(k_node == NULL){
                k_node = P->dead; r->cell.k_nodes[q] = k_node; 
            }else{
                r->cell.k_nodes[q] = k_node; k_node->cell.i_nodes[q] = r; 
            }
        }
        r->cell.ik_f = 1; 
    } 
}

void BST_initialize_grid(BST* P, Grid* G, Measurement M){
    double zeros[DIM] = {0.0};
    int max[DIM] = {INT_MAX, INT_MAX, INT_MAX};
    TreeNode* dead_node = TreeNode_create(-1, Cell_create(0, zeros, zeros, zeros, zeros, max, 0, -1, -1, -1)); 
    
    for(int i = 0; i < DIM; i++){
        dead_node->cell.i_nodes[i] = dead_node; dead_node->cell.k_nodes[i] = dead_node; 
    }
    P->dead = dead_node; 

    int current_state[DIM]; double current_state_vec[DIM]; uint64_t key; double x; TreeNode* new_node; 
    for (int i = (int) round(-3*M.cov[0][0])/G->dx[0]; i <= (int) round(3*M.cov[0][0])/G->dx[0]; i++){current_state[0] = i; current_state_vec[0] = i*G->dx[0]; 
        for (int j = (int) round(-3*M.cov[1][1])/G->dx[1]; j <= (int) round(3*M.cov[1][1])/G->dx[1]; j++){current_state[1] = j; current_state_vec[1] = j*G->dx[1];
            for (int k = (int) round(-3*M.cov[2][2])/G->dx[2]; k <= (int) round(3*M.cov[2][2])/G->dx[2]; k++){current_state[2] = k; current_state_vec[2] = k*G->dx[2];
                key = state_conversion(current_state);
                x = gauss_probability(current_state_vec, (double *)M.cov);
                new_node = TreeNode_create(key, Cell_create(x, zeros, zeros, zeros, zeros, current_state, 0, 0, 0, 0));
                P->root = BST_insert_recursive(P->root, new_node);
            }
        }
    }

    BST_initialize_ik_nodes(P->root, P);
}


void BST_create_neighbors(TreeNode* r, BST* P, Grid G){
    if (r == NULL){
        return;
    }
    BST_create_neighbors(r->left, P, G);
    BST_create_neighbors(r->right, P, G);

    if (r->cell.prob >= G.thresh){
        Cell c;
        TreeNode* new_node;
        double current_v[DIM];
        memcpy(current_v, r->cell.v, DIM * sizeof(double));
        int current_state[DIM]; 
        memcpy(current_state, r->cell.state, DIM * sizeof(int));
        int new_state[DIM]; uint64_t new_key;
        double zeros[DIM] = {0.0};
        for (int q = 0; q < DIM; q++){
            memcpy(new_state, current_state, DIM * sizeof(int));
            // Checking Forward Faces
            if(current_v[q] > 0){
                if(r->cell.k_nodes[q] == P->dead){
                    new_state[q] += 1;
                    new_key = state_conversion(new_state);
                    new_node = TreeNode_create(new_key, Cell_create(0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0, 0)); 
                    P->root = BST_insert_recursive(P->root, new_node);

                    // Checking Edges
                    for (int e = 0; e < DIM; e++){
                        memcpy(new_state, current_state, DIM * sizeof(int));
                        new_state[q] += 1;
                        if(e != q){
                            if(current_v[e] > 0){
                                new_state[e] += 1;
                                new_key = state_conversion(new_state);
                                new_node = TreeNode_create(new_key, Cell_create(0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0, 0)); 
                                P->root = BST_insert_recursive(P->root, new_node);

                            }else if (current_v[e] < 0){
                                new_state[e] -= 1;
                                new_key = state_conversion(new_state);
                                new_node = TreeNode_create(new_key, Cell_create(0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0, 0)); 
                                P->root = BST_insert_recursive(P->root, new_node);

                            }
                        }
                    }
                }else{
                    // Checking Edges
                    for (int e = 0; e < DIM; e++){
                        memcpy(new_state, current_state, DIM * sizeof(int));
                        new_state[q] += 1;
                        if(e != q){
                            if(current_v[e] > 0){
                                if(r->cell.k_nodes[q]->cell.k_nodes[e] == P->dead){
                                    new_state[e] += 1;
                                    new_key = state_conversion(new_state);
                                    new_node = TreeNode_create(new_key, Cell_create(0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0, 0)); 
                                    P->root = BST_insert_recursive(P->root, new_node);

                                }
                            }else if (current_v[e] < 0){
                                if(r->cell.k_nodes[q]->cell.i_nodes[e] == P->dead){
                                    new_state[e] -= 1;
                                    new_key = state_conversion(new_state);
                                    new_node = TreeNode_create(new_key, Cell_create(0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0, 0)); 
                                    P->root = BST_insert_recursive(P->root, new_node);

                                }
                            }
                        }
                    }
                }
                // Checking Backward Faces
            }else if(current_v[q] < 0){
                if(r->cell.i_nodes[q] == P->dead){
                    new_state[q] -= 1;
                    new_key = state_conversion(new_state);
                    new_node = TreeNode_create(new_key, Cell_create(0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0, 0)); 
                    P->root = BST_insert_recursive(P->root, new_node);

                    // Checking Edges
                    for (int e = 0; e < DIM; e++){
                        memcpy(new_state, current_state, DIM * sizeof(int));
                        new_state[q] -= 1;
                        if(e != q){
                            if(current_v[e] > 0){
                                new_state[e] += 1;
                                new_key = state_conversion(new_state);
                                new_node = TreeNode_create(new_key, Cell_create(0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0, 0)); 
                                P->root = BST_insert_recursive(P->root, new_node);

                            }else if (current_v[e] < 0){
                                new_state[e] -= 1;
                                new_key = state_conversion(new_state);
                                new_node = TreeNode_create(new_key, Cell_create(0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0, 0)); 
                                P->root = BST_insert_recursive(P->root, new_node);

                            }
                        }
                    }
                }else{
                    // Checking Edges
                    for (int e = 0; e < DIM; e++){
                        memcpy(new_state, current_state, DIM * sizeof(int));
                        new_state[q] -= 1;
                        if(e != q){
                            if(current_v[e] > 0){
                                if(r->cell.i_nodes[q]->cell.k_nodes[e] == P->dead){
                                    new_state[e] += 1;
                                    new_key = state_conversion(new_state);
                                    new_node = TreeNode_create(new_key, Cell_create(0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0, 0)); 
                                    P->root = BST_insert_recursive(P->root, new_node);

                                }
                            }else if (current_v[e] < 0){
                                if(r->cell.i_nodes[q]->cell.i_nodes[e] == P->dead){
                                    new_state[e] -= 1;
                                    new_key = state_conversion(new_state);
                                    new_node = TreeNode_create(new_key, Cell_create(0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0, 0)); 
                                    P->root = BST_insert_recursive(P->root, new_node);

                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void BST_check_cfl_condition(TreeNode* r, BST* P){
    if (r == NULL){
        return;
    }

    BST_check_cfl_condition(r->left, P);
    BST_check_cfl_condition(r->right, P);

    P->cfl_min_dt = fmin(P->cfl_min_dt,r->cell.cfl_dt);
}

void BST_get_dcu(TreeNode* r, Grid* G){
    if (r == NULL){
        return;
    }
    BST_get_dcu(r->left, G);
    BST_get_dcu(r->right, G);

    r->cell.dcu = 0; 
    double zeros[DIM] = {0.0, 0.0, 0.0};
    memcpy(r->cell.ctu, zeros, DIM*sizeof(double)); 
    r->cell.ctu[0] = 0.0; r->cell.ctu[1] = 0.0; r->cell.ctu[2] = 0.0; 
    for(int q = 0; q < DIM; q++){
        TreeNode* i_node = r->cell.i_nodes[q];

        double dcu_p = r->cell.w[q] * r->cell.prob + r->cell.u[q] * r->cell.k_nodes[q]->cell.prob;
        double dcu_m = i_node->cell.w[q]*i_node->cell.prob + i_node->cell.u[q]*r->cell.prob;

        r->cell.dcu -= (G->dt/G->dx[q])*(dcu_p-dcu_m);
    }
}

void BST_update_ctu(TreeNode* r, BST* P, Grid* G){
    if (r == NULL){
        return;
    }
    BST_update_ctu(r->left, P, G);
    BST_update_ctu(r->right, P, G);

    for(int a = 0; a < DIM; a++){
        TreeNode* i_node = r->cell.i_nodes[a];
        TreeNode* j_node; TreeNode* p_node;
        if(i_node!=P->dead){
            double F = G->dt*(r->cell.prob-i_node->cell.prob)/(2*G->dx[a]);
            for(int b = 0; b < DIM; b++){
                if (b!=a){
                    j_node = r->cell.i_nodes[b];
                    p_node = i_node->cell.i_nodes[b];

                    r->cell.ctu[b]      -= i_node->cell.w[a] * r->cell.w[b] * F;
                    j_node->cell.ctu[b] -= i_node->cell.w[a] * j_node->cell.u[b] * F;
                    i_node->cell.ctu[b] -= i_node->cell.u[a] * i_node->cell.w[b] * F;
                    p_node->cell.ctu[b] -= i_node->cell.u[a] * p_node->cell.u[b] * F;
                }
            }

            // High-Resolution Correction Terms
            double th;
            if (i_node->cell.v[a]>0){
                TreeNode* i_i_node = i_node->cell.i_nodes[a];
                th = (i_node->cell.prob-i_i_node->cell.prob)/(r->cell.prob-i_node->cell.prob);
            }else{
                th = (r->cell.k_nodes[a]->cell.prob-r->cell.prob)/(r->cell.prob-i_node->cell.prob);
            }

            i_node->cell.ctu[a] += fabs(i_node->cell.v[a])*(G->dx[a]/G->dt - fabs(i_node->cell.v[a]))*F*mc(th);
        }
    }
}

void BST_godunov_method(BST* P, Grid* G){
    BST_get_dcu(P->root, G);
    BST_update_ctu(P->root, P, G);
}

void BST_update_prob(TreeNode* r, Grid* G){
    if (r == NULL){
        return;
    }

    BST_update_prob(r->left, G);
    BST_update_prob(r->right, G);

    r->cell.prob += r->cell.dcu;
    for(int q = 0; q < DIM; q++){
        r->cell.prob -= (G->dt/G->dx[q])*(r->cell.ctu[q]-r->cell.i_nodes[q]->cell.ctu[q]);
    }
}


void BST_mark_cells(TreeNode* r, BST* P, Grid* G){
    if (r == NULL){
        return;
    }
    BST_mark_cells(r->left, P, G);
    BST_mark_cells(r->right, P, G);

    r->cell.ik_f = 0; bool DELETE = true;
    if (r->cell.prob < G->thresh){

        for(int q = 0; q < DIM; q++){
            // Looking at Backwards Node

            if(r->cell.i_nodes[q] != P->dead){
                if ((r->cell.i_nodes[q]->cell.v[q]>0)&&(r->cell.i_nodes[q]->cell.prob >= G->thresh)){
                    DELETE = false;
                    break;
                }else{
                    for (int e = 0; e < DIM; e++){
                        if(e!=q){
                            if ((r->cell.i_nodes[q]->cell.i_nodes[e]->cell.v[q]>0)&&(r->cell.i_nodes[q]->cell.i_nodes[e]->cell.v[e]>0)&&(r->cell.i_nodes[q]->cell.i_nodes[e]->cell.prob >= G->thresh)){
                                DELETE = false;
                                break;
                            }

                            if ((r->cell.i_nodes[q]->cell.k_nodes[e]->cell.v[q]>0)&&(r->cell.i_nodes[q]->cell.k_nodes[e]->cell.v[e]<0)&&(r->cell.i_nodes[q]->cell.k_nodes[e]->cell.prob >= G->thresh)){
                                DELETE = false;
                                break;
                            }
                        }
                    }
                }
            }
            // Looking at Forwards Node
            if(r->cell.k_nodes[q] != P->dead){
                if ((r->cell.k_nodes[q]->cell.v[q]<0)&&(r->cell.k_nodes[q]->cell.prob >= G->thresh)){
                    DELETE = false;
                    break;
                }else{
                    for (int e = 0; e < DIM; e++){
                        if(e!=q){
                            if ((r->cell.k_nodes[q]->cell.i_nodes[e]->cell.v[q]<0)&&(r->cell.k_nodes[q]->cell.i_nodes[e]->cell.v[e]>0)&&(r->cell.k_nodes[q]->cell.i_nodes[e]->cell.prob >= G->thresh)){
                                DELETE = false;
                                break;
                            }

                            if ((r->cell.k_nodes[q]->cell.k_nodes[e]->cell.v[q]<0)&&(r->cell.k_nodes[q]->cell.k_nodes[e]->cell.v[e]<0)&&(r->cell.k_nodes[q]->cell.k_nodes[e]->cell.prob >= G->thresh)){
                                DELETE = false;
                                break;
                            }
                        }
                    }
                }
            }
        }

        if(DELETE){
            r->cell.del_f = 1;
        }
    }
}

void BST_delete_cells(TreeNode* r, BST* P, Grid* G){
    if (r == NULL){
        return;
    }
    BST_delete_cells(r->left, P, G);
    BST_delete_cells(r->right, P, G);

    if (r->cell.del_f == 1){
        P->root = BST_delete_node(P->root, r->key);
    }
}

void BST_get_tree_info(TreeNode* r, BST* P, Grid* G){
    if (r == NULL){
        return;
    }
    BST_get_tree_info(r->left, P, G);
    BST_get_tree_info(r->right, P, G);

    P->max_key = fmax(P->max_key, r->key);
    if(r->cell.prob >= G->thresh){
        P->a_count++; 
    }
    P->tot_count++; 
}

void BST_measurement_update(TreeNode* r, Grid G, Measurement M){
    if (r == NULL){
        return;
    }

    BST_measurement_update(r->left, G, M);
    BST_measurement_update(r->right, G, M);

    double current_state_vec[DIM];
    for(int i = 0; i < DIM; i++){
        current_state_vec[i] = r->cell.state[i]*G.dx[i] - M.mean[i];
    }
    double x = gauss_probability(current_state_vec, (double *)M.cov);

    r->cell.prob *= exp(-x/2);
}

/*==============================================================================
Python Wrapper Functions
==============================================================================*/
BST initialize_gbees(const char* FILE_NAME, Grid* G, double* measure_time){
    //================================ Read in initial discrete measurement info ===============================//
    printf("Reading in initial discrete measurement...\n\n");

    Measurement M = Measurement_create();

    FILE *measurement_file = fopen(FILE_NAME, "r");

    char line[256];
    char *token;
    int count = 0;
    fgets(line, sizeof(line), measurement_file); // skip label line
    fgets(line, sizeof(line), measurement_file);
    token = strtok(line, " ");
    while (token != NULL && count < DIM) {
        M.mean[count++] = strtod(token, NULL);
        token = strtok(NULL, " ");
    }
    count = 0; 

    fgets(line, sizeof(line), measurement_file); // skip blank line
    fgets(line, sizeof(line), measurement_file); // skip label line
    for(int i = 0; i < DIM; i++){
        fgets(line, sizeof(line), measurement_file); 
        token = strtok(line, " ");
        while (token != NULL && count < DIM) {
            M.cov[i][count++] = strtod(token, NULL);
            token = strtok(NULL, " ");
        }
        count = 0; 
    }
    fgets(line, sizeof(line), measurement_file); // skip blank line
    fgets(line, sizeof(line), measurement_file); // skip label line
    fgets(line, sizeof(line), measurement_file); 
    M.T = strtod(line, NULL); *measure_time = M.T; 
    fclose(measurement_file);
    //==========================================================================================================//

    //========================================== Read in user inputs ===========================================//
    printf("Reading in user inputs...\n\n");
 
    memcpy(G->center, M.mean, DIM * sizeof(double)); 
    for(int i = 0; i < DIM; i++){
        G->center[i] = M.mean[i]; 
        G->dx[i] = pow(M.cov[i][i],0.5)/2.0; 
    }
    //==========================================================================================================//
    BST P = BST_create(); 

    printf("Initializing distribution...\n\n");

    BST_initialize_grid(&P, G, M); 
    BST_normalize_tree(P.root, &P, G); 

    return P;
}

void record_data(TreeNode* r, const char* FILE_NAME, Grid G, const double t){
    FILE* file = fopen(FILE_NAME, "w");
    fprintf(file, "%lf\n", t);
    write_file(file, r, G);
    fclose(file);
}

void grow_tree(BST* P, Grid G){
    BST_create_neighbors(P->root, P, G);
    BST_initialize_ik_nodes(P->root, P);
    P->root = BST_balance(P->root);
}

void initialize_vuw(TreeNode* r, Grid G, Traj T){
    if(r == NULL){ 
        return;
    }
    initialize_vuw(r->left, G, T);
    initialize_vuw(r->right, G, T);

    if(r->cell.new_f==0){
        double x[DIM];
        for(int i = 0; i < DIM; i++){
            x[i] = G.dx[i]*r->cell.state[i]+G.center[i];
        }

        double v1 = T.sigma*(x[1]-(x[0]+(G.dx[0]/2.0)));
        double v2 = -(x[1]+(G.dx[1]/2.0))-x[0]*x[2];
        double v3 = -T.b*(x[2]+(G.dx[2]/2.0))+x[0]*x[1]-T.b*T.r;
        double v[DIM] = {v1,v2,v3};
        double u[DIM] = {fmin(v1,0.0),fmin(v2,0.0),fmin(v3,0.0)};
        double w[DIM] = {fmax(v1,0.0),fmax(v2,0.0),fmax(v3,0.0)};
        memcpy(r->cell.v, v, DIM * sizeof(double)); 
        memcpy(r->cell.u, u, DIM * sizeof(double)); 
        memcpy(r->cell.w, w, DIM * sizeof(double)); 
        r->cell.new_f = 1;

        double sum = 0;
        for(int q = 0; q < DIM; q++) {
            sum += fabs(r->cell.v[q]) / G.dx[q];
        }

        r->cell.cfl_dt = 1.0/sum;
    }
}

void update_prob(BST* P, Grid* G, double rt){
    BST_check_cfl_condition(P->root, P); 
    G->dt = fmin(P->cfl_min_dt, rt);
    BST_godunov_method(P, G); 
    BST_update_prob(P->root, G);
    BST_normalize_tree(P->root, P, G); 
    P->cfl_min_dt = INT_MAX;
}

void prune_tree(BST* P, Grid* G){
    BST_mark_cells(P->root, P, G);
    BST_delete_cells(P->root, P, G);
    BST_normalize_tree(P->root, P, G);
    BST_initialize_ik_nodes(P->root, P);
}

void get_tree_info(TreeNode* r, BST* P, Grid* G){
    P->a_count = 0; 
    P->tot_count = 0; 
    P->max_key = 0; 
    BST_get_tree_info(P->root, P, G);
}

void measurement_update(const char* FILE_NAME, BST* P, Grid G, double* measure_time){
    Measurement M = Measurement_create();

    FILE *measurement_file = fopen(FILE_NAME, "r");

    char line[256];
    char *token;
    int count = 0;
    fgets(line, sizeof(line), measurement_file); // skip label line
    fgets(line, sizeof(line), measurement_file);
    token = strtok(line, " ");
    while (token != NULL && count < DIM) {
        M.mean[count++] = strtod(token, NULL);
        token = strtok(NULL, " ");
    }
    count = 0; 

    fgets(line, sizeof(line), measurement_file); // skip blank line
    fgets(line, sizeof(line), measurement_file); // skip label line
    for(int i = 0; i < DIM; i++){
        fgets(line, sizeof(line), measurement_file); 
        token = strtok(line, " ");
        while (token != NULL && count < DIM) {
            M.cov[i][count++] = strtod(token, NULL);
            token = strtok(NULL, " ");
        }
        count = 0; 
    }
    fgets(line, sizeof(line), measurement_file); // skip blank line
    fgets(line, sizeof(line), measurement_file); // skip label line
    fgets(line, sizeof(line), measurement_file); 
    M.T = strtod(line, NULL); *measure_time = M.T; 
    fclose(measurement_file);

    BST_measurement_update(P->root, G, M);
}