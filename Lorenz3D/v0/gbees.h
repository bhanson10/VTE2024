#ifndef GBEES_H
#define GBEES_H

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

Grid Grid_create(double thresh, double dt, double* center, double* dx){
    Grid G; 
    G.thresh = thresh; 
    G.dt = dt; 
    for(int i = 0; i < DIM; i++){
        G.center[i] = center[i]; 
        G.dx[i] = dx[i]; 
    }
    return G;
}

typedef struct Traj {
    double sigma;
    double b;
    double r;
} Traj;

Traj Traj_create(double sigma, double b, double r){
    Traj T; 
    T.sigma = sigma; 
    T.b = b; 
    T.r = r; 
    return T;
}

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

Measurement Measurement_create(const char* M_PATH) {
    FILE *m_file = fopen(M_PATH, "r");
    Measurement M;
    char line[256];
    char *token;
    int count = 0;
    fgets(line, sizeof(line), m_file); // skip label line
    fgets(line, sizeof(line), m_file);
    token = strtok(line, " ");
    while (token != NULL && count < DIM) {
        M.mean[count++] = strtod(token, NULL);
        token = strtok(NULL, " ");
    }
    count = 0; 

    fgets(line, sizeof(line), m_file); // skip blank line
    fgets(line, sizeof(line), m_file); // skip label line
    for(int i = 0; i < DIM; i++){
        fgets(line, sizeof(line), m_file); 
        token = strtok(line, " ");
        while (token != NULL && count < DIM) {
            M.cov[i][count++] = strtod(token, NULL);
            token = strtok(NULL, " ");
        }
        count = 0; 
    }
    fgets(line, sizeof(line), m_file); // skip blank line
    fgets(line, sizeof(line), m_file); // skip label line
    fgets(line, sizeof(line), m_file); 
    M.T = strtod(line, NULL);
    fclose(m_file); 
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

TreeNode* insert_recursive(TreeNode* r, TreeNode* new_node){
    if(r == NULL){
        r = new_node; 
        return r; 
    }

    if (new_node->key < r->key){
        r->left = insert_recursive(r->left, new_node);
    }else if (new_node->key > r->key){
        r->right = insert_recursive(r->right, new_node);
    }else{
        return r;
    }
    return r;
}

TreeNode* search_recursive(TreeNode* r, uint64_t k){ 
    if ((r == NULL)||(r->key == k)){
        return r;
    }else if (k < r->key){
        return search_recursive(r->left, k);
    }else{
        return search_recursive(r->right,k);
    }
}

TreeNode* min_value_node(TreeNode* node){ 
    TreeNode* current = node;
    while(current->left != NULL){
        current = current->left;
    }
    return current;
}

TreeNode* delete_node(TreeNode* r, uint64_t k){ 
    if(r==NULL){
        return NULL;
    }else if(k < r->key){
        r->left = delete_node(r->left, k);
    }else if(k > r->key){
        r->right = delete_node(r->right,k);
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
            r->right = delete_node(r->right, temp->key);
        }
    }
    return r;
}

int get_height(TreeNode* r){
    if (r == NULL){
        return 0;
    }

    return 1 + fmax(get_height(r->left), get_height(r->right));
}

int get_difference(TreeNode* r){
    int l_height = get_height(r->left);
    int r_height = get_height(r->right);
    int b_factor = l_height - r_height;
    return b_factor;
}

TreeNode* rr_rotate(TreeNode* parent){
    TreeNode* t;
    t = parent->right;
    parent->right = t->left;
    t->left = parent;
    return t;
}

TreeNode* ll_rotate(TreeNode* parent){
    TreeNode* t;
    t = parent->left;
    parent->left = t->right;
    t->right = parent;
    return t;
}

TreeNode* lr_rotate(TreeNode* parent){
    TreeNode* t;
    t = parent->left;
    parent->left = rr_rotate(t);
    return ll_rotate(parent);
}

TreeNode* rl_rotate(TreeNode* parent){
    TreeNode* t;
    t = parent->right;
    parent->right = ll_rotate(t);
    return rr_rotate(parent);
}

TreeNode* balance(TreeNode* r){
    int bal_factor = get_difference(r);
    while(abs(bal_factor) > 1){
        if (bal_factor > 1){
            if (get_difference(r->left) > 0){
                r = ll_rotate(r);
            }else{
                r = lr_rotate(r);
            }
        }else if (bal_factor < -1){
            if (get_difference(r->right) > 0){
                r = rl_rotate(r);
            }else{
                r = rr_rotate(r);
            }
        }
        bal_factor = get_difference(r);
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

void get_sum(TreeNode* r, double* prob_sum){
    if (r == NULL){
        return;
    }
    get_sum(r->left, prob_sum);
    get_sum(r->right, prob_sum);

    *prob_sum += r->cell.prob;
}

void divide_sum(TreeNode* r, BST* P, double prob_sum){
    if (r == NULL){
        return;
    }
    divide_sum(r->left, P, prob_sum);
    divide_sum(r->right, P, prob_sum);

    r->cell.prob /= prob_sum;
}

void normalize_tree(BST* P){
    double prob_sum = 0;
    get_sum(P->root, &prob_sum);
    divide_sum(P->root, P, prob_sum);
}

void initialize_vuw(TreeNode* r, BST* P, Grid* G, Traj T){
    if(r == NULL){ 
        return;
    }
    initialize_vuw(r->left, P, G, T);
    initialize_vuw(r->right, P, G, T);
    
    if(r->cell.new_f==0){
        double x[DIM];
        for(int i = 0; i < DIM; i++){
            x[i] = G->dx[i]*r->cell.state[i]+G->center[i];
        }

        double v[DIM] = {T.sigma*(x[1]-(x[0]+(G->dx[0]/2.0))), -(x[1]+(G->dx[1]/2.0))-x[0]*x[2], -T.b*(x[2]+(G->dx[2]/2.0))+x[0]*x[1]-T.b*T.r}; 
        for(int i = 0; i < DIM; i++){
            r->cell.v[i] = v[i];
            r->cell.u[i] = fmin(v[i], 0.0); 
            r->cell.w[i] = fmax(v[i], 0.0); 
        }
        r->cell.new_f = 1;

        double sum = 0;
        for(int q = 0; q < DIM; q++) {
            sum += fabs(r->cell.v[q]) / G->dx[q];
        }

        r->cell.cfl_dt = 1.0/sum;
    }
}

void initialize_ik_nodes(TreeNode* r, BST* P){
    if(r == NULL){
        return;
    }
    initialize_ik_nodes(r->left, P); 
    initialize_ik_nodes(r->right, P);

    if(r->cell.ik_f == 0){
        int l_state[DIM]; 
        memcpy(l_state, r->cell.state, DIM * sizeof(int)); 
        for(int q = 0; q < DIM; q++){
            // Initializing i, k nodes
            int i_state[DIM]; memcpy(i_state, l_state, DIM * sizeof(int)); i_state[q] = i_state[q] - 1; uint64_t i_key = state_conversion(i_state); 
            int k_state[DIM]; memcpy(k_state, l_state, DIM * sizeof(int)); k_state[q] = k_state[q] + 1; uint64_t k_key = state_conversion(k_state); 
            TreeNode* i_node = search_recursive(P->root, i_key); TreeNode* k_node = search_recursive(P->root, k_key); 

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

void initialize_grid(BST* P, Grid* G, Measurement M, Traj T){
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
                P->root = insert_recursive(P->root, new_node);
            }
        }
    }
    initialize_vuw(P->root, P, G, T);
    initialize_ik_nodes(P->root, P);
}


void create_neighbors(TreeNode* r, BST* P, Grid G){
    if (r == NULL){
        return;
    }
    create_neighbors(r->left, P, G);
    create_neighbors(r->right, P, G);

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
                    P->root = insert_recursive(P->root, new_node);

                    // Checking Edges
                    for (int e = 0; e < DIM; e++){
                        memcpy(new_state, current_state, DIM * sizeof(int));
                        new_state[q] += 1;
                        if(e != q){
                            if(current_v[e] > 0){
                                new_state[e] += 1;
                                new_key = state_conversion(new_state);
                                new_node = TreeNode_create(new_key, Cell_create(0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0, 0)); 
                                P->root = insert_recursive(P->root, new_node);

                            }else if (current_v[e] < 0){
                                new_state[e] -= 1;
                                new_key = state_conversion(new_state);
                                new_node = TreeNode_create(new_key, Cell_create(0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0, 0)); 
                                P->root = insert_recursive(P->root, new_node);

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
                                    P->root = insert_recursive(P->root, new_node);

                                }
                            }else if (current_v[e] < 0){
                                if(r->cell.k_nodes[q]->cell.i_nodes[e] == P->dead){
                                    new_state[e] -= 1;
                                    new_key = state_conversion(new_state);
                                    new_node = TreeNode_create(new_key, Cell_create(0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0, 0)); 
                                    P->root = insert_recursive(P->root, new_node);

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
                    P->root = insert_recursive(P->root, new_node);

                    // Checking Edges
                    for (int e = 0; e < DIM; e++){
                        memcpy(new_state, current_state, DIM * sizeof(int));
                        new_state[q] -= 1;
                        if(e != q){
                            if(current_v[e] > 0){
                                new_state[e] += 1;
                                new_key = state_conversion(new_state);
                                new_node = TreeNode_create(new_key, Cell_create(0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0, 0)); 
                                P->root = insert_recursive(P->root, new_node);

                            }else if (current_v[e] < 0){
                                new_state[e] -= 1;
                                new_key = state_conversion(new_state);
                                new_node = TreeNode_create(new_key, Cell_create(0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0, 0)); 
                                P->root = insert_recursive(P->root, new_node);

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
                                    P->root = insert_recursive(P->root, new_node);

                                }
                            }else if (current_v[e] < 0){
                                if(r->cell.i_nodes[q]->cell.i_nodes[e] == P->dead){
                                    new_state[e] -= 1;
                                    new_key = state_conversion(new_state);
                                    new_node = TreeNode_create(new_key, Cell_create(0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0, 0)); 
                                    P->root = insert_recursive(P->root, new_node);

                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void check_cfl_condition(TreeNode* r, BST* P){
    if (r == NULL){
        return;
    }

    check_cfl_condition(r->left, P);
    check_cfl_condition(r->right, P);

    P->cfl_min_dt = fmin(P->cfl_min_dt,r->cell.cfl_dt);
}

void get_dcu(TreeNode* r, Grid* G){
    if (r == NULL){
        return;
    }
    get_dcu(r->left, G);
    get_dcu(r->right, G);

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

void update_ctu(TreeNode* r, BST* P, Grid* G){
    if (r == NULL){
        return;
    }
    update_ctu(r->left, P, G);
    update_ctu(r->right, P, G);

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

void godunov_method(BST* P, Grid* G){
    get_dcu(P->root, G);
    update_ctu(P->root, P, G);
}

void update_prob(TreeNode* r, Grid* G){
    if (r == NULL){
        return;
    }

    update_prob(r->left, G);
    update_prob(r->right, G);

    r->cell.prob += r->cell.dcu;
    for(int q = 0; q < DIM; q++){
        r->cell.prob -= (G->dt/G->dx[q])*(r->cell.ctu[q]-r->cell.i_nodes[q]->cell.ctu[q]);
    }
}


void mark_cells(TreeNode* r, BST* P, Grid* G){
    if (r == NULL){
        return;
    }
    mark_cells(r->left, P, G);
    mark_cells(r->right, P, G);

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

void delete_cells(TreeNode* r, BST* P, Grid* G){
    if (r == NULL){
        return;
    }
    delete_cells(r->left, P, G);
    delete_cells(r->right, P, G);

    if (r->cell.del_f == 1){
        P->root = delete_node(P->root, r->key);
    }
}

void get_tree_info_recursive(TreeNode* r, BST* P, Grid* G){
    if (r == NULL){
        return;
    }
    get_tree_info_recursive(r->left, P, G);
    get_tree_info_recursive(r->right, P, G);

    P->max_key = fmax(P->max_key, r->key);
    if(r->cell.prob >= G->thresh){
        P->a_count++; 
    }
    P->tot_count++; 
}

void get_tree_info(BST* P, Grid* G){
    P->a_count = 0; 
    P->tot_count = 0; 
    P->max_key = 0;
    get_tree_info_recursive(P->root, P, G); 
}

void record_data(TreeNode* r, const char* FILE_NAME, Grid G, const double t){
    FILE* file = fopen(FILE_NAME, "w");
    fprintf(file, "%lf\n", t);
    write_file(file, r, G);
    fclose(file);
}

void grow_tree(BST* P, Grid G, Traj T){
    create_neighbors(P->root, P, G);
    initialize_vuw(P->root, P, &G, T);
    initialize_ik_nodes(P->root, P);
    P->root = balance(P->root);
}

void prune_tree(BST* P, Grid* G){
    mark_cells(P->root, P, G);
    delete_cells(P->root, P, G);
    normalize_tree(P);
    initialize_ik_nodes(P->root, P);
}

void measurement_update_recursive(TreeNode* r, Grid G, Measurement M){
    if (r == NULL){
        return;
    }

    measurement_update_recursive(r->left, G, M);
    measurement_update_recursive(r->right, G, M);

    double current_state_vec[DIM];
    for(int i = 0; i < DIM; i++){
        current_state_vec[i] = r->cell.state[i]*G.dx[i] - M.mean[i];
    }
    double x = gauss_probability(current_state_vec, (double *)M.cov);

    r->cell.prob *= exp(-x/2);
}

char* concat_m(const char* str1, const char* str2, int num1) {
    int num1_len = snprintf(NULL, 0, "%d", num1);
    char* str3 = ".txt"; 
    size_t total_len = strlen(str1) + strlen(str2) + num1_len + strlen(str3) + 1; 
    char* result = (char*)malloc(total_len * sizeof(char));
    if (result == NULL) {
        printf("Memory allocation failed\n");
        return NULL;
    }
    strcpy(result, str1);
    strcat(result, str2);
    char num1_str[num1_len + 1];
    sprintf(num1_str, "%d", num1);
    strcat(result, num1_str);
    strcat(result, str3);

    return result;
}

char* concat_p(const char* str1, const char* str2, int num1, const char* str3, int num2) {
    int num1_len = snprintf(NULL, 0, "%d", num1);
    int num2_len = snprintf(NULL, 0, "%d", num2);
    char* str4 = ".txt"; 
    size_t total_len = strlen(str1) + strlen(str2) + num1_len + strlen(str3) + num2_len + strlen(str4) + 1; 
    char* result = (char*)malloc(total_len * sizeof(char));
    if (result == NULL) {
        printf("Memory allocation failed\n");
        return NULL;
    }
    strcpy(result, str1);
    strcat(result, str2);
    char num1_str[num1_len + 1];
    sprintf(num1_str, "%d", num1);
    strcat(result, num1_str);
    strcat(result, str3);
    char num2_str[num2_len + 1];
    sprintf(num2_str, "%d", num2);
    strcat(result, num2_str);
    strcat(result, str4);

    return result;
}

#endif //GBEES_H