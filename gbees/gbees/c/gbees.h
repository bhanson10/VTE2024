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
#include <time.h>

#define TOL 1E-8

/*==============================================================================
                            STRUCTURE DEFINITIONS
==============================================================================*/
typedef struct {
    double *mean;
    double **cov;
    double T;
} Meas;

Meas Meas_create(int dim, const char* M_DIR, const char* M_FILE) {
    char M_PATH[256];
    snprintf(M_PATH, sizeof(M_PATH), "%s%s", M_DIR, M_FILE);

    FILE *m_file = fopen(M_PATH, "r");
    if (m_file == NULL) {
        fprintf(stderr, "Error: could not open file %s\n", M_PATH);
        exit(EXIT_FAILURE);
    }

    Meas M;
    M.mean = malloc(dim * sizeof(double));
    M.cov = malloc(dim * sizeof(double *));
    if (M.mean == NULL || M.cov == NULL) {
        fprintf(stderr, "Error: memory allocation failure during measurement creation\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < dim; i++) {
        M.cov[i] = malloc(dim * sizeof(double));
        if (M.cov[i] == NULL) {
            fprintf(stderr, "Error: memory allocation failure during measurement creation\n");
            exit(EXIT_FAILURE);
        }
    }

    char line[256];
    char *token;
    int count = 0;

    fgets(line, sizeof(line), m_file); // skip label line
    fgets(line, sizeof(line), m_file); // mean vector
    token = strtok(line, " ");
    while (token != NULL && count < dim) { // read mean vector
        M.mean[count++] = strtod(token, NULL);
        token = strtok(NULL, " ");
    }
    count = 0;

    // Read covariance matrix
    fgets(line, sizeof(line), m_file); // skip blank line
    fgets(line, sizeof(line), m_file); // skip label line
    for (int i = 0; i < dim; i++) { // read covariance matrix
        fgets(line, sizeof(line), m_file);
        token = strtok(line, " ");
        while (token != NULL && count < dim) {
            M.cov[i][count++] = strtod(token, NULL);
            token = strtok(NULL, " ");
        }
        count = 0;
    }

    fgets(line, sizeof(line), m_file); // skip blank line
    fgets(line, sizeof(line), m_file); // skip label line
    fgets(line, sizeof(line), m_file); // read T value
    M.T = strtod(line, NULL);

    fclose(m_file);
    return M;
}

void Meas_free(Meas *M, int dim) {
    if (M->mean != NULL) {
        free(M->mean);
    }
    if (M->cov != NULL) {
        for (int i = 0; i < dim; i++) {
            if (M->cov[i] != NULL) {
                free(M->cov[i]);
            }
        }
        free(M->cov);
    }
}

typedef struct TreeNode TreeNode;

typedef struct {
    double prob;
    double *v;
    double *u;
    double *w;
    double *ctu;
    int *state;
    TreeNode **i_nodes;
    TreeNode **k_nodes;
    double dcu;
    double cfl_dt;
    int new_f;
    int ik_f;
    int del_f;
} Cell;

Cell Cell_create(int dim, double prob, double* v, double* u, double* w, double* ctu, int* state, double dcu, int new_f, int ik_f, int del_f) {
    Cell c;
    c.prob = prob;

    // Allocate memory for arrays
    c.v = (double *)malloc(dim * sizeof(double));
    c.u = (double *)malloc(dim * sizeof(double));
    c.w = (double *)malloc(dim * sizeof(double));
    c.ctu = (double *)malloc(dim * sizeof(double));
    c.state = (int *)malloc(dim * sizeof(int));
    c.i_nodes = (TreeNode **)malloc(dim * sizeof(TreeNode *));
    c.k_nodes = (TreeNode **)malloc(dim * sizeof(TreeNode *));

    // Check for memory allocation failure
    if (c.v == NULL || c.u == NULL || c.w == NULL || c.ctu == NULL || c.state == NULL || c.i_nodes == NULL || c.k_nodes == NULL) {
        if (c.v) free(c.v);
        if (c.u) free(c.u);
        if (c.w) free(c.w);
        if (c.ctu) free(c.ctu);
        if (c.state) free(c.state);
        if (c.i_nodes) free(c.i_nodes);
        if (c.k_nodes) free(c.k_nodes);
        fprintf(stderr, "Error: memory allocation failure during cell creation\n");
        exit(EXIT_FAILURE);
    }

    memcpy(c.v, v, dim * sizeof(double));
    memcpy(c.u, u, dim * sizeof(double));
    memcpy(c.w, w, dim * sizeof(double));
    memcpy(c.ctu, ctu, dim * sizeof(double));
    memcpy(c.state, state, dim * sizeof(int));
    for (int i = 0; i < dim; i++) {
        c.i_nodes[i] = NULL;
        c.k_nodes[i] = NULL;
    }
    c.dcu = dcu;
    c.new_f = new_f;
    c.ik_f = ik_f;
    c.del_f = del_f;

    return c;
}

void Cell_free(Cell *c) {
    // Free the allocated memory
    free(c->v);
    free(c->u);
    free(c->w);
    free(c->ctu);
    free(c->state);
    free(c->i_nodes);
    free(c->k_nodes);

    // Set pointers to NULL to avoid dangling pointers
    c->v = NULL;
    c->u = NULL;
    c->w = NULL;
    c->ctu = NULL;
    c->state = NULL;
    c->i_nodes = NULL;
    c->k_nodes = NULL;
}

typedef struct Grid {
    int dim; 
    double thresh;
    double dt;
    double *center;
    double *dx;
} Grid;

Grid Grid_create(int dim, double thresh, double* center, double* dx){
    Grid G; 
    G.dim = dim; 
    G.thresh = thresh; 
    G.dt = 0; 
    G.center = malloc(dim * sizeof(double));
    G.dx = malloc(dim * sizeof(double *));
    if (G.center == NULL || G.dx == NULL) {
        fprintf(stderr, "Error: memory allocation failure during grid creation\n");
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i < dim; i++){
        G.center[i] = center[i]; 
        G.dx[i] = dx[i]; 
    }
    return G;
}

void Grid_free(Grid* G) {
    // Free the allocated memory
    free(G->center);
    free(G->dx);

    // Set pointers to NULL to avoid dangling pointers
    G->center = NULL;
    G->dx = NULL;
}

typedef struct Traj {
    double *coef;
} Traj;

Traj Traj_create(int n, double* coef){
    Traj T; 
    T.coef = malloc(n * sizeof(double));
    if (T.coef == NULL) {
        fprintf(stderr, "Error: memory allocation failure during trajectory creation\n");
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i < n; i++){
        T.coef[i] = coef[i]; 
    }
    return T;
}

typedef struct TreeNode {
    uint64_t key;
    Cell cell;
    TreeNode* left;
    TreeNode* right;
} TreeNode;

TreeNode* TreeNode_create(uint64_t k, Cell c) {
    TreeNode* node = (TreeNode*)malloc(sizeof(TreeNode));
    if (node == NULL) {
        fprintf(stderr, "Error: memory allocation failure during TreeNode creation\n");
        exit(EXIT_FAILURE);
    }
    node->key = k;
    node->cell = c;
    node->left = NULL;
    node->right = NULL;
    return node;
}

void TreeNode_free(TreeNode* r) {
    if (r == NULL) {
        return;
    }
    // Free the left and right subtrees
    TreeNode_free(r->left);
    TreeNode_free(r->right);
    
    // Free the Cell within the node
    Cell_free(&(r->cell));
    
    // Free the node itself
    free(r);
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

void BST_free(BST* P){
    // Free the allocated memory
    free(P->dead);
    free(P->root);

    // Set pointers to NULL to avoid dangling pointers
    P->dead = NULL;
    P->root = NULL;
}

/*==============================================================================
                        NON-MEMBER FUNCTION DEFINITIONS
==============================================================================*/
uint64_t rosenberg_pair(const int* state, int d, int m) { // Recursive Rosenberg pairing
    if (d == 1) {
        return state[0];
    }

    int* new_state = (int*)malloc((d-1) * sizeof(int));
    if (new_state == NULL) {
        fprintf(stderr, "Memory allocation failure during Rosenberg pairing\n");
        exit(EXIT_FAILURE);
    }
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

uint64_t state_conversion(int dim, const int* state) { // Convert from n-dimensional state to unique key via Rosenberg pairing
    int* shift_state = (int*)malloc(dim * sizeof(int));
    if (shift_state == NULL) {
        fprintf(stderr, "Memory allocation failure during state conversion\n");
        exit(EXIT_FAILURE);
    }

    int m;
    uint64_t key;

    // Perform negative shift
    for (int i = 0; i < dim; i++) {
        if (state[i] < 0) {
            shift_state[i] = -2 * state[i] - 1;
        } else {
            shift_state[i] = 2 * state[i];
        }
    }

    // Determine the maximum value in shift_state
    m = shift_state[0];
    for (int i = 1; i < dim; i++) {
        if (shift_state[i] > m) {
            m = shift_state[i];
        }
    }

    // Compute the key using the Rosenberg pairing function
    key = rosenberg_pair(shift_state, dim, m);

    // Free the allocated memory
    free(shift_state);
    
    return key;
}

double mc(double th){ // MC flux limiter
    double min1 = fmin((1 + th)/2.0, 2.0);
    return fmax(0.0, fmin(min1, 2*th)); 
}

double gauss_probability(int dim, double* x, Meas M){ // MC Calculate gaussian probability at state x given mean and covariance
    lapack_int n = dim;
    double* mat = (double *)malloc(dim * dim * sizeof(double));
    if (mat == NULL) {
        fprintf(stderr, "Error: memory allocation failure during gauss probability\n");
        exit(EXIT_FAILURE);
    }
    for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim; j++){
            mat[i * dim + j] = M.cov[i][j];
        }
    }
    lapack_int ipiv[dim];
    lapack_int info;
    info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, mat, n, ipiv);
    info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, mat, n, ipiv);
    double y[dim];
    cblas_dgemv(CblasRowMajor, CblasNoTrans, dim, dim, 1.0, mat, dim, x, 1, 0.0, y, 1);
    free(mat); 
    return exp(-cblas_ddot(dim, x, 1, y, 1)/2.0);
}

/*==============================================================================
                        MEMBER FUNCTION DEFINITIONS
==============================================================================*/
typedef double* (*AdvectionFunc)(double*, double*, Grid*, Traj);

TreeNode* insert_recursive(TreeNode* r, TreeNode* new_node){
    if(r == NULL){
        return new_node;
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
        fprintf(myfile, "%.10e", r->cell.prob);
        for (int i = 0; i < G.dim; i++) {
            fprintf(myfile, " %.10e", G.dx[i] * r->cell.state[i] + G.center[i]);
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

void initialize_vuw(double* (*ADV)(double*, double*, Grid*, Traj), TreeNode* r, BST* P, Grid* G, Traj T){
    if(r == NULL){ 
        return;
    }
    initialize_vuw(ADV, r->left, P, G, T);
    initialize_vuw(ADV, r->right, P, G, T);
    
    if(r->cell.new_f==0){
        double x[G->dim];
        for(int i = 0; i < G->dim; i++){
            x[i] = G->dx[i]*r->cell.state[i]+G->center[i];
        }

        double* v = (*ADV)(x, G->dx, G, T); 

        double sum = 0;
        for(int i = 0; i < G->dim; i++){
            r->cell.v[i] = v[i];
            r->cell.u[i] = fmin(v[i], 0.0); 
            r->cell.w[i] = fmax(v[i], 0.0); 
            sum += fabs(r->cell.v[i]) / G->dx[i];
        }
        r->cell.new_f = 1;
        r->cell.cfl_dt = 1.0/sum;
    }
}

void initialize_ik_nodes(TreeNode* r, BST* P, Grid* G){
    if(r == NULL){
        return;
    }
    initialize_ik_nodes(r->left, P, G); 
    initialize_ik_nodes(r->right, P, G);

    if(r->cell.ik_f == 0){
        for(int i = 0; i < G->dim; i++){
            // Initializing i, k nodes
            int i_state[G->dim]; memcpy(i_state, r->cell.state, G->dim * sizeof(int)); i_state[i] = i_state[i] - 1; uint64_t i_key = state_conversion(G->dim, i_state); 
            int k_state[G->dim]; memcpy(k_state, r->cell.state, G->dim * sizeof(int)); k_state[i] = k_state[i] + 1; uint64_t k_key = state_conversion(G->dim, k_state); 
            TreeNode* i_node = search_recursive(P->root, i_key); TreeNode* k_node = search_recursive(P->root, k_key); 

            if(i_node == NULL){
                i_node = P->dead; r->cell.i_nodes[i] = i_node; 
            }else{
                r->cell.i_nodes[i] = i_node; i_node->cell.k_nodes[i] = r; 
            }

            if(k_node == NULL){
                k_node = P->dead; r->cell.k_nodes[i] = k_node; 
            }else{
                r->cell.k_nodes[i] = k_node; k_node->cell.i_nodes[i] = r; 
            }
        }
        r->cell.ik_f = 1; 
    } 
}

void recursive_loop(int level, int dim, double* dx, int* current_state, double* current_state_vec, double* zeros, Meas M, BST* P, Grid* G){
    if (level == dim) {
        uint64_t key = state_conversion(dim, current_state);
        double x = gauss_probability(dim, current_state_vec, M);
        TreeNode* new_node = TreeNode_create(key, Cell_create(G->dim, x, zeros, zeros, zeros, zeros, current_state, 0, 0, 0, 0));
        P->root = insert_recursive(P->root, new_node);
        return;
    }

    int start = (int) round(-3 * pow(M.cov[level][level], 0.5) / dx[level]);
    int end = (int) round(3 * pow(M.cov[level][level], 0.5) / dx[level]);
    for (int i = start; i <= end; i++) {
        current_state[level] = i;
        current_state_vec[level] = i * dx[level];
        recursive_loop(level + 1, dim, dx, current_state, current_state_vec, zeros, M, P, G);
    }
}

void initialize_grid(double* (*ADV)(double*, double*, Grid*, Traj), BST* P, Grid* G, Meas M, Traj T){
    double zeros[G->dim];
    int max[G->dim]; 
    for(int i = 0; i < G->dim; i++){
        zeros[i] = 0.0; 
        max[i]   = INT_MAX;
    }
    TreeNode* dead_node = TreeNode_create(-1, Cell_create(G->dim, 0, zeros, zeros, zeros, zeros, max, 0, -1, -1, -1)); 
    
    for(int i = 0; i < G->dim; i++){
        dead_node->cell.i_nodes[i] = dead_node; dead_node->cell.k_nodes[i] = dead_node; 
    }
    P->dead = dead_node; 

    int current_state[G->dim]; double current_state_vec[G->dim]; uint64_t key; double x; TreeNode* new_node; 
    recursive_loop(0, G->dim, G->dx, current_state, current_state_vec, zeros, M, P, G);
    initialize_vuw(ADV, P->root, P, G, T);
    initialize_ik_nodes(P->root, P, G);
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
        double current_v[G.dim];
        memcpy(current_v, r->cell.v, G.dim * sizeof(double));
        int current_state[G.dim]; 
        memcpy(current_state, r->cell.state, G.dim * sizeof(int));
        int new_state[G.dim]; uint64_t new_key;
        double zeros[G.dim];
        for(int i = 0; i < G.dim; i++){
            zeros[i] = 0.0; 
        }
        for (int i = 0; i < G.dim; i++){
            memcpy(new_state, current_state, G.dim * sizeof(int));
            // Checking Forward Faces
            if(current_v[i] > 0){
                if(r->cell.k_nodes[i] == P->dead){
                    new_state[i] += 1;
                    new_key = state_conversion(G.dim, new_state);
                    new_node = TreeNode_create(new_key, Cell_create(G.dim, 0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0, 0)); 
                    P->root = insert_recursive(P->root, new_node);

                    // Checking Edges
                    for (int j = 0; j < G.dim; j++){
                        memcpy(new_state, current_state, G.dim * sizeof(int));
                        new_state[i] += 1;
                        if(j != i){
                            if(current_v[j] > 0){
                                new_state[j] += 1;
                                new_key = state_conversion(G.dim, new_state);
                                new_node = TreeNode_create(new_key, Cell_create(G.dim, 0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0, 0)); 
                                P->root = insert_recursive(P->root, new_node);

                            }else if (current_v[j] < 0){
                                new_state[j] -= 1;
                                new_key = state_conversion(G.dim, new_state);
                                new_node = TreeNode_create(new_key, Cell_create(G.dim, 0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0, 0)); 
                                P->root = insert_recursive(P->root, new_node);

                            }
                        }
                    }
                }else{
                    // Checking Edges
                    for (int j = 0; j < G.dim; j++){
                        memcpy(new_state, current_state, G.dim * sizeof(int));
                        new_state[i] += 1;
                        if(j != i){
                            if(current_v[j] > 0){
                                if(r->cell.k_nodes[i]->cell.k_nodes[j] == P->dead){
                                    new_state[j] += 1;
                                    new_key = state_conversion(G.dim, new_state);
                                    new_node = TreeNode_create(new_key, Cell_create(G.dim, 0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0, 0)); 
                                    P->root = insert_recursive(P->root, new_node);

                                }
                            }else if (current_v[j] < 0){
                                if(r->cell.k_nodes[i]->cell.i_nodes[j] == P->dead){
                                    new_state[j] -= 1;
                                    new_key = state_conversion(G.dim, new_state);
                                    new_node = TreeNode_create(new_key, Cell_create(G.dim, 0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0, 0)); 
                                    P->root = insert_recursive(P->root, new_node);

                                }
                            }
                        }
                    }
                }
                // Checking Backward Faces
            }else if(current_v[i] < 0){
                if(r->cell.i_nodes[i] == P->dead){
                    new_state[i] -= 1;
                    new_key = state_conversion(G.dim, new_state);
                    new_node = TreeNode_create(new_key, Cell_create(G.dim, 0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0, 0)); 
                    P->root = insert_recursive(P->root, new_node);

                    // Checking Edges
                    for (int j = 0; j < G.dim; j++){
                        memcpy(new_state, current_state, G.dim * sizeof(int));
                        new_state[i] -= 1;
                        if(j != i){
                            if(current_v[j] > 0){
                                new_state[j] += 1;
                                new_key = state_conversion(G.dim, new_state);
                                new_node = TreeNode_create(new_key, Cell_create(G.dim, 0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0, 0)); 
                                P->root = insert_recursive(P->root, new_node);

                            }else if (current_v[j] < 0){
                                new_state[j] -= 1;
                                new_key = state_conversion(G.dim, new_state);
                                new_node = TreeNode_create(new_key, Cell_create(G.dim, 0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0, 0)); 
                                P->root = insert_recursive(P->root, new_node);

                            }
                        }
                    }
                }else{
                    // Checking Edges
                    for (int j = 0; j < G.dim; j++){
                        memcpy(new_state, current_state, G.dim * sizeof(int));
                        new_state[i] -= 1;
                        if(j != i){
                            if(current_v[j] > 0){
                                if(r->cell.i_nodes[i]->cell.k_nodes[j] == P->dead){
                                    new_state[j] += 1;
                                    new_key = state_conversion(G.dim, new_state);
                                    new_node = TreeNode_create(new_key, Cell_create(G.dim, 0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0, 0)); 
                                    P->root = insert_recursive(P->root, new_node);

                                }
                            }else if (current_v[j] < 0){
                                if(r->cell.i_nodes[i]->cell.i_nodes[j] == P->dead){
                                    new_state[j] -= 1;
                                    new_key = state_conversion(G.dim, new_state);
                                    new_node = TreeNode_create(new_key, Cell_create(G.dim, 0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0, 0)); 
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

void get_dcu(TreeNode* r, BST* P, Grid* G){
    if (r == NULL){
        return;
    }
    get_dcu(r->left, P, G);
    get_dcu(r->right, P, G);

    r->cell.dcu = 0; 
    TreeNode* i_node; TreeNode* k_node; 
    for(int i = 0; i < G->dim; i++){
        r->cell.ctu[i] = 0.0; 
        i_node = r->cell.i_nodes[i]; k_node = r->cell.k_nodes[i];

        double dcu_p = 0; double dcu_m = 0; 

        if(k_node != P->dead){
            dcu_p = r->cell.w[i] * r->cell.prob + r->cell.u[i] * k_node->cell.prob;
        }
        if(i_node != P->dead){
            dcu_m = i_node->cell.w[i]*i_node->cell.prob + i_node->cell.u[i]*r->cell.prob;
        }

        r->cell.dcu -= (G->dt/G->dx[i])*(dcu_p-dcu_m);
    }
}

void update_ctu(TreeNode* r, BST* P, Grid* G){
    if (r == NULL){
        return;
    }
    update_ctu(r->left, P, G);
    update_ctu(r->right, P, G);

    for(int i = 0; i < G->dim; i++){
        TreeNode* i_node = r->cell.i_nodes[i];
        TreeNode* j_node; TreeNode* p_node;
        if(i_node!=P->dead){
            double F = G->dt*(r->cell.prob-i_node->cell.prob)/(2*G->dx[i]);
            for(int j = 0; j < G->dim; j++){
                if (j!=i){
                    j_node = r->cell.i_nodes[j];
                    p_node = i_node->cell.i_nodes[j];

                    r->cell.ctu[j]      -= i_node->cell.w[i] * r->cell.w[j] * F;
                    i_node->cell.ctu[j] -= i_node->cell.u[i] * i_node->cell.w[j] * F;

                    if(j_node!=P->dead){
                        j_node->cell.ctu[j] -= i_node->cell.w[i] * j_node->cell.u[j] * F;
                    }
                    if(p_node!=P->dead){
                        p_node->cell.ctu[j] -= i_node->cell.u[i] * p_node->cell.u[j] * F;
                    }
                }
            }

            // High-Resolution Correction Terms
            double th;
            if (i_node->cell.v[i]>0){
                TreeNode* i_i_node = i_node->cell.i_nodes[i];
                th = (i_node->cell.prob-i_i_node->cell.prob)/(r->cell.prob-i_node->cell.prob);
            }else{
                th = (r->cell.k_nodes[i]->cell.prob-r->cell.prob)/(r->cell.prob-i_node->cell.prob);
            }

            i_node->cell.ctu[i] += fabs(i_node->cell.v[i])*(G->dx[i]/G->dt - fabs(i_node->cell.v[i]))*F*mc(th);
        }
    }
}

void godunov_method(BST* P, Grid* G){
    get_dcu(P->root, P, G);
    update_ctu(P->root, P, G);
}

void update_prob(TreeNode* r, Grid* G){
    if (r == NULL){
        return;
    }

    update_prob(r->left, G);
    update_prob(r->right, G);

    r->cell.prob += r->cell.dcu;
    for(int i = 0; i < G->dim; i++){
        r->cell.prob -= (G->dt/G->dx[i])*(r->cell.ctu[i]-r->cell.i_nodes[i]->cell.ctu[i]);
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

        for(int i = 0; i < G->dim; i++){
            // Looking at Backwards Node

            if(r->cell.i_nodes[i] != P->dead){
                if ((r->cell.i_nodes[i]->cell.v[i]>0)&&(r->cell.i_nodes[i]->cell.prob >= G->thresh)){
                    DELETE = false;
                    break;
                }else{
                    for (int j = 0; j < G->dim; j++){
                        if(j!=i){
                            if ((r->cell.i_nodes[i]->cell.i_nodes[j]->cell.v[i]>0)&&(r->cell.i_nodes[i]->cell.i_nodes[j]->cell.v[j]>0)&&(r->cell.i_nodes[i]->cell.i_nodes[j]->cell.prob >= G->thresh)){
                                DELETE = false;
                                break;
                            }

                            if ((r->cell.i_nodes[i]->cell.k_nodes[j]->cell.v[i]>0)&&(r->cell.i_nodes[i]->cell.k_nodes[j]->cell.v[j]<0)&&(r->cell.i_nodes[i]->cell.k_nodes[j]->cell.prob >= G->thresh)){
                                DELETE = false;
                                break;
                            }
                        }
                    }
                }
            }
            // Looking at Forwards Node
            if(r->cell.k_nodes[i] != P->dead){
                if ((r->cell.k_nodes[i]->cell.v[i]<0)&&(r->cell.k_nodes[i]->cell.prob >= G->thresh)){
                    DELETE = false;
                    break;
                }else{
                    for (int j = 0; j < G->dim; j++){
                        if(j!=i){
                            if ((r->cell.k_nodes[i]->cell.i_nodes[j]->cell.v[i]<0)&&(r->cell.k_nodes[i]->cell.i_nodes[j]->cell.v[j]>0)&&(r->cell.k_nodes[i]->cell.i_nodes[j]->cell.prob >= G->thresh)){
                                DELETE = false;
                                break;
                            }

                            if ((r->cell.k_nodes[i]->cell.k_nodes[j]->cell.v[i]<0)&&(r->cell.k_nodes[i]->cell.k_nodes[j]->cell.v[j]<0)&&(r->cell.k_nodes[i]->cell.k_nodes[j]->cell.prob >= G->thresh)){
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

void grow_tree(double* (*ADV)(double*, double*, Grid*, Traj), BST* P, Grid G, Traj T){
    create_neighbors(P->root, P, G);
    initialize_vuw(ADV, P->root, P, &G, T);
    initialize_ik_nodes(P->root, P, &G);
    P->root = balance(P->root);
}

void prune_tree(BST* P, Grid* G){
    mark_cells(P->root, P, G);
    delete_cells(P->root, P, G);
    normalize_tree(P);
    initialize_ik_nodes(P->root, P, G);
}

void measurement_update_recursive(TreeNode* r, Grid G, Meas M){
    if (r == NULL){
        return;
    }

    measurement_update_recursive(r->left, G, M);
    measurement_update_recursive(r->right, G, M);

    double current_state_vec[G.dim];
    for(int i = 0; i < G.dim; i++){
        current_state_vec[i] = r->cell.state[i]*G.dx[i] - M.mean[i];
    }
    double x = gauss_probability(G.dim, current_state_vec, M);

    r->cell.prob *= exp(-x/2);
}

char* concat_m(const char* str1, const char* str2, int num1) {
    int num1_len = snprintf(NULL, 0, "%d", num1);
    char* str3 = ".txt"; 
    size_t total_len = strlen(str1) + strlen(str2) + num1_len + strlen(str3) + 1; 
    char* result = (char*)malloc(total_len * sizeof(char));
    if (result == NULL) {
        fprintf(stderr, "Error: memory allocation failure during concat_m\n");
        exit(EXIT_FAILURE);
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
        fprintf(stderr, "Error: memory allocation failure during concat_p\n");
        exit(EXIT_FAILURE);
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

void run_gbees(double* (*ADV)(double*, double*, Grid*, Traj), Grid G, Meas M, Traj T, char* P_DIR, char* M_DIR, int NUM_DIST, int NUM_MEAS, int DEL_STEP, int OUTPUT_FREQ, int DIM, bool OUTPUT, bool RECORD, bool MEASURE){
    char* P_PATH;
    double RECORD_TIME = M.T/(NUM_DIST-1);      

    BST P = BST_create(); 

    printf("Initializing distribution...\n\n");

    initialize_grid(ADV, &P, &G, M, T); 
    normalize_tree(&P); 

    printf("Entering time marching...\n\n");

    clock_t start = clock(); 
    double tt = 0;
    for(int nm = 0; nm < NUM_MEAS; nm++){

        get_tree_info(&P, &G); 
        printf("Timestep: %d-0, Program time: %f s, Sim. time: 0", nm, ((double)(clock()-start))/CLOCKS_PER_SEC); 
        printf(" TU, Active/Total Cells: %d/%d, Max key %%: %e\n", P.a_count, P.tot_count, (double)(P.max_key)/(pow(2,64)-1)*100); 
        if(RECORD){P_PATH = concat_p(P_DIR, "/P", nm, "/pdf_", 0); record_data(P.root, P_PATH, G, tt);};

        double mt = 0; int record_count = 1; int step_count = 1; double rt;
        while(fabs(mt - M.T) > TOL) { // time between discrete measurements

            rt = 0;
            while (rt < RECORD_TIME) { // time between PDF recordings
            
                grow_tree(ADV, &P, G, T);
                check_cfl_condition(P.root, &P);
                G.dt = fmin(P.cfl_min_dt, RECORD_TIME - rt);
                rt += G.dt;
                godunov_method(&P, &G);
                update_prob(P.root, &G);
                normalize_tree(&P);

                if (step_count % DEL_STEP == 0) { // deletion procedure
                    prune_tree(&P, &G);
                }

                if ((OUTPUT) && (step_count % OUTPUT_FREQ == 0)) { // print size to terminal
                    get_tree_info(&P, &G);
                    printf("Timestep: %d-%d, Program time: %f s, Sim. time: %f", nm, step_count, ((double)(clock()-start))/CLOCKS_PER_SEC, tt + mt + rt); 
                    printf(" TU, Active/Total Cells: %d/%d, Max key %%: %e\n", P.a_count, P.tot_count, (double)(P.max_key)/(pow(2,64)-1)*100); 
                }

                step_count += 1; P.cfl_min_dt = INT_MAX;
            }
            if (step_count % OUTPUT_FREQ != 0) { // print size to terminal
                get_tree_info(&P, &G);
                printf("Timestep: %d-%d, Program time: %f s, Sim. time: %f", nm, step_count - 1, ((double)(clock()-start))/CLOCKS_PER_SEC, tt + mt + rt); 
                printf(" TU, Active/Total Cells: %d/%d, Max key %%: %e\n", P.a_count, P.tot_count, (double)(P.max_key)/(pow(2,64)-1)*100); 
            }

            if (RECORD) { // record PDF
                printf("\nRECORDING PDF AT: %f TU...\n\n", tt + mt + rt);
                P_PATH = concat_p(P_DIR, "/P", nm, "/pdf_", record_count); 
                record_data(P.root, P_PATH, G, tt + mt + rt);
                record_count += 1;
            }

            mt += rt;
        }

        tt += mt;
        if ((MEASURE) && (nm < NUM_MEAS - 1)) { // peform discrete measurement update

            printf("\nPERFORMING BAYESIAN UPDATE AT: %f TU..\n\n.", tt);

            // freeing previous measurement 
            Meas_free(&M, G.dim);                                                    

            // reading new measurement
            char* M_FILE_NAME = "/measurement";                                             
            char* M_FILE_EXT = ".txt"; 
            int length = snprintf(NULL, 0, ".%s%d%s", M_FILE_NAME, nm+1, M_FILE_EXT) + 1;
            char* M_FILE = (char *)malloc(length);
            if (M_FILE == NULL) {
                fprintf(stderr, "Error: memory allocation failure during discrete measurement update\n");
                exit(EXIT_FAILURE);
            }
            snprintf(M_FILE, length, ".%s%d%s", M_FILE_NAME, nm+1, M_FILE_EXT);
            Meas M = Meas_create(DIM, M_DIR, M_FILE);                                       
            free(M_FILE); 

            // performing discrete update
            measurement_update_recursive(P.root, G, M);                                     
            normalize_tree(&P);
            prune_tree(&P, &G);
        }
    }
    
    return;
}
#endif //GBEES_H