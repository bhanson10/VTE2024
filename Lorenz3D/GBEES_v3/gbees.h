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
#include <Python.h>

#define DIM 3

/*==============================================================================
Non-member Function DEFINITIONS
==============================================================================*/
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
typedef struct Node Node;

typedef struct Cell { 
    double prob;
    double v[DIM];
    double u[DIM];
    double w[DIM];
    double ctu[DIM];
    int state[DIM];
    Node* i_nodes[DIM];
    Node* k_nodes[DIM];
    double dcu;
    double cfl_dt;
    int new_f;
    int del_f;
} Cell;

Cell Cell_create(double prob, double* v, double* u, double* w, double* ctu, int* state, double dcu, int new_f, int del_f) {
    Cell c;
    c.prob = prob;
    for(int i = 0; i < DIM; i++){
        c.v[i] = v[i]; c.u[i] = u[i]; c.w[i] = w[i]; c.ctu[i] = ctu[i]; c.state[i] = state[i]; 
    }
    c.dcu = dcu; 
    c.new_f = new_f; 
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

typedef struct Node {
    Cell cell;
    struct Node* next;
} Node;

Node* Node_create(Cell c) {
    Node* node = (Node*)malloc(sizeof(Node));
    node->cell = c;
    node->next = NULL;
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

typedef struct LL{
    Node* head;
    Node* tail;
    Node* dead;
    int a_count;
    int tot_count;
    double cfl_min_dt;
    PyObject* adv_func;
 } LL;

LL LL_create() {
    LL P;
    P.head = NULL;
    P.tail = NULL; 
    P.dead = NULL; 
    P.a_count = 0;
    P.tot_count = 0;
    P.cfl_min_dt = INT_MAX;

    return P; 
}

void insert(Node* new_node, LL* P){
    if(P->head == NULL){
        P->head = new_node; 
        P->tail = new_node; 
    }else{
        P->tail->next = new_node; 
        P->tail = new_node; 
    }
}

Node* search(LL* P, int* state) {
    Node* r = P->head;
    while (r != NULL) {
        int isEqual = 1; 
        for (int i = 0; i < DIM; i++) {
            if (r->cell.state[i] != state[i]) {
                isEqual = 0; 
                break;
            }
        }
        if (isEqual) {
            return r; 
        }
        r = r->next;
    }
    return NULL; 
}

void get_sum(LL* P, double* prob_sum){
    Node* r = P->head;
    while (r != NULL) {
        *prob_sum += r->cell.prob;
        r = r->next;
    }
    return;
}

void divide_sum(LL* P, double prob_sum){
    Node* r = P->head;
    while (r != NULL) {
        r->cell.prob /= prob_sum;
        r = r->next;
    }
    return;
}

void normalize_tree(LL* P, Grid* G){
    double prob_sum = 0;
    get_sum(P, &prob_sum);
    divide_sum(P, prob_sum);
}

void initialize_vuw(LL* P, Grid* G, Traj T){
    Node* r = P->head;
    while (r != NULL){
        if(r->cell.new_f==0){
            double x[DIM];
            for(int i = 0; i < DIM; i++){
                x[i] = G->dx[i]*r->cell.state[i]+G->center[i];
            }

            PyObject* p_x = PyList_New(DIM);
            PyObject* p_dx = PyList_New(DIM);
            PyObject* p_T = PyList_New(DIM);

            for(int i = 0; i < DIM; i++){
                PyList_SetItem(p_x, i, PyFloat_FromDouble(x[i]));
                PyList_SetItem(p_dx, i, PyFloat_FromDouble(G->dx[i]));
            }
            PyList_SetItem(p_T, 0, PyFloat_FromDouble(T.sigma));
            PyList_SetItem(p_T, 1, PyFloat_FromDouble(T.b));
            PyList_SetItem(p_T, 2, PyFloat_FromDouble(T.r));

            PyObject* pArgs = PyTuple_New(3);
            PyTuple_SetItem(pArgs, 0, p_x);  
            PyTuple_SetItem(pArgs, 1, p_dx);  
            PyTuple_SetItem(pArgs, 2, p_T); 

            PyObject* p_v = PyObject_CallObject(P->adv_func, pArgs); // call adv_func.py

            double v[DIM]; int count = 0; 
            for (Py_ssize_t i = 0; i < PyList_Size(p_v); i++) {
                v[count] = PyFloat_AsDouble(PyList_GetItem(p_v, i));
                count++; 
            }

            double u[DIM] = {fmin(v[0],0.0),fmin(v[1],0.0),fmin(v[2],0.0)};
            double w[DIM] = {fmax(v[0],0.0),fmax(v[1],0.0),fmax(v[2],0.0)};
            memcpy(r->cell.v, v, DIM * sizeof(double)); 
            memcpy(r->cell.u, u, DIM * sizeof(double)); 
            memcpy(r->cell.w, w, DIM * sizeof(double)); 
            r->cell.new_f = 1;

            double sum = 0;
            for(int q = 0; q < DIM; q++) {
                sum += fabs(r->cell.v[q]) / G->dx[q];
            }

            r->cell.cfl_dt = 1.0/sum;
        }
        r = r->next;
    }
    return; 
}

void initialize_ik_nodes(Node* r, LL* P){
    int l_state[DIM]; 
    memcpy(l_state, r->cell.state, DIM * sizeof(int)); 
    for(int q = 0; q < DIM; q++){
        // Initializing i, k nodes
        int i_state[DIM]; memcpy(i_state, l_state, DIM * sizeof(int)); i_state[q] = i_state[q] - 1;
        int k_state[DIM]; memcpy(k_state, l_state, DIM * sizeof(int)); k_state[q] = k_state[q] + 1; 
        Node* i_node = search(P, i_state); Node* k_node = search(P, k_state); 

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
}

void initialize_grid(LL* P, Grid* G, Measurement M, Traj T){
    double zeros[DIM] = {0.0};
    int max[DIM] = {INT_MAX, INT_MAX, INT_MAX};
    Node* dead_node = Node_create(Cell_create(0, zeros, zeros, zeros, zeros, max, 0, -1, -1)); 
    
    for(int i = 0; i < DIM; i++){
        dead_node->cell.i_nodes[i] = dead_node; dead_node->cell.k_nodes[i] = dead_node; 
    }
    P->dead = dead_node; 

    int current_state[DIM]; double current_state_vec[DIM]; uint64_t key; double x; Node* new_node; 
    for (int i = (int) round(-3*M.cov[0][0])/G->dx[0]; i <= (int) round(3*M.cov[0][0])/G->dx[0]; i++){current_state[0] = i; current_state_vec[0] = i*G->dx[0]; 
        for (int j = (int) round(-3*M.cov[1][1])/G->dx[1]; j <= (int) round(3*M.cov[1][1])/G->dx[1]; j++){current_state[1] = j; current_state_vec[1] = j*G->dx[1];
            for (int k = (int) round(-3*M.cov[2][2])/G->dx[2]; k <= (int) round(3*M.cov[2][2])/G->dx[2]; k++){current_state[2] = k; current_state_vec[2] = k*G->dx[2];
                x = gauss_probability(current_state_vec, (double *)M.cov);
                new_node = Node_create(Cell_create(x, zeros, zeros, zeros, zeros, current_state, 0, 0, 0));
                initialize_ik_nodes(new_node, P); 
                insert(new_node, P);
            }
        }
    }

    // Initialize the Python interpreter
    Py_Initialize();
    PyRun_SimpleString("import sys");
    PyRun_SimpleString("sys.path.append(\".\")");
    PyObject* pName = PyUnicode_DecodeFSDefault("advection");
    PyObject* pModule = PyImport_Import(pName);
    P->adv_func = PyObject_GetAttrString(pModule, "advection");

    initialize_vuw(P, G, T);
}

void create_neighbors(LL* P, Grid G){
    Node* r = P->head; 
    while ((r != NULL)||(r->cell.new_f != 0)){

        if (r->cell.prob >= G.thresh){
            Cell c;
            Node* new_node;
            double current_v[DIM];
            memcpy(current_v, r->cell.v, DIM * sizeof(double));
            int current_state[DIM]; 
            memcpy(current_state, r->cell.state, DIM * sizeof(int));
            int new_state[DIM]; 
            double zeros[DIM] = {0.0};
            for (int q = 0; q < DIM; q++){
                memcpy(new_state, current_state, DIM * sizeof(int));
                // Checking Forward Faces
                if(current_v[q] > 0){
                    if(r->cell.k_nodes[q] == P->dead){
                        new_state[q] += 1;
                        new_node = Node_create(Cell_create(0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0)); 
                        initialize_ik_nodes(new_node, P); 
                        insert(new_node, P);

                        // Checking Edges
                        for (int e = 0; e < DIM; e++){
                            memcpy(new_state, current_state, DIM * sizeof(int));
                            new_state[q] += 1;
                            if(e != q){
                                if(current_v[e] > 0){
                                    new_state[e] += 1;
                                    new_node = Node_create(Cell_create(0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0)); 
                                    initialize_ik_nodes(new_node, P); 
                                    insert(new_node, P);

                                }else if (current_v[e] < 0){
                                    new_state[e] -= 1;
                                    new_node = Node_create(Cell_create(0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0)); 
                                    initialize_ik_nodes(new_node, P); 
                                    insert(new_node, P);
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
                                        new_node = Node_create(Cell_create(0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0)); 
                                        initialize_ik_nodes(new_node, P); 
                                        insert(new_node, P);
                                    }
                                }else if (current_v[e] < 0){
                                    if(r->cell.k_nodes[q]->cell.i_nodes[e] == P->dead){
                                        new_state[e] -= 1;
                                        new_node = Node_create(Cell_create(0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0)); 
                                        initialize_ik_nodes(new_node, P); 
                                        insert(new_node, P);
                                    }
                                }
                            }
                        }
                    }
                // Checking Backward Faces
                }else if(current_v[q] < 0){
                    if(r->cell.i_nodes[q] == P->dead){
                        new_state[q] -= 1;
                        new_node = Node_create(Cell_create(0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0)); 
                        initialize_ik_nodes(new_node, P); 
                        insert(new_node, P);

                        // Checking Edges
                        for (int e = 0; e < DIM; e++){
                            memcpy(new_state, current_state, DIM * sizeof(int));
                            new_state[q] -= 1;
                            if(e != q){
                                if(current_v[e] > 0){
                                    new_state[e] += 1;
                                    new_node = Node_create(Cell_create(0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0)); 
                                    initialize_ik_nodes(new_node, P); 
                                    insert(new_node, P);

                                }else if (current_v[e] < 0){
                                    new_state[e] -= 1;
                                    new_node = Node_create(Cell_create(0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0)); 
                                    initialize_ik_nodes(new_node, P); 
                                    insert(new_node, P);
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
                                        new_node = Node_create(Cell_create(0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0)); 
                                        initialize_ik_nodes(new_node, P); 
                                        insert(new_node, P);
                                    }
                                }else if (current_v[e] < 0){
                                    if(r->cell.i_nodes[q]->cell.i_nodes[e] == P->dead){
                                        new_state[e] -= 1;
                                        new_node = Node_create(Cell_create(0, zeros, zeros, zeros, zeros, new_state, 0, 0, 0)); 
                                        initialize_ik_nodes(new_node, P); 
                                        insert(new_node, P);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        r = r->next; 
    }
    return; 
}

void check_cfl_condition(LL* P){
    Node* r = P->head;
    while (r != NULL) {
        P->cfl_min_dt = fmin(P->cfl_min_dt,r->cell.cfl_dt);
        r = r->next;
    }
    return; 
}

/*
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
    get_dcu(P->head, G);
    update_ctu(P->head, P, G);
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

    bool DELETE = true;
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
        P->head = delete_node(P->head, r->key);
    }
}
*/

void get_tree_info(LL* P, Grid* G){
    P->a_count = 0; 
    P->tot_count = 0; 
    Node* r = P->head;
    while (r != NULL) {
        if(r->cell.prob >= G->thresh){
            P->a_count++; 
        }
        P->tot_count++; 
        r = r->next; 
    }
    return;
}

void record_data(LL* P, const char* FILE_NAME, Grid G, const double t){
    FILE* file = fopen(FILE_NAME, "w");
    fprintf(file, "%lf\n", t);
    Node* r = P->head;
    while (r != NULL) {
        if (r->cell.prob >= G.thresh) {
            fprintf(file, "%lf", r->cell.prob);
            for (int i = 0; i < DIM; i++) {
                fprintf(file, " %lf", G.dx[i] * r->cell.state[i] + G.center[i]);
            }
            fprintf(file, "\n");
        }
        r = r->next;
    }
    fclose(file);
    return; 
}

void grow_tree(LL* P, Grid G, Traj T){
    create_neighbors(P, G);
    initialize_vuw(P, &G, T);
}

/*
void prune_tree(BST* P, Grid* G){
    mark_cells(P->head, P, G);
    delete_cells(P->head, P, G);
    normalize_tree(P, G);
    initialize_ik_nodes(P->head, P);
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
*/

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