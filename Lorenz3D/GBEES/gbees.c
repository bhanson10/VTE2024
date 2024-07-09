#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <stdio.h>
#include <lapacke.h>
#include <cblas.h>
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

// double mc(double th){
//     double min1 = fmin((1 + th)/2.0, 2.0);
//     return fmax(0.0, fmin(min1, 2*th)); 
// }

// int get_size(TreeNode* r){
//     if (r == NULL){
//         return 0;
//     }

//     return 1 + get_size(r->left) + get_size(r->right);
// }

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

Cell Cell_create(double prob, double* v, double* u, double* w, double* ctu, double dcu, int new_f, int ik_f, int del_f) {
    Cell c;
    c.prob = prob;
    for(int i = 0; i < DIM; i++){
        c.v[i] = v[i]; c.u[i] = u[i]; c.w[i] = w[i]; c.ctu[i] = ctu[i];
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
    double del[DIM];
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
    P.a_count = 0;
    P.a_count = 0;
    P.tot_count = 1;
    P.max_key = -1;
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

// TreeNode* BST_search_recursive(TreeNode* r, uint64_t k){ 
//     if ((r == NULL)||(r->key == k)){
//         return r;
//     }else if (k < r->key){
//         return search_recursive(r->left, k);
//     }else{
//         return search_recursive(r->right,k);
//     }
// }

// TreeNode* min_value_node(TreeNode* node){ 
//     TreeNode* current = node;
//     while(current->left != NULL){
//         current = current->left;
//     }
//     return current;
// }

// TreeNode* BST_delete_node(TreeNode* r, uint64_t k){ 
//     if(r==NULL){
//         return NULL;
//     }else if(k < r->key){
//         r->left = BST_delete_node(r->left, k);
//     }else if(k > r->key){
//         r->right = BST_delete_node(r->right,k);
//     }else{
//         if(r->left == NULL){
//             TreeNode* temp = r->right;
//             free(r);
//             return temp;
//         } else if (r->right == NULL){
//             TreeNode* temp = r->left;
//             free(r);
//             return temp;
//         }else{
//             TreeNode* temp = min_value_node(r->right);
//             r->key = temp->key; r->cell = temp->cell;
//             r->right = BST_delete_node(r->right, temp->key);
//         }
//     }
//     return r;
// }

// int get_height(TreeNode* r){
//     if (r == NULL){
//         return 0;
//     }

//     return 1 + fmax(get_height(r->left), get_height(r->right));
// }

// int get_difference(TreeNode* r){
//     int l_height = get_height(r->left);
//     int r_height = get_height(r->right);
//     int b_factor = l_height - r_height;
//     return b_factor;
// }

// TreeNode* BST_rr_rotate(TreeNode* parent){
//     TreeNode* t;
//     t = parent->right;
//     parent->right = t->left;
//     t->left = parent;
//     return t;
// }

// TreeNode* BST_ll_rotate(TreeNode* parent){
//     TreeNode* t;
//     t = parent->left;
//     parent->left = t->right;
//     t->right = parent;
//     return t;
// }

// TreeNode* BST_lr_rotate(TreeNode* parent){
//     TreeNode* t;
//     t = parent->left;
//     parent->left = BST_rr_rotate(t);
//     return BST_ll_rotate(parent);
// }

// TreeNode* BST_rl_rotate(TreeNode* parent){
//     TreeNode* t;
//     t = parent->right;
//     parent->right = BST_ll_rotate(t);
//     return BST_rr_rotate(parent);
// }

// TreeNode* BST_balance(TreeNode* r){
//     int bal_factor = BST_get_difference(r);
//     while(abs(bal_factor) > 1){
//         if (bal_factor > 1){
//             if (BST_get_difference(r->left) > 0){
//                 r = BST_ll_rotate(r);
//             }else{
//                 r = BST_lr_rotate(r);
//             }
//         }else if (bal_factor < -1){
//             if (BST_get_difference(r->right) > 0){
//                 r = BST_rl_rotate(r);
//             }else{
//                 r = BST_rr_rotate(r);
//             }
//         }
//         bal_factor = BST_get_difference(r);
//     }
//     return r;
// }

// void write_file(FILE* myfile, Grid G, TreeNode* r){
//     if (r == NULL){
//         return;
//     }

//     write_file(myfile, G, r->left);
//     write_file(myfile, G, r->right);

//     if (r->cell.prob >= G.thresh) {
//         fprintf(myfile, "%lf", r->cell.prob);
//         for (int i = 0; i < DIM; i++) {
//             fprintf(myfile, " %lf", G.del[i] * r->cell.state[i] + G.center[i]);
//         }
//         fprintf(myfile, "\n");
//     }
// }

void BST_initialize_grid(Grid G, Traj T, Measurement M, BST* P){
    double zeros[DIM] = {0.0};
    TreeNode* dead_node = TreeNode_create(-1, Cell_create(0, zeros, zeros, zeros, zeros, 0, -1, -1, -1)); 
    
    for(int i = 0; i < DIM; i++){
        dead_node->cell.i_nodes[i] = dead_node; dead_node->cell.k_nodes[i] = dead_node; 
    }
    P->dead = dead_node; 

    int current_state[DIM]; double current_state_vec[DIM]; uint64_t key; TreeNode* new_node; 
    for (int i = (int) round(-3*M.cov[0][0])/G.del[0]; i <= (int) round(3*M.cov[0][0])/G.del[0]; i++){current_state[0] = i; current_state_vec[0] = i*G.del[0]; 
        for (int j = (int) round(-3*M.cov[1][1])/G.del[1]; j <= (int) round(3*M.cov[1][1])/G.del[1]; j++){current_state[1] = j; current_state_vec[1] = j*G.del[1];
            for (int k = (int) round(-3*M.cov[2][2])/G.del[2]; k <= (int) round(3*M.cov[2][2])/G.del[2]; k++){current_state[2] = k; current_state_vec[2] = k*G.del[2];
                key = state_conversion(current_state);
                double x = gauss_probability(current_state_vec, (double *)M.cov);
                new_node = TreeNode_create(key, Cell_create(x, zeros, zeros, zeros, zeros, 0, 0, 0, 0)); 
                P->root = BST_insert_recursive(P->root, new_node);
            }
        }
    }
}

/*==============================================================================
Python Wrapper Functions
==============================================================================*/

// BST initialize_gbees(char FILE_PATH[], const int NM){
void initialize_gbees(const char* FILE_PATH){

    printf("Reading in discrete measurements...\n\n");

    Measurement M = Measurement_create();

    char filepath[256];
    FILE *file;
    char line[256];
    char *token;
    int count = 0;

    // Construct the full file path
    snprintf(filepath, sizeof(filepath), "%s%s", FILE_PATH, "/measurements0.txt");

    FILE *measurement_file = fopen(filepath, "r");

    fgets(line, sizeof(line), measurement_file); // skip label line
    fgets(line, sizeof(line), measurement_file);
    token = strtok(line, " ");
        while (token != NULL && count < DIM) {
            // Convert token to double
            M.mean[count++] = strtod(token, NULL);
            token = strtok(NULL, " ");
        }

    // Close the file
    fclose(measurement_file);

    printf("%f, %f, %f\n", M.mean[0], M.mean[1], M.mean[2]);
    /*
    printf("Reading in discrete measurements...\n\n");

    Measurement M = Measurement_create();

    char filepath[256];
    char line[256];
    FILE *file;
    char *token;
    double values[100]; 
    int count = 0;

    snprintf(filepath, sizeof(filepath), "%s%s", FILE_PATH, "/measurements0.txt");
    FILE* measurement_file = fopen(filepath, "r");
    
    fgets(line, sizeof(line), measurement_file);

    if (fgets(line, sizeof(line), file) != NULL) {
        token = strtok(line, " ");
        while (token != NULL && count < 100) {
            values[count++] = strtod(token, NULL);
            token = strtok(NULL, " ");
        }

        // Print the parsed doubles
        printf("Parsed doubles from %s:\n", filepath);
        for (int i = 0; i < count; i++) {
            printf("%f\n", values[i]);
        }
    }

    fclose(measurement_file);
    */
}

int main(){

    char* FILE_PATH = "./Data";
    initialize_gbees(FILE_PATH);

    return 0;
}