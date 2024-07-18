#include "gbees.h"
#include <time.h>

int main(){
    //=================================== Read in initial discrete measurement =================================//
    printf("\nReading in initial discrete measurement...\n\n");

    char* P_DIR = "./Data/PDFs"; // Saved PDFs path
    char* P_PATH;

    const int NM = 1; // Number of measurements
    char* M_DIR = "./Data/Measurements"; // Measurement path
    char* M_FILE = "/measurement0.txt"; 
    char M_PATH[100]; strcpy(M_PATH, M_DIR); strcat(M_PATH, M_FILE);
    Measurement M = Measurement_create(M_PATH);
    //==========================================================================================================//

    //========================================== Read in user inputs ===========================================//
    printf("Reading in user inputs...\n\n");

    double del[DIM] = {pow(M.cov[0][0],0.5)/2, pow(M.cov[1][1],0.5)/2, pow(M.cov[2][2],0.5)/2, pow(M.cov[3][3],0.5)/2, pow(M.cov[4][4],0.5)/2, pow(M.cov[5][5],0.5)/2};
    Grid G = Grid_create(1E-7, 0, M.mean, del);     

    Traj T = Traj_create(2.528017528540000E-5); // PCR3BP trajectory attributes (mu)

    bool OUTPUT = true;     // Write info to terminal
    bool RECORD = true;     // Write PDFs to .txt file
    bool MEASURE = true;    // Take discrete measurement updates
    int output_freq = 20;   // Number of steps per output to terminal
    int del_step = 10;      // Number of steps per deletion procedure
    int num_dist = 17;      // Number of distributions recorded per measurement
    double record_time = M.T/(num_dist-1); // Time between recording PDF
    //==========================================================================================================//

    //================================================= GBEES ==================================================//
    BST P = BST_create(); 

    printf("Initializing distribution...\n\n");

    initialize_grid(&P, &G, M, T); 
    normalize_tree(&P); 

    printf("Entering time marching...\n\n");

    clock_t start = clock(); 
    double tt = 0;
    for(int nm = 0; nm < NM; nm++){

        get_tree_info(&P, &G); 
        printf("Timestep: %d-0, Program time: %f s, Sim. time: 0", nm, ((double)(clock()-start))/CLOCKS_PER_SEC); 
        printf(" TU, Active/Total Cells: %d/%d, Max key %%: %f\n", P.a_count, P.tot_count, (double)(P.max_key)/(pow(2,64)-1)*100); 
        if(RECORD){P_PATH = concat_p(P_DIR, "/P", nm, "/pdf_", 0); record_data(P.root, P_PATH, G, tt);};

        double mt = 0; int record_count = 1; int step_count = 1; double rt;
        while(mt < M.T) { // Time between discrete measurements

            rt = 0;
            while (rt < record_time) { // Time between PDF recordings

                grow_tree(&P, G, T);
                check_cfl_condition(P.root, &P);
                G.dt = fmin(P.cfl_min_dt, record_time - rt);
                rt += G.dt;
                godunov_method(&P, &G);
                update_prob(P.root, &G);
                normalize_tree(&P);

                if (step_count % del_step == 0) { // Deletion procedure
                    prune_tree(&P, &G);
                }

                if ((OUTPUT) && (step_count % output_freq == 0)) { // Print size to terminal
                    get_tree_info(&P, &G);
                    printf("Timestep: %d-%d, Program time: %f s, Sim. time: %f", nm, step_count, ((double)(clock()-start))/CLOCKS_PER_SEC, tt + mt + rt); 
                    printf(" TU, Active/Total Cells: %d/%d, Max key %%: %e\n", P.a_count, P.tot_count, (double)(P.max_key)/(pow(2,64)-1)*100); 
                }

                step_count += 1; P.cfl_min_dt = INT_MAX;
            }
            if ((!OUTPUT) || (step_count % output_freq != 0)){
                get_tree_info(&P, &G);
                printf("Timestep: %d-%d, Program time: %f s, Sim. time: %f", nm, step_count - 1, ((double)(clock()-start))/CLOCKS_PER_SEC, tt + mt + rt); 
                printf(" TU, Active/Total Cells: %d/%d, Max key %%: %e\n", P.a_count, P.tot_count, (double)(P.max_key)/(pow(2,64)-1)*100); 
            }

            if (RECORD) { // Record PDF
                printf("\nRECORDING PDF AT: %f TU...\n\n", tt + mt + rt);
                P_PATH = concat_p(P_DIR, "/P", nm, "/pdf_", record_count); 
                record_data(P.root, P_PATH, G, tt + mt + rt);
                record_count += 1;
            }

            mt += rt;
        }

        tt += mt;
        if ((MEASURE) && (nm < NM - 1)) { // Perform discrete measurement update
            printf("\nPERFORMING BAYESIAN UPDATE AT: %f TU..\n\n.", tt);
            Measurement Ml = Measurement_create(M_PATH);

            measurement_update_recursive(P.root, G, Ml);
            normalize_tree(&P);
            prune_tree(&P, &G);
        }
    }

    return 0;
}
