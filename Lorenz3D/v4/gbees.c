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

    double del[DIM] = {pow(M.cov[0][0],0.5)/2, pow(M.cov[1][1],0.5)/2, pow(M.cov[2][2],0.5)/2};
    Grid G = Grid_create(2E-5, 0, M.mean, del);     

    Traj T = Traj_create(4.0, 1.0, 48.0); // Lorenz3D trajectory attributes (sigma, beta, r)

    bool OUTPUT = true;     // Write info to terminal
    bool RECORD = true;     // Write PDFs to .txt file
    bool MEASURE = true;    // Take discrete measurement updates
    int output_freq = 1;   // Number of steps per output to terminal
    int del_step = 25;      // Number of steps per deletion procedure
    int num_dist = 6;       // Number of distributions recorded per measurement
    double record_time = M.T/(num_dist-1); // Time between recording PDF
    //==========================================================================================================//

    //================================================= GBEES ==================================================//
    LL P = LL_create(); 

    printf("Initializing distribution...\n\n");

    initialize_grid(&P, &G, M, T); 
    normalize_tree(&P, &G); 

    printf("Entering time marching...\n\n");

    clock_t start = clock(); 
    double tt = 0;
    for(int nm = 0; nm < NM; nm++){

        get_tree_info(&P, &G); 
        printf("Timestep: %d-0, Program time: %f s, Sim. time: 0 TU, Active/Total Cells: %d/%d\n", nm, ((double)(clock()-start))/CLOCKS_PER_SEC, P.a_count, P.tot_count);
        if(RECORD){P_PATH = concat_p(P_DIR, "/P", nm, "/pdf_", 0); record_data(&P, P_PATH, G, tt);};

        double mt = 0; int record_count = 1; int step_count = 1; double rt;
        while(mt < M.T) { // Time between discrete measurements

            rt = 0;
            while (rt < record_time) { // Time between PDF recordings

                grow_tree(&P, G, T);
                check_cfl_condition(&P);
                G.dt = fmin(P.cfl_min_dt, record_time - rt);
                rt += G.dt;
                godunov_method(&P, &G);
                update_prob(&P, &G);
                normalize_tree(&P, &G);

                /*
                if (step_count % del_step == 0) { // Deletion procedure
                    prune_tree(&P, &G);
                }
                */
                
                if ((OUTPUT) && (step_count % output_freq == 0)) { // Print size to terminal
                    get_tree_info(&P, &G);
                    printf("Timestep: %d-%d, Program time: %f s, Sim. time: %f TU, Active/Total Cells: %d/%d\n", nm, step_count, ((double)(clock()-start))/CLOCKS_PER_SEC, tt + mt + rt, P.a_count, P.tot_count);
                }

                step_count += 1; P.cfl_min_dt = INT_MAX; 

                return 0; 
            }
            /*
            if (step_count % output_freq != 0) {
                get_tree_info(&P, &G);
                printf("Timestep: %d-%d, Program time: %f s, Sim. time: %f TU, Active/Total Cells: %d/%d\n", nm, step_count - 1, ((double)(clock()-start))/CLOCKS_PER_SEC, tt + mt + rt, P.a_count, P.tot_count);
            }

            if (RECORD) { // Record PDF
                printf("\nRECORDING PDF AT: %f TU...\n\n", tt + mt + rt);
                P_PATH = concat_p(P_DIR, "/P", nm, "/pdf_", record_count); 
                record_data(P.head, P_PATH, G, tt + mt + rt);
                record_count += 1;
            }

            mt += rt;
            */
        }
        
        /*
        tt += mt;
        if ((MEASURE) && (nm < NM - 1)) { // Perform discrete measurement update
            printf("\nPERFORMING BAYESIAN UPDATE AT: %f TU..\n\n.", tt);
            Measurement Ml = Measurement_create(M_PATH);

            measurement_update_recursive(P.head, G, Ml);
            normalize_tree(&P, &G);
            prune_tree(&P, &G);
        }
        */
    }

    // Py_Finalize();


    return 0;
}
