#include <stdint.h>
#include <math.h>
#include <stdlib.h>

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

uint64_t state_conversion(const int* state, const int DIM) {
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
