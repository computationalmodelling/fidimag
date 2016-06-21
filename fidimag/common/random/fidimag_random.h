#ifndef __FIDIMAG_RANDOM__
#define __FIDIMAG_RANDOM__



//#include<omp.h>

#define WIDE_PI 3.1415926535897932384626433832795L

//=================================================
//random number, mt19937
typedef struct {
    unsigned int MT[624];
    unsigned int matrix[2];// = { 0, 0x9908b0dfU};
    int index_t;
    int seed;
    
} mt19937_state;

#define	MT19973_RAND_MAX 4294967295u

mt19937_state *create_mt19937_state(void);
void finalize_mt19937_state(mt19937_state *state);

void initial_rng_mt19973(mt19937_state *state, int seed);
//inline unsigned int rand_int(mt19937_state *state);
double random_double(mt19937_state *state);
void gauss_random_vector(mt19937_state *state, double *x, int n);

//================================================
void random_spin_uniform(double *spin, int n);

#endif
