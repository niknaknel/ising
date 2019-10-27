/**
 * Ising Model
 * Author: Annika Nel 19907281
**/

#include "ising.h"

/* Structs */
typedef struct {
    int L;
    int N;
    int XNN;
    int YNN;
    int *s;
    double T;
    double beta;
    double prob[5];
    int *M;
    int *E;
    long *seed;
} Lattice;

/* Functions */
void run(Lattice *lattice, double temp, int steps);
void initialize(Lattice *lattice, double temp, int steps);
Lattice new(int L, int state, int time_seed);
void sweep(Lattice *lattice);
void display_lattice(Lattice *lattice);
long *gen_seed(int time_seed);
////////

int z = 0;

int main(int argc, char *argv[])
{
    double T;
    int L, init_state, steps;

    /* Handle arguments */
    T = atof(argv[1]);
    init_state = atoi(argv[2]);
    L = 10;
    steps = 100;

    Lattice lat = new(L, init_state, TIME_SEED);
    printf("Mps\n");
    run(&lat, T, steps);

    for (int i = 0; i < steps * lat.N + 1; i++) {
        printf("%.5f\n", lat.M[i]/(double) lat.N);
    }
}

void run(Lattice *lattice, double temp, int steps)
{
    int i, *E_0, *M_0;

    initialize(lattice, temp, steps);

    /* save first position of array */
    E_0 = lattice->E;
    M_0 = lattice->M;

    /* sweep */
    for (i = 0; i < steps; i++) sweep(lattice);

    /* Move pointer back to first position */
    lattice->E = E_0;
    lattice->M = M_0;

}

void initialize(Lattice *lattice, double temp, int steps)
{
    int i, sum, nn, N, XNN, YNN, m;
    int *E, *M;
    double beta;

    N = lattice->N;
    XNN = lattice->XNN;
    YNN = lattice->YNN;
    beta = 1/temp;

    lattice->T = temp;
    lattice->beta = beta;

    /* Precalculate probabilities */
    for (i = 2; i < 5; i += 2) lattice->prob[i] = exp(-2*beta*i);

    /* initialise M and E arrays */
    E = malloc(sizeof(int) * steps * N);
    M = malloc(sizeof(int) * steps * N);

    /* Calculate initial energy */
    sum = 0;
    for (i = 0; i < N; i++) {
        if ((nn = i+XNN) >= N)   nn -= N;
        sum += lattice->s[nn];
        if ((nn = i+YNN) >= N)   nn -= N;
        sum += lattice->s[nn];
    }
    E[0] = -J*sum;

    /* Calculate initial magnetisation */
    m = 0;
    for (i = 0; i < N; i++) m += lattice->s[i];
    M[0] = m;

    lattice->E = E;
    lattice->M = M;
}

Lattice new(int L, int state, int time_seed)
{
    int i, N;
    int *s;

    Lattice lattice;
    N = L*L;
    s = malloc(sizeof(int) * N);

    /* initialise basics */
    lattice.L = L;
    lattice.N = N;
    lattice.XNN = 1;
    lattice.YNN = L;

    /* Generate random seed */
    lattice.seed = gen_seed(time_seed);

    /* Initialize lattice */
    switch(state)
    {
        case T_ZERO_POS:
            for (i = 0; i < N; i++) s[i] = 1;
            break;

        case T_ZERO_NEG:
            for (i = 0; i < N; i++) s[i] = -1;
            break;

        case T_INF:
            for (i = 0; i < N; i++) {
                s[i] = double_ran0(lattice.seed) < 0.5 ? -1 : 1;
            }
            break;

        default:
            for (i = 0; i < N; i++) s[i] = 1;
            break;
    }

    lattice.s = s;

    return lattice;

}

void sweep(Lattice *lattice)
{
    int i, k, nn, sum, delta, N, XNN, YNN, E_old, M_old;
    int *s, *E, *M;
    double *prob;
    long *SEED;

    // get values from lattice
    N = lattice->N;
    XNN = lattice->XNN;
    YNN = lattice->YNN;
    s = lattice->s;
    E = lattice->E;
    M = lattice->M;
    prob = lattice->prob;
    SEED = lattice->seed;

    // run N flip attempts
    for (k = 0; k < N; k++) {
        // save current M & E and increment pointer
        E_old = *E;
        E++;
        M_old = *M;
        M++;
        z++;

        /* Choose a site */
        i = (int) (N * ran0(SEED));

        /* Calculate the sum of the neighbouring spins */
        if ((nn = i+XNN) >= N)   nn -= N;
        sum = s[nn];
        if ((nn = i-XNN) <  0)   nn += N;
        sum += s[nn];
        if ((nn = i+YNN) >= N)   nn -= N;
        sum += s[nn];
        if ((nn = i-YNN) <  0)   nn += N;
        sum += s[nn];

        /* Calculate the change in energy */
        delta = sum * s[i];

        /* Decide whether to flip the spin */
        if (2*delta <= 0 || double_ran0(SEED) <  prob[delta]) {
            // printf("flip! %d\n", i);
            s[i] = -s[i];

            /* Update energy and magnetisation */
            *E = E_old + 2*delta;
            *M = M_old + 2*s[i];

        } else {
            *E = E_old;
            *M = M_old;
        }
    }

    /* move pointers */
    lattice->E = E;
    lattice->M = M;
}

void display_lattice(Lattice *lattice)
{
    int i, j, k;

    for (j = 0; j < lattice->L; j++) {
        for (i = 0; i < lattice->L; i++) {
            k = (j * lattice->L) + i;
            printf("%2d ", lattice->s[k]);
        }
        printf("\n");
    }
    printf("\n");

}

long *gen_seed(int time_seed)
{
    long* r = malloc(sizeof(long));

    if (time_seed) {
        *r = time(0);
    } else {
        *r = 123123;
    }

    return r;
}
