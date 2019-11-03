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

/* Utility Functions */
void write_result(char *file_name, double *M, int t_max);
double *read_result(char *file_name, int t_max);
void run(Lattice *lattice, double temp, int t_max);
void initialize(Lattice *lattice, double temp, int t_max);
Lattice new(int L, int state, int time_seed);
void sweep(Lattice *lattice);
void display_lattice(Lattice *lattice);
long *gen_seed(int time_seed);
void free_lattice(Lattice *lattice);

/* Output Functions */

void phase_diagram();
int correlation_time(double *Mps, int t_eq, int t_max);
double chi(int t, double *Mps, int t_max);
int equilibration_time(int L, double temp);
void test();
////////

int main(int argc, char *argv[])
{
    double elapsed;
    clock_t t0, t1;

    t0 = clock();
    test();
    t1 = clock();

    elapsed = (double)(t1 - t0) / CLOCKS_PER_SEC;
    printf("time elapsed: %f\n", elapsed);
}

void phase_diagram()
{
    FILE * fp;
    int i;
    fp = fopen ("./out/test01.csv","w");

    for(i = 0; i < 25;i++){
        fprintf (fp, "This is line %d\n",i + 1);
    }

    /* close the file*/
    fclose (fp);
}

void test()
{
    int i, t_eq, t_max, tries, *M, L;
    double *Mps, T;

    t_max = 1000;
    t_eq = t_max;
    T = 1.2;
    L = 50;
    tries = 0;
    Mps = malloc(sizeof(double) * t_max);


    Lattice lat = new(L, T_INF, TIME_SEED);
    run(&lat, T, t_max);

    for (i = 0; i < t_max; i++) {
        Mps[i] = lat.M[i] / (double) lat.N;
    }

    while ((t_eq >= t_max) & (tries < 5)) {
        t_eq = equilibration_time(L, T);
        tries++;
    }

    if (tries >= 5) {
        printf("Couldn't find equilibrium!\n");
    } else {
        int t_corr = correlation_time(Mps, t_eq, t_max);
        printf("tau=%d\n", t_corr);
    }

//    free(M);

}

int correlation_time(double *Mps, int t_eq, int t_max)
{
    int i, tau, t_shift_max;
    double *M_eq, *X, xnorm;
    FILE * fp;

    t_shift_max = t_max - t_eq;
    printf("teq: %d = %f, tmax: %d\n", t_eq, Mps[t_eq], t_max);
    M_eq = &Mps[t_eq];

    X = malloc(sizeof(double) * t_shift_max);

    X[0] = chi(0, M_eq, t_shift_max);

    for (i = 0; i < t_shift_max; i++) {
        X[i] = chi(i, M_eq, t_shift_max);
    }

    fp = fopen ("./out/chi.csv","w");
    for (i = 0; i < t_shift_max; i++) {
        xnorm = X[i]/X[0];
        fprintf (fp, "%.5f\n", xnorm);
    }
    fclose (fp);

    tau = 0;
    printf("x0/x0 = %f\n", X[1]/X[0]);
    while ((X[tau]/X[0] > 1/M_E) & (tau < t_shift_max)) {
        tau++;
    }

    // free malloc'd variables
    free(X);
//    free(M_eq);
    free(fp);

    return tau;

}

double chi(int t, double *Mps, int t_max)
{
    double coeff, chi, sum1, sum2, sum3;
    int tp;

    // WORKED WHEN WAS JUST (t_max - t) for total M?????
    coeff = 1 / (double) (t_max - t);

    // calculate first term
    sum1 = 0, sum2 = 0, sum3 = 0;
    for (tp = 0; tp <= t_max - t; tp++) {
        sum1 += Mps[tp] * Mps[tp + t];
        sum2 += Mps[tp];
        sum3 += Mps[tp + t];
    }

    chi = coeff*sum1 - coeff*sum2*coeff*sum3;
    return chi;
}

int equilibration_time(int L, double temp)
{
    int t_eq, t_max, diff_count;
    double diff;
    Lattice lat1, lat2;
    t_max = 1000;

    lat1 = new(L, T_ZERO_POS, TIME_SEED);
    lat2 = new(L, T_INF, TIME_SEED);

    run(&lat1, temp, t_max);
    run(&lat2, temp, t_max);

//    printf("M1,M2\n");

    diff_count = 0;
    t_eq = 0;

    while ((t_eq < t_max) & (diff_count < 5)) {

//        printf("%.5f,%.5f\n", lat1.M[t_eq]/(double) lat1.N, lat2.M[t_eq]/(double) lat2.N);
        if (temp < Tc) {
            // compensate for alternate sign equilibrium
            diff = abs(abs(lat2.M[t_eq]) - abs(lat1.M[t_eq])) / (double) lat1.N;
        } else {
            diff = abs(lat2.M[t_eq] - lat1.M[t_eq]) / (double) lat1.N;
        }

        if (diff < THRESHOLD) {
            diff_count++;
        } else {
            diff_count = 0;
        }

        t_eq++;
    }

    free_lattice(&lat1);
    free_lattice(&lat2);

    return t_eq;

}

void write_result(char *file_name, double *M, int t_max)
{
    FILE *fp;
    int i;

    fp = fopen (file_name,"w");

    for (i = 0; i < t_max; i++) {
        fprintf (fp, "%.5f\n", M[i]);
    }

    fclose (fp);
}

double *read_result(char *file_name, int t_max)
{
    FILE *fp;
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    int i;
    double *M;

    M = malloc(sizeof(double) * t_max);
    fp = fopen (file_name,"r");

    i = 0;
    while ((read = getline(&line, &len, fp)) != -1) {
        M[i] = atof(line);
        i++;
    }

    fclose(fp);
    if (line) free(line);

    return M;
}

void run(Lattice *lattice, double temp, int t_max)
{
    int i, *E_0, *M_0;

    initialize(lattice, temp, t_max);

    /* save first position of array */
    E_0 = lattice->E;
    M_0 = lattice->M;

    /* sweep */
    for (i = 0; i < t_max; i++) sweep(lattice);

    /* Move pointer back to first position */
    lattice->E = E_0;
    lattice->M = M_0;

}

void initialize(Lattice *lattice, double temp, int t_max)
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
    E = malloc(sizeof(int) * t_max);
    M = malloc(sizeof(int) * t_max);

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

    int dE, dM;
    dE = 0;
    dM = 0;

    // get values from lattice
    N = lattice->N;
    XNN = lattice->XNN;
    YNN = lattice->YNN;
    s = lattice->s;
    E = lattice->E;
    M = lattice->M;
    prob = lattice->prob;
    SEED = lattice->seed;

    // save E & M values
    E_old = *E;
    E++;
    M_old = *M;
    M++;

    // run N flip attempts
    for (k = 0; k < N; k++) {
        // save current M & E and increment pointer
//        E_old = *E;
//        E++;
//        M_old = *M;
//        M++;

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
//            *E = E_old + 2*delta;
//            *M = M_old + 2*s[i];
            dE += 2*delta;
            dM += 2*s[i];

        } //else {
//            *E = E_old;
//            *M = M_old;
        //}
    }

    *E = E_old + dE;
    *M = M_old + dM;

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

void free_lattice(Lattice *lattice)
{
    free(lattice->s);
    free(lattice->M);
    free(lattice->E);
    free(lattice->seed);
}
