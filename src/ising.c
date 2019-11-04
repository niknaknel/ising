/**
 * Ising Model
 * Author: Annika Nel 19907281
**/

#include "ising.h"

int main(int argc, char *argv[])
{
    double elapsed;
    clock_t t0, t1;

    t0 = clock();
    phase_diagram();
    t1 = clock();

    elapsed = (double)(t1 - t0) / CLOCKS_PER_SEC;
    printf("time elapsed: %f\n", elapsed);
}

void test()
{
    int t_eq, t_max, tries, *M, L;
    double T;

    t_max = 5000;
    t_eq = t_max;
    T = 2.6;
    L = 50;
    tries = 0;

    Lattice lat = new(L, T_ZERO_POS, TIME_SEED);
    run(&lat, T, t_max);

    while ((t_eq >= t_max) & (tries < 5)) {
        t_eq = equilibration_time(L, T);
        tries++;
    }

    if (tries >= 5) {
        printf("Couldn't find equilibrium!\n");
    } else {
        int t_corr = correlation_time(lat.M, t_eq, t_max);
        printf("tau=%d\n", t_corr);
    }

}

void plot_specific_heat()
{
    FILE * fp;
    int i, L;
    double T, pop_mean;

    L = 100;
    fp = fopen ("./out/specific_heat.csv","w");
    fprintf(fp, "T,Cv\n");

    for (T = 0.2; T < 5.2; T+=0.2) {
        printf("T = %.2f: ", T);
        pop_mean = 0;

        #pragma omp parallel for shared(pop_mean)
        for (i = 0; i < 10; i++) {
            printf("%d, ", i+1);

            int t_eq, t_max, t_corr;
            double Cv;
            Lattice lat;

            t_eq = equilibration_time(L, T);
            t_max = t_eq*2;

            lat = new(L, T_ZERO_POS, TIME_SEED);
            run(&lat, T, t_max);

            t_corr = correlation_time(lat.M, t_eq, t_max);

            if (t_eq + N_SAMPLES*t_corr > t_max) run(&lat, T, t_eq + N_SAMPLES*t_corr);

            Cv = specific_heat(&lat, t_eq, t_corr);

            #pragma omp atomic
            pop_mean += Cv / 10.0;
        }
        printf("\n");
        fprintf(fp, "%.1f, %.5f\n", T, pop_mean);
    }

    fclose(fp);
}

void autocorrelation()
{
    FILE * fp;
    int i, L;
    double T, pop_mean;

    L = 100;
    fp = fopen ("./out/tau.csv","w");
    fprintf(fp, "T,tau\n");

    for (T = 0.2; T < 5.2; T+=0.2) {
        printf("T = %.2f: ", T);
        pop_mean = 0;

        #pragma omp parallel for shared(pop_mean)
        for (i = 0; i < 10; i++) {
            printf("%d, ", i+1);
            int t_eq = equilibration_time(L, T);
            int t_max = t_eq*2;
            Lattice lat = new(L, T_ZERO_POS, TIME_SEED);
            run(&lat, T, t_max);
            int tau = correlation_time(lat.M, t_eq, t_max);
            #pragma omp atomic
            pop_mean += tau / 10.0;
        }
        printf("\n");
        fprintf(fp, "%.1f, %.5f\n", T, pop_mean);
    }

    fclose(fp);
}


void phase_diagram()
{
    FILE * fp;
    int i, L;
    double T, pop_mean;
    Tuple stats;

    L = 100;
    fp = fopen ("./out/phase/0mag_sample1.csv","w");
    fprintf(fp, "T,Mps\n");

    for (T = 0.2; T < 5.2; T+=0.2) {
        printf("T = %.2f: ", T);
        pop_mean = 0;

        #pragma omp parallel for shared(pop_mean)
        for (i = 0; i < 10; i++) {
            printf("%d, ", i+1);
            stats = sample_magnetization(L, T, i+1);

            #pragma omp atomic
            pop_mean += stats.mean / 10.0;
        }
        printf("\n");
        fprintf(fp, "%.1f, %.5f\n", T, pop_mean);
    }

    fclose(fp);
}

Tuple sample_magnetization(int L, double temp, int try)
{
    int i, t_max, t_max_nu, t_eq, t_corr;
    double mean, var, *M_sample;
    Lattice lat;
    Tuple stats;

    // initialize constants
    M_sample = malloc(sizeof(double) * N_SAMPLES);

    // find equilibration time
    t_eq = equilibration_time(L ,temp);

    // find correlation time
    t_max = (int) 2 * t_eq;
    lat = new(L, T_ZERO_POS, TIME_SEED);
    run(&lat, temp, t_max); // run for correlation
    t_corr = correlation_time(lat.M, t_eq, t_max);

    // rerun simulation with new max time
    t_max_nu = t_eq + 2*N_SAMPLES*t_corr;
    run(&lat, temp, t_max_nu);

    // sample magnetization and calculate mean and variance
    for (i = 0; i < N_SAMPLES; i++) {
        M_sample[i] = lat.M[t_eq + i*t_corr] / (double) lat.N;
    }

    mean = 0, var = 0;
    for (i = 0; i < N_SAMPLES; i++) mean += M_sample[i]/N_SAMPLES;
    for (i = 0; i < N_SAMPLES; i++) var += pow(M_sample[i] - mean, 2)/(N_SAMPLES-1); // sample variance vs population variance ???
    stats.mean = (temp < Tc) ? fabs(mean) : mean;
    stats.stddev = sqrt(var);

    // write result
//    char *fname = malloc(sizeof(char) * 30);
//    sprintf(fname, "out/phase/mag_%.1f_%02d.csv", temp, try);
//    write_result(fname, lat.Mps, t_max_nu);

    free_lattice(&lat);

    return stats;
}

int correlation_time(int *M, int t_eq, int t_max)
{
    int i, tau, t_shift_max;
    double *X, xnorm;
    FILE * fp;

    t_shift_max = t_max - t_eq;
    X = malloc(sizeof(double) * t_shift_max);

    X[0] = chi(0, &M[t_eq], t_shift_max);

    for (i = 0; i < t_shift_max; i++) {
        X[i] = chi(i, &M[t_eq], t_shift_max);
    }

    tau = 0;
    while ((X[tau]/X[0] > 1/M_E) & (tau < t_shift_max)) {
        tau++;
    }

    /* write to file */
    //    fp = fopen ("./out/chi.csv","w");
//    for (i = 0; i < 250 ; i++) {
//        xnorm = X[i]/X[0];
//        fprintf (fp, "%.5f\n", xnorm);
//    }
//    fclose (fp);

    // free malloc'd variables
    free(X);

    return tau;
}

double chi(int t, int *M, int t_max)
{
    double coeff, chi, sum1, sum2, sum3;
    int tp;

    coeff = 1.0 / (t_max - t);

    // calculate first term
    sum1 = 0, sum2 = 0, sum3 = 0;
    for (tp = 0; tp < t_max - t; tp++) {
        sum1 += M[tp] * M[tp + t];
        sum2 += M[tp];
        sum3 += M[tp + t];
    }

    chi = coeff*sum1 - coeff*sum2*coeff*sum3;

    return chi;
}

int equilibration_time(int L, double temp)
{
    int t_eq, t_max, diff_count;
    double diff;
    Lattice lat1, lat2;
    t_max = 10000;

    lat1 = new(L, T_ZERO_POS, TIME_SEED);
    lat2 = new(L, T_INF, TIME_SEED);

    run(&lat1, temp, t_max);
    run(&lat2, temp, t_max);

    diff_count = 0;
    t_eq = 0;

    while ((t_eq < t_max) & (diff_count < 5)) {

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

double specific_heat(Lattice *lat, int t_eq, int t_corr)
{
    int i;
    double mean, var, Cv;

    mean = 0, var = 0;

    /* calculate mean */
    for (i = 0; i < N_SAMPLES; i++) mean += lat->E[t_eq + i*t_corr]/N_SAMPLES;

    /* calculate variance */
    for (i = 0; i < N_SAMPLES; i++) var += pow(lat->E[t_eq + i*t_corr] - mean, 2) / (N_SAMPLES-1); // sample vs population?

    Cv = (lat->beta / lat->T) * var;
    Cv = Cv / (lat->N); // per site

    return Cv;
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

void run(Lattice *lat, double temp, int t_max)
{
    int i, *E_0, *M_0;
    double *Mps;

    initialize(lat, temp, t_max);

    /* save first position of array */
    E_0 = lat->E;
    M_0 = lat->M;

    /* sweep */
    for (i = 0; i < t_max; i++) sweep(lat);

    /* Move pointer back to first position */
    lat->E = E_0;
    lat->M = M_0;

    /* Find M per site and save to struct */
    Mps = malloc(sizeof(double) * t_max);
    for (i = 0; i < t_max; i++) Mps[i] = lat->M[i] / (double) lat->N;

    lat->Mps = Mps;

}

void initialize(Lattice *lat, double temp, int t_max)
{
    int i, sum, nn, N, XNN, YNN, m;
    int *E, *M;
    double beta;

    // free lattice if 2nd + call to run
    if (lat->M != NULL) free(lat->M);
    if (lat->Mps != NULL) free(lat->Mps);
    if (lat->E != NULL) free(lat->E);

    N = lat->N;
    XNN = lat->XNN;
    YNN = lat->YNN;
    beta = 1/temp;

    lat->T = temp;
    lat->beta = beta;

    /* Precalculate probabilities */
    for (i = 2; i < 5; i += 2) lat->prob[i] = exp(-2*beta*i);

    /* initialise M and E arrays */
    E = malloc(sizeof(int) * (t_max+1));
    M = malloc(sizeof(int) * (t_max+1));
    E[t_max] = INT_MIN;
    M[t_max] = INT_MIN;

    /* Calculate initial energy */
    sum = 0;
    for (i = 0; i < N; i++) {
        if ((nn = i+XNN) >= N)   nn -= N;
        sum += lat->s[nn];
        if ((nn = i+YNN) >= N)   nn -= N;
        sum += lat->s[nn];
    }
    E[0] = -J*sum;

    /* Calculate initial magnetisation */
    m = 0;
    for (i = 0; i < N; i++) m += lat->s[i];
    M[0] = m;

    lat->E = E;
    lat->M = M;
}

Lattice new(int L, int state, int time_seed)
{
    int i, N;
    int *s;

    Lattice lat;
    N = L*L;
    s = malloc(sizeof(int) * N);

    /* initialise basics */
    lat.L = L;
    lat.N = N;
    lat.XNN = 1;
    lat.YNN = L;
    lat.M = NULL;
    lat.E = NULL;
    lat.Mps = NULL;

    /* Generate random seed */
    lat.seed = gen_seed(time_seed);

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
                s[i] = double_ran0(lat.seed) < 0.5 ? -1 : 1;
            }
            break;

        default:
            for (i = 0; i < N; i++) s[i] = 1;
            break;
    }

    lat.s = s;

    return lat;
}

void sweep(Lattice *lat)
{
    int i, k, nn, sum, delta, N, XNN, YNN, E_old, M_old;
    int *s, *E, *M;
    double *prob;
    long *SEED;

    int dE, dM;
    dE = 0;
    dM = 0;

    // get values from lattice (for readability)
    N = lat->N;
    XNN = lat->XNN;
    YNN = lat->YNN;
    s = lat->s;
    E = lat->E;
    M = lat->M;
    prob = lat->prob;
    SEED = lat->seed;

    // save E & M values and increment pointer
    E_old = *E;
    E++;
    M_old = *M;
    M++;

    // run N flip attempts
    for (k = 0; k < N; k++) {
        /* Choose a site */
        i = (int) ((N-1) * ran0(SEED));

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
            dE += 2*delta;
            dM += 2*s[i];

        }
    }

    *E = E_old + dE;
    *M = M_old + dM;

    /* move pointers */
    lat->E = E;
    lat->M = M;
}

void display_lattice(Lattice *lat)
{
    int i, j, k;

    for (j = 0; j < lat->L; j++) {
        for (i = 0; i < lat->L; i++) {
            k = (j * lat->L) + i;
            printf("%2d ", lat->s[k]);
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

void free_lattice(Lattice *lat)
{
    free(lat->s);
    free(lat->M);
    free(lat->E);
    free(lat->Mps);
    free(lat->seed);
}
