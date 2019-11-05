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
    t1 = clock();

    elapsed = (double)(t1 - t0) / CLOCKS_PER_SEC;
    printf("time elapsed: %f\n", elapsed);
}

void test()
{
    Lattice lat;
    int i, j;
    int sizes[1] = {100};
    double temps[25] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0};

    for (i = 0; i < 25; i++) {
        printf("T = %.2f: ", temps[i]);
        for (j = 0; j < 1; j++) {
            printf("L%d, ", sizes[j]);
            int t_eq, t_max;
            Lattice lat;

            t_eq = equilibration_time(&lat, sizes[j], temps[i], FALSE);
            t_max = t_eq*2 < T_MAX ? t_eq*2 : T_MAX;

            correlation_time(&lat, t_eq, t_max, TRUE);
        }
        printf("\n");
    }

    free_lattice(&lat);
}

void plot_specific_heat(int L)
{
    FILE * fp;
    int i;
    double T, pop_mean;

    char *fname = malloc(sizeof(char) * 50);
    sprintf(fname, "./out/specific_heat/cvVt_L%d.csv", L);
    fp = fopen(fname,"w");

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

            t_eq = equilibration_time(&lat, L, T, FALSE);
            t_max = t_eq*2 < T_MAX ? t_eq*2 : T_MAX;

            t_corr = correlation_time(&lat, t_eq, t_max, FALSE);
            Cv = specific_heat(&lat, t_eq, t_corr);

            #pragma omp atomic
            pop_mean += Cv / 10.0;
        }
        printf("\n");
        fprintf(fp, "%.1f, %.5f\n", T, pop_mean);
    }

    fclose(fp);
    free(fname);
}

void autocorrelation(int L)
{
    FILE * fp;
    int i;
    double T, pop_mean;

    char *fname = malloc(sizeof(char) * 50);
    sprintf(fname, "./out/autocorr/tauVt_L%d.csv", L);
    fp = fopen(fname,"w");

    fprintf(fp, "T,tau\n");

    for (T = 0.2; T < 5.2; T+=0.2) {
        printf("T = %.2f: ", T);
        pop_mean = 0;

        #pragma omp parallel for shared(pop_mean)
        for (i = 0; i < 10; i++) {
            printf("%d, ", i+1);
            Lattice lat;
            int t_eq = equilibration_time(&lat, L, T, FALSE);
            int t_max = t_eq*2 < T_MAX ? t_eq*2 : T_MAX;
            int tau = correlation_time(&lat, t_eq, t_max, FALSE);

            #pragma omp atomic
            pop_mean += tau / 10.0;
        }
        printf("\n");
        fprintf(fp, "%.1f, %.5f\n", T, pop_mean);
    }

    fclose(fp);
    free(fname);
}


void phase_diagram(int L)
{
    FILE * fp;
    int i;
    double T, pop_mean;
    Tuple stats;

    char *fname = malloc(sizeof(char) * 50);
    sprintf(fname, "./out/phase/mag_sample_L%d.csv", L);
    fp = fopen(fname,"w");

    fprintf(fp, "T,Mps\n");

    for (T = 0.2; T < 5.2; T+=0.2) {
        printf("T = %.2f: ", T);
        pop_mean = 0;

        #pragma omp parallel for shared(pop_mean)
        for (i = 0; i < 10; i++) {
            printf("%d, ", i+1);
            stats = sample_magnetization(L, T);

            #pragma omp atomic
            pop_mean += stats.mean / 10.0;
        }
        printf("\n");
        fprintf(fp, "%.1f, %.5f\n", T, pop_mean);
    }

    free(fname);
    fclose(fp);
}

Tuple sample_magnetization(int L, double temp)
{
    int i, t_max, t_max_nu, t_eq, t_corr;
    double mean, var, *M_sample;
    Lattice lat;
    Tuple stats;

    M_sample = malloc(sizeof(double) * N_SAMPLES);

    // find equilibration time
    t_eq = equilibration_time(&lat, L, temp, FALSE);

    // find correlation time
    t_max = t_eq*2 < T_MAX ? t_eq*2 : T_MAX;
    t_corr = correlation_time(&lat, t_eq, t_max, FALSE);

    t_max_nu = t_eq + 2*N_SAMPLES*t_corr;
    if (t_max_nu > T_MAX) run(&lat, temp, t_max_nu); // rerun if necessary

    // sample magnetization and calculate mean and variance
    for (i = 0; i < N_SAMPLES; i++) {
        M_sample[i] = lat.Mps[t_eq + i*t_corr];
    }

    mean = 0, var = 0;
    for (i = 0; i < N_SAMPLES; i++) mean += M_sample[i]/N_SAMPLES;
    for (i = 0; i < N_SAMPLES; i++) var += pow(M_sample[i] - mean, 2)/(N_SAMPLES-1); // sample variance vs population variance ???
    stats.mean = (temp < Tc) ? fabs(mean) : mean;
    stats.stddev = sqrt(var);

    free_lattice(&lat);

    return stats;
}

int correlation_time(Lattice *lat, int t_eq, int t_max, int write_to_file)
{
    int i, tau, t_shift_max;
    double *X;

    t_shift_max = t_max - t_eq;
    X = malloc(sizeof(double) * t_shift_max);

    X[0] = chi(0, &(lat->M[t_eq]), t_shift_max);

    for (i = 0; i < t_shift_max; i++) {
        X[i] = chi(i, &(lat->M[t_eq]), t_shift_max);
    }

    tau = 0;
    while ((X[tau]/X[0] > 1/M_E) & (tau < t_shift_max)) {
        tau++;
    }

    /* write to file */
    if (write_to_file & (tau > 0)) {
        FILE *fp;
        double xnorm;
        char *fname = malloc(sizeof(char) * 50);

        sprintf(fname, "./out/autocorr/corr_T%.2f_L%d_tau%d.csv", lat->T, lat->L, tau);
        fp = fopen(fname,"w");

        for (i = 0; i < (N_SAMPLES/2)*tau ; i++) {
            xnorm = X[i]/X[0];
            fprintf(fp, "%.5f\n", xnorm);
        }
        fclose(fp);
        free(fname);
    }

    free(X);

    return tau;
}

double chi(int t, int *M, int t_max)
{
    double coeff, chi, sum1, sum2, sum3;
    int tp;

    coeff = 1.0 / (t_max - t);

    // calculate sums
    sum1 = 0, sum2 = 0, sum3 = 0;
    for (tp = 0; tp < t_max - t; tp++) {
        sum1 += M[tp] * M[tp + t];
        sum2 += M[tp];
        sum3 += M[tp + t];
    }

    chi = coeff*sum1 - coeff*sum2*coeff*sum3;

    return chi;
}

int equilibration_time(Lattice *lat, int L, double temp, int write_to_file)
{
    int t_eq, diff_count, coeff;
    double diff;
    Lattice lat1, lat2;

    /* instantiate lattices with different initial states */
    lat1 = new(L, T_ZERO_POS, TIME_SEED);
    lat2 = new(L, T_INF, TIME_SEED);

    /* run simulation up until T_MAX for both lattices */
    run(&lat1, temp, T_MAX);
    run(&lat2, temp, T_MAX);

    diff_count = 0;
    t_eq = 0;

    /** Find equilibrium.
     *  Loop will end once the difference in magnetisation for the two lattices
     *  has been withing the given threshold for DIFF_MAX consecutive time steps.
     **/
    while ((t_eq < T_MAX) & (diff_count < DIFF_MAX)) {

        if (temp < Tc) {
            // compensate for alternate sign equilibrium
            diff = fabs(fabs(lat2.Mps[t_eq]) - fabs(lat1.Mps[t_eq]));

        } else {
            diff = (L / 200.0) * fabs(lat2.Mps[t_eq] - lat1.Mps[t_eq]);
        }

        coeff = abs(temp - Tc) < 0.15 ? 10/sqrt(L) : 1;

        /* check if difference is within specified threshold */
        if (diff < coeff * THRESHOLD) {
            diff_count++;
        } else {
            diff_count = 0;
        }

        t_eq++;
    }

    t_eq -= DIFF_MAX; // move back to start of equilibrium

    /* write to file */
    if (write_to_file) {
        FILE *fp;
        int i, t_max, flip;
        char *fname = malloc(sizeof(char) * 50);

        sprintf(fname, "./out/equil/eq_T%.2f_L%d_teq%d.csv", temp, L, t_eq);
        fp = fopen(fname,"w");

        fprintf(fp, "M1, M2\n");
        t_max = 2*t_eq < T_MAX ? 2*t_eq : T_MAX;
        t_max = t_max > 500 ? t_max : 500;
        flip = lat2.Mps[t_max-1] > 0 ? 1 : -1;

        for (i = 0; i < t_max; i++) {
            fprintf(fp, "%.5f, %.5f\n", lat1.Mps[i], flip*lat2.Mps[i]);
        }
        fclose(fp);
        free(fname);
    }

    /* copy lat1 to given lat argument and free lat2*/
    copy_lattice(lat, &lat1);
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

void write_spins(int L, double temp, int init_state)
{
    FILE *fp;
    Lattice lat;
    char *file_name;
    int i, j;

    file_name = malloc(sizeof(char) * 50);
    sprintf(file_name, "out/spins/spins_L%d_T%.1f.txt", L, temp);

    fp = fopen(file_name,"w");

    lat = new(L, init_state, TIME_SEED);
    initialize(&lat, temp, T_MAX/5);

    /* sweep */
    for (i = 0; i < T_MAX/5; i++) {
        sweep(&lat);

        fprintf(fp, "[");
        for(j = 0; j < lat.N - 1; j++) fprintf(fp, "%d,", lat.s[j]);
        fprintf(fp, "%d]\n", lat.s[lat.N - 1]);
    }

    fclose(fp);
}

void write_result(char *file_name, double *M, int t_max)
{
    FILE *fp;
    int i;

    fp = fopen (file_name,"w");

    for (i = 0; i < t_max; i++) {
        fprintf (fp, "%.5f\n", M[i]);
    }

    fclose(fp);
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

void copy_lattice(Lattice *new, Lattice *lat)
{
    new->L = lat->L;
    new->N = lat->N;
    new->XNN = lat->XNN;
    new->YNN = lat->YNN;
    new->s = lat->s;
    new->T = lat->T;
    new->beta = lat->beta;
    *new->prob = lat->prob[0];
    new->M = lat->M;
    new->E = lat->E;
    new->Mps = lat->Mps;
    new->seed = lat->seed;
}

void free_lattice(Lattice *lat)
{
    free(lat->s);
    free(lat->M);
    free(lat->E);
    free(lat->Mps);
    free(lat->seed);
}
