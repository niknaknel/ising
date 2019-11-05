/**
 * Ising Header
 * Author: Annika Nel 19907281
**/

#ifndef PHYS344_ISING_H
#define PHYS344_ISING_H
#endif //PHYS344_ISING_H

#define T_ZERO_NEG 0  // uniform 1s initial state
#define T_ZERO_POS 1  // uniform -1s initial state
#define T_INF 2       // random initial state
#define J 1
#define Tc 2.2        // Critical temperature
#define STATIC_SEED 0
#define TIME_SEED 1
#define TYPE_INT 0
#define TYPE_DOUBLE 1
#define N_SAMPLES 10.0
#define T_MAX 10000
#define THRESHOLD 0.04
#define DIFF_MAX 50
#define TRUE 1
#define FALSE 0

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <limits.h>

#include "ran0.h"
#include "double_ran0.h"

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
    double *Mps;
    long *seed;
} Lattice;

typedef struct {
    double mean;
    double stddev;
} Tuple;

/* Utility Functions */
void write_result(char *file_name, double *M, int t_max);
double *read_result(char *file_name, int t_max);
void run(Lattice *lat, double temp, int t_max);
void initialize(Lattice *lat, double temp, int t_max);
Lattice new(int L, int state, int time_seed);
void sweep(Lattice *lat);
void display_lattice(Lattice *lat);
long *gen_seed(int time_seed);
void free_lattice(Lattice *lat);
void copy_lattice(Lattice *new, Lattice *lat);

/* Output Functions */
void plot_specific_heat(int L);
void phase_diagram(int L);
void autocorrelation(int L);
Tuple sample_magnetization(int L, double temp);
int correlation_time(Lattice *lat, int t_eq, int t_max, int write_to_file);
double chi(int t, int *M, int t_max);
int equilibration_time(Lattice *lat, int L, double temp, int write_to_file);
double specific_heat(Lattice *lat, int t_eq, int t_corr);
void write_spins(int L, double temp, int init_state);
void test();
