//
// Created by annika on 2019/10/27.
//

#ifndef PHYS344_ISING_H
#define PHYS344_ISING_H

#endif //PHYS344_ISING_H

#define T_ZERO_NEG 0  // uniform 1s initial state
#define T_ZERO_POS 1  // uniform -1s initial state
#define T_INF 2       // random initial state
#define J 1
#define Tc 2.2        // Critical temperature
#define THRESHOLD 0.1
#define STATIC_SEED 0
#define TIME_SEED 1
#define TYPE_INT 0
#define TYPE_DOUBLE 1

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
void run(Lattice *lattice, double temp, int t_max);
void initialize(Lattice *lattice, double temp, int t_max);
Lattice new(int L, int state, int time_seed);
void sweep(Lattice *lattice);
void display_lattice(Lattice *lattice);
long *gen_seed(int time_seed);
void free_lattice(Lattice *lattice);

/* Output Functions */
Tuple sample_magnetization(int L, double temp, int try);
void phase_diagram();
void autocorrelation();
int corr_hack(double temp);
int correlation_time(int *M, int t_eq, int t_max);
double chi(int t, int *M, int t_max);
int equilibration_time(int L, double temp);
void test();
