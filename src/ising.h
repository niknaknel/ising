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
#define THRESHOLD 0.2
#define STATIC_SEED 0
#define TIME_SEED 1
#define TYPE_INT 0
#define TYPE_DOUBLE 1

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#include "ran0.h"
#include "double_ran0.h"
