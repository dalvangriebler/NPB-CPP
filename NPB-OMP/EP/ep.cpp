/*--------------------------------------------------------------------

    Information on NAS Parallel Benchmarks is available at:

    http://www.nas.nasa.gov/Software/NPB/

    Authors: P. O. Frederickson
           D. H. Bailey
           A. C. Woo

    CPP and OpenMP version:
            Dalvan Griebler <dalvangriebler@gmail.com>
            Júnior Löff <loffjh@gmail.com>

--------------------------------------------------------------------*/



#include "npbparams.hpp"
#include <iostream>
#include <../common/npb-CPP.hpp>

/* parameters */
#define	MK		16
#define	MM		(M - MK)
#define	NN		(1 << MM)
#define	NK		(1 << MK)
#define	NQ		10
#define EPSILON		1.0e-8
#define	A		1220703125.0
#define	S		271828183.0
#define	TIMERS_ENABLED	FALSE

/* global variables */
/* common /storage/ */
static double x[(2*NK)+1];
#pragma omp threadprivate(x)
static double q[NQ];

/*--------------------------------------------------------------------
      program EMBAR
c-------------------------------------------------------------------*/
/*
c   This is the serial version of the APP Benchmark 1,
c   the "embarassingly parallel" benchmark.
c
c   M is the Log_2 of the number of complex pairs of uniform (0, 1) random
c   numbers.  MK is the Log_2 of the size of each batch of uniform random
c   numbers.  MK can be set for convenience on a given system, since it does
c   not affect the results.
*/
int main(int argc, char **argv) {
    double Mops, t1, sx, sy, tm, an, gc;
    double dum[3] = { 1.0, 1.0, 1.0 };
    int np,i, k, nit, k_offset, j;
    int nthreads = 1;
    boolean verified;
    char size[13+1];	/* character*13 */

    /*
    c   Because the size of the problem is too large to store in a 32-bit
    c   integer for some classes, we put it into a string (for printing).
    c   Have to strip off the decimal point put in there by the floating
    c   point print statement (internal file)
    */

    printf("\n\n NAS Parallel Benchmarks 4.0 OpenMP C++ version"" - EP Benchmark\n");
    printf("\n\n Developed by: Dalvan Griebler <dalvan.griebler@acad.pucrs.br>\n");
    sprintf(size, "%12.0f", pow(2.0, M+1));
    for (j = 13; j >= 1; j--) {
        if (size[j] == '.') size[j] = ' ';
    }
    printf(" Number of random numbers generated: %13s\n", size);

    verified = FALSE;

    /*
    c   Compute the number of "batches" of random number pairs generated
    c   per processor. Adjust if the number of processors does not evenly
    c   divide the total number
    */
    np = NN;

    /*
    c   Call the random number generator functions and initialize
    c   the x-array to reduce the effects of paging on the timings.
    c   Also, call all mathematical functions that are used. Make
    c   sure these initializations cannot be eliminated as dead code.
    */
    vranlc(0, &(dum[0]), dum[1], &(dum[2]));
    dum[0] = randlc(&(dum[1]), dum[2]);
    for (i = 0; i < 2*NK; i++) x[i] = -1.0e99;
    Mops = log(sqrt(fabs(max(1.0, 1.0))));



    timer_clear(1);
    timer_clear(2);
    timer_clear(3);

    timer_start(1);

    vranlc(0, &t1, A, x);

    /*   Compute AN = A ^ (2 * NK) (mod 2^46). */

    t1 = A;

    for ( i = 1; i <= MK+1; i++) {
        an = randlc(&t1, t1);
    }

    an = t1;
    gc = 0.0;
    sx = 0.0;
    sy = 0.0;

    for ( i = 0; i <= NQ - 1; i++) {
        q[i] = 0.0;
    }

    /*
    c   Each instance of this loop may be performed independently. We compute
    c   the k offsets separately to take into account the fact that some nodes
    c   have more numbers to generate than others
    */
    k_offset = -1;

    #pragma omp parallel copyin(x)
    {
        double t1, t2, t3, t4, x1, x2;
        int kk, i, ik, l;
        double qq[NQ];		/* private copy of q[0:NQ-1] */

        for (i = 0; i < NQ; i++) qq[i] = 0.0;

        #pragma omp for reduction(+:sx,sy)
        for (k = 1; k <= np; k++) {
            kk = k_offset + k;
            t1 = S;
            t2 = an;

            /*  Find starting seed t1 for this kk. */

            for (i = 1; i <= 100; i++) {
                ik = kk / 2;
                if (2 * ik != kk) t3 = randlc(&t1, t2);
                if (ik == 0) break;
                t3 = randlc(&t2, t2);
                kk = ik;
            }

            /*      Compute uniform pseudorandom numbers. */

            if (TIMERS_ENABLED == TRUE) timer_start(3);
            vranlc(2*NK, &t1, A, x);
            if (TIMERS_ENABLED == TRUE) timer_stop(3);

            /*
            c       Compute Gaussian deviates by acceptance-rejection method and
            c       tally counts in concentric square annuli.  This loop is not
            c       vectorizable.
            */
            if (TIMERS_ENABLED == TRUE) timer_start(2);

            for ( i = 1; i <= NK; i++) {
                x1 = 2.0 * x[2*i-1] - 1.0;
                x2 = 2.0 * x[2*i] - 1.0;
                t1 = pow2(x1) + pow2(x2);
                if (t1 <= 1.0) {
                    t2 = sqrt(-2.0 * log(t1) / t1);
                    t3 = (x1 * t2);				/* Xi */
                    t4 = (x2 * t2);				/* Yi */
                    l = max(fabs(t3), fabs(t4));
                    qq[l] += 1.0;				/* counts */
                    sx = sx + t3;				/* sum of Xi */
                    sy = sy + t4;				/* sum of Yi */
                }
            }
            if (TIMERS_ENABLED == TRUE) timer_stop(2);
        }
        #pragma omp critical
        {
            for (i = 0; i <= NQ - 1; i++) q[i] += qq[i];
        }
#if defined(_OPENMP)
        #pragma omp master
        nthreads = omp_get_num_threads();
#endif /* _OPENMP */
    } /* end of parallel region */
    for (i = 0; i <= NQ-1; i++) {
        gc = gc + q[i];
    }

    timer_stop(1);
    tm = timer_read(1);



    nit = 0;
    if (M == 24) {
        if((fabs((sx- (-3.247834652034740e3))/-3.247834652034740e3) <= EPSILON) && (fabs((sy- (-6.958407078382297e3))/-6.958407078382297e3) <= EPSILON)) {
            verified = TRUE;
        }
    } else if (M == 25) {
        if ((fabs((sx- (-2.863319731645753e3))/-2.863319731645753e3) <= EPSILON) && (fabs((sy- (-6.320053679109499e3))/-6.320053679109499e3) <= EPSILON)) {
            verified = TRUE;
        }
    } else if (M == 28) {
        //if ((fabs((sx- (-4.295875165629892e3))/sx) <= EPSILON) && (fabs((sy- (-1.580732573678431e4))/sy) <= EPSILON)) {
        if ((fabs((sx- (-4.295875165629892e3))/-4.295875165629892e3) <= EPSILON) && (fabs((sy- (-1.580732573678431e4))/-1.580732573678431e4) <= EPSILON)) {
            verified = TRUE;
        }
    } else if (M == 30) {
        if ((fabs((sx- (4.033815542441498e4))/4.033815542441498e4) <= EPSILON) && (fabs((sy- (-2.660669192809235e4))/-2.660669192809235e4) <= EPSILON)) {
            verified = TRUE;
        }
    } else if (M == 32) {
        if ((fabs((sx- (4.764367927995374e4))/4.764367927995374e4) <= EPSILON) && (fabs((sy- (-8.084072988043731e4))/-8.084072988043731e4) <= EPSILON)) {
            verified = TRUE;
        }
    } else if (M == 36) {
        if ((fabs((sx- (1.982481200946593e5))/1.982481200946593e5) <= EPSILON) && (fabs((sy- (-1.020596636361769e5))/-1.020596636361769e5) <= EPSILON)) {
            verified = TRUE;
        }
    } else if (M == 40) {
        if ((fabs((sx- (-5.319717441530e5))/-5.319717441530e5) <= EPSILON) && (fabs((sy- (-3.688834557731e5))/-3.688834557731e5) <= EPSILON)) {
            verified = TRUE;
        }
    }

    Mops = pow(2.0, M+1)/tm/1000000.0;

    printf("EP Benchmark Results: \n" "CPU Time = %10.4f\n" "N = 2^%5d\n" "No. Gaussian Pairs = %15.0f\n"
           "Sums = %25.15e %25.15e\n" "Counts:\n", tm, M, gc, sx, sy);
    for (i = 0; i  <= NQ-1; i++) {
        printf("%3d %15.0f\n", i, q[i]);
    }

    c_print_results((char*)"EP", CLASS, M+1, 0, 0, nit, nthreads, tm, Mops, (char*)"Random numbers generated",
                    verified, (char*)NPBVERSION, (char*)COMPILETIME, (char*)CS1, (char*)CS2, (char*)CS3, (char*)CS4, (char*)CS5, (char*)CS6, (char*)CS7);

    if (TIMERS_ENABLED == TRUE) {
        printf("Total time:     %f", timer_read(1));
        printf("Gaussian pairs: %f", timer_read(2));
        printf("Random numbers: %f", timer_read(3));
    }
    return 0;
}
