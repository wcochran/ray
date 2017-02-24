#ifndef LINEAR_H
#define LINEAR_H

/*
 * solveLinearSystem()
 * Solve linear system of equations using Gaussian row
 * reduction to place system in upper-echelon form, and then
 * backsolve.
 *
 * AX = B
 * A is RxR, X is RxS, B is RxS
 * M = [A|B] (augmented Rx(R+S) matrix)
 * system is assumed to *not* be singular
 */
int solveLinearSystem(int R, int S, double *M[], double *X[]);

#endif /* LINEAR_H */
