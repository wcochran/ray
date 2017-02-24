#ifndef UNIMODALROOT_H
#define UNIMODALROOT_H

double unimodalRoot(void *data, double (*g)(void *data, double t),
		    double a, double b, double ga, double gb, double tmin);

#endif /* UNIMODALROOT_H */
