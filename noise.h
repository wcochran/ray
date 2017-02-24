#ifndef NOISE_H
#define NOISE_H

/*
 * Single octave solid noise.
 */
double noise3(double x, double y, double z);

/*
 * Turbulence:
 *   sum 1/f*noise(f*x)
 * 1 <= f <= freq
 * Choose a power of 2 for freq.
 */
double turbulence3(double x, double y, double z, double freq);

#endif /* NOISE_H */
