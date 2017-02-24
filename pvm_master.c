#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "pvm3.h"
#include "pnmio.h"

#define MAX_SLAVES 18

#define VERBOSE

int main(int argc, char *argv[]) {
  int mytid;
  int n, numSlaves, tids[MAX_SLAVES];
  int status;
  pnm_image *image;
  FILE *f;
  char *base;
  static char fname[200];

  if (argc < 3) {
    fprintf(stderr, "usage: %s <num_slaves> <prog_name> <prog_args...>\n",
	    argv[0]);
    exit(-1);
  }

  if ((numSlaves = atoi(argv[1])) <= 0 || numSlaves > MAX_SLAVES) {
    fprintf(stderr, "bogus number of slaves!\n");
    exit(-1);
  }

  /*
   * Open file that will hold resulting image.
   */
  if ((base = strrchr(argv[2], '/')) == NULL)
    base = argv[2];
  sprintf(fname, "%s.ppm", base);

  if ((f = fopen(fname, "wb")) == NULL) {
    perror(fname);
    exit(-1);
  }
  
  /* 
   * Enroll in PVM.
   */
  if ((mytid = pvm_mytid()) < 0) {
    fprintf(stderr, "Unable to enroll master in PVM!\n");
    exit(-1);
  }

  /*
   * Spawn slaves.
   * Each slave exec's the program specified by argv[2] and we pass
   * our remaining args to it.
   */
  status = pvm_spawn(argv[2], &argv[3], PvmTaskDefault, "", numSlaves, tids);
  if (status < 0) {
    char *err = "unknown";
    switch(status) {
    case PvmBadParam: err = "invalid argument value"; break;
    case PvmNoHost:   err = "bogus host"; break;
    case PvmNoFile:   err = "executable not found"; break;
    case PvmNoMem:    err = "not enough memory on host";
    case PvmSysErr:   err = "pvmd not responding";
    case PvmOutOfRes: err = "out of resources";
    }
    fprintf(stderr, "pvm_spawn() failed [%d]: %s!\n", status, err);
    pvm_exit();
    exit(-1);
  }

  if (status != numSlaves) {
    fprintf(stderr, "status[%d] != numSlaves[%d]!\n", status, numSlaves);
    numSlaves = status;
    fprintf(stderr, "going with %d slaves...\n", numSlaves);
  }

#ifdef VERBOSE
  printf("master task id = %d\n", mytid);
  for (n = 0; n < numSlaves; n++)
    printf("slave #%d task id = %d\n", n, tids[n]);
  fflush(stdout);
#endif  

  /* 
   * Tell each slave who their master is, what their task number is,
   * and how many total slaves there are.
   */
  for (n = 0; n < numSlaves; n++) {
    pvm_initsend(PvmDataRaw);
    pvm_pkint(&mytid, 1, 1);
    pvm_pkint(&n, 1, 1);
    pvm_pkint(&numSlaves, 1, 1);
    if ((status = pvm_send(tids[n], 1)) < 0) {
      int i;
      fprintf(stderr, "pvm_send(slave=%d) error=%d!\n", n, status);
      for (i = 0; i < n; i++)
	pvm_kill(tids[i]);
      pvm_exit(); 
      exit(-1);
    }
#ifdef VERBOSE
    printf("initial data sent to slave %d\n", n);
    fflush(stdout);
#endif  
  }
  
  /*
   * Wait for slaves to respond. Each slave returns a message
   * containing its slave number, the dimensions of the image, 
   * and its portion of computed pixels.
   * When we hear back from the first minion we will allocate
   * an image buffer.
   */
  for (n = 0; n < numSlaves; n++) {
    int r,c, w,h;    
    int slave;

#ifdef VERBOSE
    printf("waiting for %s slave %d...\n", (n == 0) ? "first" : "next", n);
    fflush(stdout);
#endif  

    pvm_recv(-1, 1);
    pvm_upkint(&slave, 1, 1);
    pvm_upkint(&w, 1, 1);
    pvm_upkint(&h, 1, 1);

    /* sanity check */
    if (slave < 0 || slave >= numSlaves ||
	w <= 0 || w > 10000 || h <= 0 || h > 10000) {
      int i;
      fprintf(stderr, "bogus slave number and/or dimensions from slave!\n");
      for (i = 0; i < numSlaves; i++)
	pvm_kill(tids[i]);
      pvm_exit();
      exit(-1);
    }

    /* allocate image */
    if (n == 0) {  
#ifdef VERBOSE
      printf("allocating %dx%d image...\n", w, h);
      fflush(stdout);
#endif  
      if ((image = allocate_ppm_image(w, h)) == NULL) {
	int i;
	fprintf(stderr, "RED ALERT: could not allocate %dx%d image!\n", w,h);
	for (i = 0; i < numSlaves; i++)
	  pvm_kill(tids[i]);
	pvm_exit();
	exit(-1);
      }
    }
      
    /* get pixels from slave */
    for (r = 0; r < h; r++)  
      for (c = 0; c < w; c++) {
	unsigned char color[3];
	if ((r + c) % numSlaves != slave) 
	  continue;      /* not this slave's pixel*/
	pvm_upkbyte(color, 3, 1);
	PPM_PIXEL_R(image,r,c) = color[0];
	PPM_PIXEL_G(image,r,c) = color[1];
	PPM_PIXEL_B(image,r,c) = color[2];
      }

#ifdef VERBOSE
    printf("all pixels received from slave #%d\n", slave);
    fflush(stdout);
#endif
  }

  /*
   * Done with PVM.
   */
  pvm_exit();

  /*
   * Write image.
   */
  write_pnm_image(image, f);
  fclose(f);

  return 0;
}
