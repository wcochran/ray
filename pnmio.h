
/* $Author: cs548 $ $Revision: 1.1 $ $Date: 2009/10/09 19:02:19 $ */

#ifndef _PNMIO_H_
#define _PNMIO_H_

#include <stdlib.h>
#include <stdio.h>

/*
 * pnmio.h
 *
 * structure and functions definitions for pnmio library
 *
 * 1/16/95 PJF
 * 1/26/95 PJF updated from suggestions by Radu
 * 2/7/95 PJF new macro for Wayne
 * 10/1/04 WOC pixel data size moved from unsigned char's to int's
 *
 */

typedef enum {PBM_ASCII=1, PGM_ASCII=2, PPM_ASCII=3,
              PBM_RAW=4, PGM_RAW=5, PPM_RAW=6} pnm_type;

#define PBM PBM_RAW
#define PGM PGM_RAW
#define PPM PPM_RAW

typedef struct _pnmi {
 pnm_type type;
 unsigned int nr, nc; /* # of rows and columns */
 unsigned int maxval; /* max value */
 unsigned int *data;  /* usually a byte is enough, but we may need more */
 } pnm_image;

#define PNM_TYPE(x) (x)->type
#define PNM_NR(x) (x)->nr
#define PNM_NC(x) (x)->nc
#define PNM_MAXVAL(x) (x)->maxval
#define PNM_DATA(x) (x)->data
#define PGM_PIXEL(x,i,j) (x)->data[(i)*PNM_NC(x)+(j)]
#define PPM_PIXEL_R(x,i,j) (x)->data[((i)*PNM_NC(x)+j)*3]
#define PPM_PIXEL_G(x,i,j) (x)->data[((i)*PNM_NC(x)+j)*3+1]
#define PPM_PIXEL_B(x,i,j) (x)->data[((i)*PNM_NC(x)+j)*3+2]

/* this macro deprecated -- WOC */
#define PPM_PIXEL_0RGB(x,i,j) ((PPM_PIXEL_R(x,i,j)<<16)*0x00ff0000)|((PPM_PIXEL_G(x,i,j)<<8)&0x0000ff00)|(PPM_PIXEL_B(x,i,j)&0x000000ff)

#define PBM_PIXEL(x,i,j) PGM_PIXEL(x,i,j)

#define IS_PBM(x) ((PNM_TYPE(x)==PBM_ASCII)||(PNM_TYPE(x)==PBM_RAW))
#define IS_PGM(x) ((PNM_TYPE(x)==PGM_ASCII)||(PNM_TYPE(x)==PGM_RAW))
#define IS_PPM(x) ((PNM_TYPE(x)==PPM_ASCII)||(PNM_TYPE(x)==PPM_RAW))

/* this macro (may be) deprecated -- WOC */
#define PNM_PIXEL_SIZE(x) ((IS_PBM(x)||IS_PGM(x)) ? 1 : 3)

pnm_image *read_pnm_image(FILE *fp);

pnm_image *read_pnm_image_from_file(char *filename);

pnm_image *allocate_pbm_image(int nr,int nc);
pnm_image *allocate_pgm_image(int nr,int nc);
pnm_image *allocate_ppm_image(int nr,int nc);

pnm_image *allocate_pnm_image(pnm_type t,int nr, int nc);

void cleanup_pnm_image(pnm_image *im);

void write_pnm_image(pnm_image *im,FILE *fp);
void write_pnm_image_to_file(pnm_image *im,char *filename);

pnm_image *dup_pnm_image(pnm_image *im);

#endif
