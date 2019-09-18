
/* $Author: cs548 $ $Revision: 1.1 $ $Date: 2009/10/09 19:02:19 $ */

/*
 * pnmio.c
 *
 * library support routines for pbm, pgm, ppm files as defined in
 * Jef Poskanzer's pbmplus routines.
 *
 * 1/16/95 PJF new
 * 1/22/95 added allocate_pXm_image (x = B, G, P), write_pnm_image(),
 *         write_pnm_image_to_file()
 * 1/26/95 fixed a bug in the code which writes PBM files in raw format.
 * 1/26/95 added code to support output of ASCII formats, plus
 *         cleanup_pnm_image(), allocate_pnm_image() per Radu
 * 1/31/95 fixed a bug in the code which *reads* PBM files in raw format.
 * 1/7/96 added clear_pnm_image(), dup_pnm_image()
 * 10/1/04 WOC pixel data size moved from unsigned char's to int's
 *
 * I'm sure there are bugs and bad assumptions in here.
 *
 */

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#include "pnmio.h"


static unsigned int getint(FILE *fp)
{
  char c;
  char buf[256],*ptr=buf;

  while(!isdigit((c=getc(fp)))) {
   if (c=='#') while(getc(fp)!='\n') ;
   };
  *ptr++=c;
  while (isdigit((c=getc(fp)))) {
   *ptr++=c;
   };

  /* this routine gobbles the character after the number too. */
  /* make sure if c=='#' we skip over the comment */
  if (c=='#') while (getc(fp)!='\n') ;
  *ptr=0;

  return atoi(buf);
}

#define CHECK(im) {\
 if (PNM_DATA(im)==NULL) {\
  perror("allocate_pnm_image failed"); \
  free(im); \
  return NULL; \
  } \
 }

pnm_image *allocate_pnm_image(pnm_type t, int nr, int nc)
{
  pnm_image *im=(pnm_image *)malloc(sizeof(pnm_image));
  unsigned int size = ((t==PPM_RAW)||(t==PPM_ASCII)) ? 3 : 1 ;

  if (!im) {
   perror("allocate_pnm_image failed");
   return NULL;
   };
  PNM_NR(im)=nr;
  PNM_NC(im)=nc;
  PNM_DATA(im) = (unsigned int *)malloc(size*nr*nc*sizeof(unsigned int));
  CHECK(im);
  PNM_TYPE(im) = t;
  if ((t != PBM_RAW) && (t != PBM_ASCII)) PNM_MAXVAL(im) = 255;
  return im;
}

pnm_image *allocate_pbm_image(int nr,int nc)
{
  return allocate_pnm_image(PBM_RAW,nr,nc);
}

pnm_image *allocate_pgm_image(int nr,int nc)
{
  return allocate_pnm_image(PGM_RAW,nr,nc);
}

pnm_image *allocate_ppm_image(int nr,int nc)
{
  return allocate_pnm_image(PPM_RAW,nr,nc);
}

pnm_image *read_pnm_image_from_file(char *filename)
{
 int l=(int)strlen(filename);
 FILE *fp;
 char cmd[256];
 pnm_image *im;

 /* .Z ? */
 if ((filename[l-1]=='Z') && (filename[l-2] == '.')) 
  sprintf(cmd,"zcat %s",filename);
 else if ((filename[l-1]=='z') && (filename[l-2]=='g') && (filename[l-3]=='.'))
  sprintf(cmd,"gunzip -c %s",filename);
 else
  sprintf(cmd,"cat %s",filename);

 fp=popen(cmd,"r");
 if (!fp) {
  perror("popen");
  return NULL;
  };

 im=read_pnm_image(fp);

 pclose(fp);
 return im;
}


pnm_image *read_pnm_image(FILE *fp)
{
  unsigned char p,pnmfiletype,*buf,mask;
  unsigned int *ptr;
  int nr,nc,i,j;
  pnm_image *im=(pnm_image *)malloc(sizeof(pnm_image));

  if ((p=getc(fp)) != 'P') {
   fprintf(stderr,"Error: input file is not a PNM image (lacks 'P' as first character)\n");
   return NULL;
   };

  pnmfiletype=getc(fp);

  nc=PNM_NC(im) = getint(fp);
  nr=PNM_NR(im) = getint(fp);

  switch(pnmfiletype) {
   case '1': /* PBM file, ASCII format */
    PNM_TYPE(im)=PBM;
    ptr=PNM_DATA(im)=(unsigned int *)malloc(sizeof(unsigned int)*nr*nc);
    /* we can't use getint here, because I have seen PBM files which string
       the bits together.  I have also seen the bits separated by whitespace. */
    for(i=0;i<nr;i++) {
     for(j=0;j<nc;j++) {
      while (!isdigit(p=getc(fp))) ;
      if (feof(fp)) {
       fprintf(stderr,"end-of-file reached reading ASCII PBM file.\n");
       return NULL;
       };
      PBM_PIXEL(im,i,j) = p-'0';
      };
     };
    break;

   case '4': /* PBM file, binary format */
    PNM_TYPE(im)=PBM;
    buf=(unsigned char *)malloc((sizeof(unsigned char)*nc+7)/8);
    PNM_DATA(im)=(unsigned int *)malloc(sizeof(unsigned int)*nr*nc);
    for(i=0;i<nr;i++) {
     if (fread(buf,sizeof(char),(nc+7)/8,fp) != (nc+7)/8) {
      perror("fread");
      return NULL;
      };
     mask=0x80;
     ptr=(unsigned int *)buf;  /* WOC -- ? */
     for(j=0;j<nc;j++) {
      if (*ptr&mask) PBM_PIXEL(im,i,j)=1;
      mask >>= 1;
      if (!mask) { ptr++; mask=0x80; };
      };
     };
    free(buf);
    break;
    
   case '2': /* PGM file, ASCII format */
    PNM_TYPE(im)=PGM;
    PNM_MAXVAL(im)=getint(fp);
    PNM_DATA(im)=(unsigned int *)malloc(nr*nc*sizeof(unsigned int));
    for(i=0;i<nr;i++) for(j=0;j<nc;j++) PGM_PIXEL(im,i,j)=getint(fp);
    break;
 
   case '5': /* PGM file, raw format */
    PNM_TYPE(im)=PGM;
    PNM_MAXVAL(im)=getint(fp);
    PNM_DATA(im)=(unsigned int *)malloc(nr*nc*sizeof(unsigned int));
    buf=(unsigned char *)malloc(nr*nc*sizeof(unsigned char));
    i=fread(buf,sizeof(unsigned char),nr*nc,fp);
    if (i != nr*nc) {
     fprintf(stderr,"fread() returned %d characters instead of %d\n",i,nr*nc);
     perror("fread");
     return NULL;
     };
    for(i=0;i<nr;i++) for(j=0;j<nc;j++) PGM_PIXEL(im,i,j)=buf[i*nc+j];
    free(buf);
    break;

   case '3': /* PPM file, ASCII format */
    PNM_TYPE(im)=PPM;
    PNM_MAXVAL(im)=getint(fp);
    PNM_DATA(im)=(unsigned int *)malloc(3*nr*nc*sizeof(unsigned int));
    for(i=0;i<nr;i++) for(j=0;j<nc;j++) {
     PPM_PIXEL_R(im,i,j)=getint(fp);
     PPM_PIXEL_G(im,i,j)=getint(fp);
     PPM_PIXEL_B(im,i,j)=getint(fp);
     };
    break;

   case '6': /* PPM file, raw format */
    PNM_TYPE(im)=PPM;
    PNM_MAXVAL(im)=getint(fp);
    buf=(unsigned char *)malloc(3*nr*nc*sizeof(unsigned char));
    PNM_DATA(im)=(unsigned int *)malloc(3*nr*nc*sizeof(unsigned int));
    if (fread(buf,sizeof(unsigned char),3*nr*nc,fp) != 3*nr*nc) {
     perror("fread");
     return NULL;
     };
    for(i=0;i<nr;i++) for(j=0;j<nc;j++) {
      PPM_PIXEL_R(im,i,j)=buf[(i*nc+j)*3];
      PPM_PIXEL_G(im,i,j)=buf[(i*nc+j)*3+1];
      PPM_PIXEL_B(im,i,j)=buf[(i*nc+j)*3+2];
    }
    free(buf);
    break;

  default: 
   fprintf(stderr,"Bogus PNM file header 'P%c'\n",pnmfiletype);
   return NULL;
  };
 
  return im;
}

void write_pnm_image(pnm_image *im,FILE *fp)
{
  int nr=PNM_NR(im),nc=PNM_NC(im),i,j;


  switch(PNM_TYPE(im)) {
   case PBM_RAW: 
    fprintf(fp,"P4\n# moo\n%d %d\n",nc,nr);
    for(i=0;i<nr;i++) {
     unsigned char c=0,mask=0x80;
     for(j=0;j<nc;j++) {
      if (PBM_PIXEL(im,i,j)) c |= mask;
      mask >>= 1;
      if (!mask) { putc(c,fp); mask=0x80; c=0; };
      };
     if (mask!=0x80) {putc(c,fp);};
     };
    break;

   case PBM_ASCII:
    fprintf(fp,"P1\n%d %d\n",nc,nr);
    for(i=0;i<nr;i++) {
     for(j=0;j<nc;j++) {
      putc(PBM_PIXEL(im,i,j)+'0',fp);
      if (!(j%65)) putc('\n',fp);
      };
     putc('\n',fp);
     };
    break;

   case PGM_RAW: 
    fprintf(fp,"P5\n%d %d\n%d\n",nc,nr,PNM_MAXVAL(im));
    for(i=0;i<nr;i++) for(j=0;j<nc;j++) {
      unsigned char pix = PGM_PIXEL(im,i,j);
      fwrite(&pix,1,1,fp);
    }
    break;

   case PGM_ASCII:
    fprintf(fp,"P2\n%d %d\n%d\n",nc,nr,PNM_MAXVAL(im));
    for(i=0;i<nr;i++) {
     for(j=0;j<nc;j++) {
      fprintf(fp,"%d ",PGM_PIXEL(im,i,j));
      if (!(j%16)) fprintf(fp,"\n");
      };
     fprintf(fp,"\n");
     };
    break;

   case PPM_RAW:
    fprintf(fp,"P6\n%d %d\n%d\n",nc,nr,PNM_MAXVAL(im));
    for(i=0;i<nr;i++) for(j=0;j<nc;j++) {
      unsigned char pix = PPM_PIXEL_R(im,i,j);
      fwrite(&pix,1,1,fp);
      pix = PPM_PIXEL_G(im,i,j);
      fwrite(&pix,1,1,fp);
      pix = PPM_PIXEL_B(im,i,j);
      fwrite(&pix,1,1,fp);
    }
    break;

   case PPM_ASCII:
    fprintf(fp,"P3\n%d %d\n%d\n",nc,nr,PNM_MAXVAL(im));
    for(i=0;i<nr;i++) {
     for(j=0;j<nc;j++) {
      fprintf(fp,"%d %d %d ",PPM_PIXEL_R(im,i,j),
                             PPM_PIXEL_G(im,i,j),PPM_PIXEL_B(im,i,j));
      if (!(j%5)) fprintf(fp,"\n");
      };
     fprintf(fp,"\n");
     };
    break;
   default:
    fprintf(stderr,"bogus value (%d) for pnm_type\n",PNM_TYPE(im));
   };
}

void write_pnm_image_to_file(pnm_image *im,char *filename)
{
 int l=strlen(filename);
 FILE *fp;
 char cmd[256];

 /* .Z ? */
 if ((filename[l-1]=='Z') && (filename[l-2] == '.')) 
  sprintf(cmd,"compress -f > %s",filename);
 else if ((filename[l-1]=='z') && (filename[l-2]=='g') && (filename[l-3]=='.'))
  sprintf(cmd,"gzip -9 > %s",filename);
 else
  sprintf(cmd,"cat > %s",filename);

 fp=popen(cmd,"w");
 if (!fp) {
  perror("popen");
  return;
  };

 write_pnm_image(im,fp);

 pclose(fp);
}

void cleanup_pnm_image(pnm_image *im)
{
  free(PNM_DATA(im));
  free(im);
}

void clear_pnm_image(pnm_image *im)
{
  int nr=PNM_NR(im), nc=PNM_NC(im);
  if (IS_PBM(im) || IS_PGM(im)) bzero(PNM_DATA(im),nr*nc*sizeof(unsigned int));
  else if (IS_PPM(im)) bzero(PNM_DATA(im),nr*nc*3*sizeof(unsigned int));
  else {
   fprintf(stderr,"Error: bogus PNM image type %d\n",PNM_TYPE(im));
   exit(-1);
   };
}

pnm_image *dup_pnm_image(pnm_image *im)
{
  pnm_image *im2;
  int nr=PNM_NR(im), nc=PNM_NC(im);

  if (IS_PBM(im)) {
   im2=allocate_pbm_image(nr,nc);
   bcopy(PNM_DATA(im),PNM_DATA(im2),nr*nc*sizeof(unsigned int));
   }
  else if (IS_PGM(im)) {
   im2=allocate_pgm_image(nr,nc);
   bcopy(PNM_DATA(im),PNM_DATA(im2),nr*nc*sizeof(unsigned int));
   }
  else if (IS_PPM(im)) {
   im2=allocate_ppm_image(nr,nc);
   bcopy(PNM_DATA(im),PNM_DATA(im2),3*nr*nc*sizeof(unsigned int));
   }
  else {
   fprintf(stderr,"Error: bogus PNM image type %d\n",PNM_TYPE(im));
   exit(-1);
   };

  /* preserve exact type */
  PNM_TYPE(im2)=PNM_TYPE(im);
  return im2;
}
