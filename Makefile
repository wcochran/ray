CC=gcc
COPTS=-g -Wall -DVERBOSE
#COPTS=-O3
YACC=yacc
#YACCOPTS=-t -v
YACCOPTS=
LEX=lex

all: ray

CLEAN=*~ *.o core *.dSYM *.yy.c *.tab.[ch] y.output ray.tar.gz

clean:
	-rm -rf $(CLEAN)

clobber:
	-rm -rf $(CLEAN) ray

#LIBS=-lm -lfl
LIBS=-lm -ll

ray: y.tab.o lex.yy.o raytrace.o pnmio.o sphere.o plane.o teapot.o bezier3.o \
            func2.o hermite.o rectgrid.o bbox.o noise.o marble.o checker.o \
            bilinear.o rectfis.o bilinearfis.o hermitefis.o ris.o \
            unimodalroot.o superellipsoid.o pointcloud.o
	$(CC) $(COPTS) $^ -o $@ $(LIBS)

y.tab.o: ray.y raytrace.h pnmio.h sphere.h
	$(YACC) -d $(YACCOPTS) $<
	$(CC) -c $(COPTS) y.tab.c

lex.yy.o: ray.l ray.y
	$(LEX) $<
	$(CC) -c $(COPTS) lex.yy.c

.c.o:
	$(CC) -c $(COPTS) $<

raytrace.o: raytrace.c raytrace.h pnmio.h
pnmio.o: pnmio.c pnmio.h
sphere.o: sphere.c sphere.h raytrace.h trim.h
plane.o: plane.c plane.h raytrace.h pnmio.h
teapot.o: teapot.c teapot.h bezier3.h bbox.h raytrace.h
bezier3.o: bezier3.c bezier3.h trim.h raytrace.h
bbox.o: bbox.c bbox.h raytrace.h
func2.o: func2.c func2.h 
hermite.o: hermite.c hermite.h func2.h 
rectgrid.o: rectgrid.c rectgrid.h raytrace.h bbox.h pnmio.h
bilinear.o: bilinear.c bilinear.h func2.h
rectfis.o: rectfis.c rectfis.h func2.h
bilinearfis.o: bilinearfis.c bilinearfis.h bilinear.h rectfis.h func2.h
hermitefis.o: hermitefis.c hermitefis.h hermite.h bilinear.h rectfis.h func2.h
ris.o: ris.c ris.h func2.h hermite.h bilinear.h
marble.o: marble.c marble.h noise.h raytrace.h
checker.o: checker.c checker.h noise.h raytrace.h
noise.o: noise.c noise.h
unimodalroot.o: unimodalroot.c unimodalroot.h
superellipsoid.o: superellipsoid.c superellipsoid.h unimodalroot.h raytrace.h bbox.h
pointcloud.o: pointcloud.c pointcloud.h bbox.h

testunimodalroot: testunimodalroot.c unimodalroot.c unimodalroot.h
	$(CC) $(COPTS) testunimodalroot.c unimodalroot.c -o $@

archive: ray.tar.gz

ray.tar.gz:
	tar czvf $@ README Makefile *.[chyl] input-samples
