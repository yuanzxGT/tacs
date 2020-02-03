include ../../Makefile.in
include ../../TACS_Common.mk

OBJS = 	plate_tri.o \
	create_plate_tri.o 

default: ${OBJS}
	${CXX} -o plate_tri plate_tri.o ${TACS_LD_FLAGS}
	${CXX} -o create_plate_tri create_plate_tri.o ${TACS_LD_FLAGS}
	

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

complex: TACS_DEF="-DTACS_USE_COMPLEX"
complex: default

clean:
	rm -f *.o create_plate_tri plate_tri

test: default
	./plate_tri
	mpirun -np 4 ./plate_tri

test_complex: default
	./plate_tri
	mpirun -np 4 ./plate_tri
