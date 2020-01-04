F90     = ifort
MPIF90  = mpif90
MKLFLAG = -mkl

Similarity.o : Similarity.f90
	$(F90) Similarity.f90 $(MKLFLAG) -c
test.o : test.f90
	$(F90) test.f90 $(MKLFLAG) -c

OBJS_TEST = Similarity.o test.o

test : $(OBJS_TEST)
	$(F90) $(OBJS_TEST) $(MKLFLAG) -o test && ./test

clean :
	rm -rf *.o *.mod