SRCS = $(shell ls *.F90)
OBJS = $(addsuffix .o, $(basename $(SRCS)))

all: $(OBJS)

.SUFFIXES:
.SUFFIXES: .o .F90

.F90.o:
	$(FC) -c $<
	ar rc libgptl.a *.o
clean:
	$(RM) *.o *.a
