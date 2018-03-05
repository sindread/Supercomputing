APPNAME = calculatePI.out
TARGETS = main.o mach2.o

CC = mpiCC
MPI = mpirun
RM=rm -f
CFLAGS = -fopenmp
LDFLAGS = -O1

build: $(TARGETS)
	$(CC) $(CFLAGS) -o $(APPNAME) $(TARGETS) $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $<

test: build
	OMP_NUM_THREADS=2 $(MPI) -np 2 $(APPNAME)

utest: build 
	./$(APPNAME) -u

vtest: build
	./$(APPNAME) -v

.PHONY: clean
clean:
	$(RM) $(APPNAME) $(TARGETS)