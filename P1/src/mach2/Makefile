APPNAME = calculatePI.out
TARGETS = main.o mach2.o

CC = g++
MPI = mpirun
RM=rm -f
CFLAGS = -fopenmp
LDFLAGS = -O1

build: $(TARGETS)
	$(CC) $(CFLAGS) -o $(APPNAME) $(TARGETS) $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $<

test: build
	./$(APPNAME)

utest: build 
	./$(APPNAME) -u

vtest: build
	./$(APPNAME) -v

.PHONY: clean
clean:
	$(RM) $(APPNAME) $(TARGETS)