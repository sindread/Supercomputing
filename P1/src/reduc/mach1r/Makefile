APPNAME = calculatePI.out
TARGETS = main.o mach1r.o

CC = mpiCC
MPI = mpirun
RM=rm -f
CFLAGS = -Wall -g
LDFLAGS = -O1

build: $(TARGETS)
	$(CC) $(CFLAGS) -o $(APPNAME) $(TARGETS) $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $<

test: build
	$(MPI) -np 2 $(APPNAME)

utest: build 
	$(MPI) -np 2 $(APPNAME) -u

.PHONY: clean
clean:
	$(RM) $(APPNAME) $(TARGETS)