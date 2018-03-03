APPNAME = calculatePI.out
TARGETS = main.o zeta1.o

CC = mpiCC
MPI = mpirun
RM=rm -f
CFLAGS = -Wall -g
LDFLAGS = -O1

all: $(TARGETS)
	$(CC) $(CFLAGS) -o $(APPNAME) $(TARGETS) $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $<

run: all
	$(MPI) -np 2 $(APPNAME)

utest: all 
	$(MPI) -np 2 $(APPNAME) -u

.PHONY: clean
clean:
	$(RM) $(APPNAME) $(TARGETS)