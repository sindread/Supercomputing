APPNAME = calculatePI.out
TARGETS = main.o mach1.o

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
	$(MPI) -np 4 $(APPNAME)

utest: all 
	./$(APPNAME) -u

vtest: all
	./$(APPNAME) -v

.PHONY: clean
clean:
	$(RM) $(APPNAME) $(TARGETS)