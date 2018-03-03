APPNAME = calculatePI.out
TARGETS = main.o zeta2.o

CC = g++
RM=rm -f
CFLAGS = -fopenmp
LDFLAGS = -O1

all: $(TARGETS)
	$(CC) $(CFLAGS) -o $(APPNAME) $(TARGETS) $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $<

run: all
	OMP_NUM_THREADS=2 ./$(APPNAME)

utest: all 
	./$(APPNAME) -u

vtest: all
	./$(APPNAME) -v

.PHONY: clean
clean:
	$(RM) $(APPNAME) $(TARGETS)