APPNAME = calculatePI.out
TARGETS = main.o zeta0.o

CC = g++
RM=rm -f
CFLAGS = -Wall -g
LDFLAGS = -O1

all: $(TARGETS)
	$(CC) $(CFLAGS) -o $(APPNAME) $(TARGETS) $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $<

run: all
	./$(APPNAME)

utest: all 
	./$(APPNAME) -u

vtest: all
	./$(APPNAME) -v

.PHONY: clean
clean:
	$(RM) $(APPNAME) $(TARGETS)