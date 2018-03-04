APPNAME = calculatePI.out
TARGETS = main.o mach0.o

CC = g++
RM=rm -f
CFLAGS = -Wall -g
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