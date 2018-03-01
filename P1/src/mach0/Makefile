APPNAME = calculatePI

CC = g++
CFLAGS = -Wall -g
LDFLAGS = -O1

TARGET = main.o zeta0.o

all: $(TARGET)
	$(CC) $(CFLAGS) -o $(APPNAME) $(TARGET) $(LDFLAGS)

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
	rm -rf $(APPNAME) $(TARGET)