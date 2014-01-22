PROGRAM = out
C_FILES = $(wildcard src/*.cc)
INCS = -Iinclude -Idlib-18.5

OBJS = $(patsubst %.cc, %.o, $(C_FILES))
CC = g++
CFLAGS = -Wall -Wextra -O2 -std=gnu++0x $(INCS)
LDFLAGS =

all: $(PROGRAM)

$(PROGRAM): depend $(OBJS)
		$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) -o $(PROGRAM)

depend: .depend

.depend: $(C_FILES)
		rm -f ./.depend > /dev/null
		$(CC) $(CFLAGS) -MM $^ > ./.depend;

include .depend

%.o: %.cc
		$(CC) $(CFLAGS) -c $< -o $@

%: %.cc
		$(CC) $(CFLAGS) -o $@ $<

clean:
		@rm -f src/*.o $(PROGRAM)

.PHONY: clean depend
