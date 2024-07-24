SRCDIR = src
INCDIR = include

CPLEX_DIR = /Applications/CPLEX_Studio2211/cplex
CPLEX_PIC_DIR = ${CPLEX_DIR}/lib/arm64_osx

SRCS = $(wildcard $(SRCDIR)/*.c)
OBJS = $(SRCS:$(SRCDIR)/%.c=%.o)
DEPS = $(wildcard $(INCDIR)/*.h)
TARGET = tsp

CC = gcc
CFLAGS = -I./${INCDIR} -I${CPLEX_DIR}/include/ilcplex

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(TARGET) -L${CPLEX_PIC_DIR}/static_pic -lcplex -lm -O3 && make clean

%.o: $(SRCDIR)/%.c $(DEPS)
	$(CC) $(CFLAGS) -c $< -o $@

debug: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(TARGET) -L${CPLEX_PIC_DIR}/static_pic -lcplex -lm -g && make clean

clean:
	rm -f $(OBJS) callgrind.out.*

clear:
	rm -f $(OBJS) $(TARGET) callgrind.out.*
