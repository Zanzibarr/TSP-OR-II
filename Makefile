OS = MAC

SRCDIR = src
INCDIR = include

ifeq (${OS}, MAC)
    CPLEX_DIR = /Applications/CPLEX_Studio2211/cplex
	CPLEX_PIC_DIR = ${CPLEX_DIR}/lib/arm64_osx
else
    CPLEX_DIR = /opt/ibm/ILOG/CPLEX_Studio2211/cplex
	CPLEX_PIC_DIR = ${CPLEX_DIR}/lib/x86-64_linux
endif

SRCS = $(wildcard $(SRCDIR)/*.c)
OBJS = $(SRCS:$(SRCDIR)/%.c=%.o)
DEPS = $(wildcard $(INCDIR)/*.h)
TARGET = tsp

CC = gcc
CFLAGS = -I./${INCDIR} -I${CPLEX_DIR}/include/ilcplex

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(TARGET) -L${CPLEX_PIC_DIR}/static_pic -lcplex -lm -O3 -g && make clean

%.o: $(SRCDIR)/%.c $(DEPS)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) callgrind.out.*

clear:
	rm -f $(OBJS) $(TARGET) callgrind.out.*
