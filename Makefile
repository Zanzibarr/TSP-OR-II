SRCDIR = src
INCDIR = include
CPLEX_DIR = /opt/ibm/ILOG/CPLEX_Studio2211/cplex/

SRCS = $(wildcard $(SRCDIR)/*.c)
OBJS = $(SRCS:$(SRCDIR)/%.c=%.o)
DEPS = $(wildcard $(INCDIR)/*.h)
TARGET = tsp

CC = gcc
CFLAGS = -I./${INCDIR} -I${CPLEX_DIR}/include/ilcplex

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(TARGET) -L${CPLEX_DIR}/lib/x86-64_linux/static_pic -lcplex -lm -g && make clean

%.o: $(SRCDIR)/%.c $(DEPS)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) callgrind.out.* *.lp

clear:
	rm -f $(OBJS) $(TARGET) callgrind.out.* *.lp
