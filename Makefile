CC = gcc
CFLAGS = -I./include

SRCDIR = src
INCDIR = include

SRCS = $(wildcard $(SRCDIR)/*.c)
OBJS = $(SRCS:$(SRCDIR)/%.c=%.o)
DEPS = $(wildcard $(INCDIR)/*.h)
TARGET = main

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(TARGET) -lm && make clean

%.o: $(SRCDIR)/%.c $(DEPS)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS)

clear:
	rm -f $(OBJS) $(TARGET)