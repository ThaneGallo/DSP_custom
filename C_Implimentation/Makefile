# Compiler
CC = gcc

# Compiler flags
CFLAGS = -Wall -Wextra -g

# Source files
SRCS = main.c dsp.c

# Object files
OBJS = $(SRCS:.c=.o)

# Header files (optional, but can be used if you have multiple header files)
DEPS = dsp.h

# Output executable
TARGET = main

# Default target
all: $(TARGET)

# Rule to link object files into the final executable
$(TARGET): $(OBJS)
	$(CC) $(OBJS) -o $(TARGET) -lm  # Link math library here

# Rule to compile .c files to .o files
%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c $< -o $@

# Clean up object files and executable
clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean  # Marking targets as phony
