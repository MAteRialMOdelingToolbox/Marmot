# -----------------------------------------------------------------------------
# Configuration Variables
# Adjust these paths to match your system installation
# -----------------------------------------------------------------------------

# Compiler to use
CXX := g++

# C++ Standard (Marmot typically requires C++20)
CXX_STD := -std=c++20

# Optimization and Warning flags
CXXFLAGS := $(CXX_STD) -O3 -Wall -Wextra -fPIC

# Installation Directories (Change these if Marmot/libs are not in standard paths)
# We assume that Marmot is installed in the Python PATH
PYTHON := python
PYTHON_PREFIX := $(shell $(PYTHON) -c "import sys; print(sys.prefix)")

MARMOT_ROOT := $(PYTHON_PREFIX)
EIGEN_DIR   := $(PYTHON_PREFIX)/include/eigen3
FASTOR_DIR  := $(PYTHON_PREFIX)/include/Fastor

# Include Paths
INCLUDES := -I$(MARMOT_ROOT)/include \
            -I$(EIGEN_DIR) \
            -I$(FASTOR_DIR)

# Library Paths
LDFLAGS := -L$(MARMOT_ROOT)/lib

# Libraries to link
# -lMarmot: The core library
# -Wl,-rpath: Tells the executable where to look for the .so file at runtime
LDLIBS := -lMarmot -Wl,-rpath,$(MARMOT_ROOT)/lib

# -----------------------------------------------------------------------------
# Build Rules
# -----------------------------------------------------------------------------

# Name of the output executable
TARGET := solver_example

# Source files (assumes your main code is in main.cpp)
SRCS := example.cpp

# Object files (automatically generated names)
OBJS := $(SRCS:.cpp=.o)

# Default target
all: $(TARGET)

# Link the executable
$(TARGET): $(OBJS)
	@echo "Linking $(TARGET)..."
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET) $(LDFLAGS) $(LDLIBS)
	@echo "Build successful! Run with: ./$(TARGET)"

# Compile source files to object files
%.o: %.cpp
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Clean up build artifacts
clean:
	rm -f $(OBJS) $(TARGET)

# Phony targets prevent conflicts with files named 'clean' or 'all'
.PHONY: all clean
