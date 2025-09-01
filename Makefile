# Simple Makefile for SVD
CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall
TARGET = svd
SHARED_TARGET = libsvd.so
SOURCE = SVD.cpp

# Try to find Boost automatically
BOOST_INCLUDE = $(shell pkg-config --cflags-only-I libboost-all 2>/dev/null || echo "-I/usr/local/include -I/opt/homebrew/include")
BOOST_LIBS = $(shell pkg-config --libs libboost-all 2>/dev/null || echo "")

# Default target
all: $(TARGET)

$(TARGET): $(SOURCE)
	$(CXX) $(CXXFLAGS) $(BOOST_INCLUDE) $(SOURCE) $(BOOST_LIBS) -o $(TARGET)

shared: $(SOURCE)
	$(CXX) $(CXXFLAGS) -fPIC -shared -DPYTHON_WRAPPER $(BOOST_INCLUDE) $(SOURCE) $(BOOST_LIBS) -o $(SHARED_TARGET)

clean:
	rm -f $(TARGET) $(SHARED_TARGET)

.PHONY: all clean shared
