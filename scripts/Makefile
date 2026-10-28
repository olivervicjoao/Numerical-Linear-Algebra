
CXX = g++

INC = -I./include

# ---------- FLAGS ----------
CXXFLAGS = -O3 -std=c++17 $(INC)

# ---------- LIBRARIES ----------
LDLIBS = -llapacke -lopenblas -lm

# ---------- SOURCES ----------
SRC_CPP := $(wildcard ./src/*.cc)

# ---------- TARGETS ----------
all: nla

nla: main.cc $(SRC_CPP)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

# ---------- CLEAN ----------
clean:
	rm -f nla

.PHONY: all clean
