# Makefile – Droplet Heating and Evaporation Simulation
# -------------------------------------------------------

CXX      = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra -Wpedantic
LDFLAGS  = -lm

# Source files shared between targets
SRC_LIB  = droplet_evaporation.cpp

# Targets
SIM_EXE  = droplet_sim
TEST_EXE = droplet_test

.PHONY: all sim test clean

all: sim test

sim: $(SIM_EXE)

test: $(TEST_EXE)

$(SIM_EXE): $(SRC_LIB) main.cpp droplet_evaporation.h
	$(CXX) $(CXXFLAGS) -o $@ $(SRC_LIB) main.cpp $(LDFLAGS)

$(TEST_EXE): $(SRC_LIB) test_droplet.cpp droplet_evaporation.h
	$(CXX) $(CXXFLAGS) -o $@ $(SRC_LIB) test_droplet.cpp $(LDFLAGS)

# Run the simulation
run_sim: $(SIM_EXE)
	./$(SIM_EXE)

# Run unit tests
run_test: $(TEST_EXE)
	./$(TEST_EXE)

# Run everything
run: run_sim run_test

clean:
	rm -f $(SIM_EXE) $(TEST_EXE) results_*.csv
