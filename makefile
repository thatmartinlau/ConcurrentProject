CXX = g++
CXXFLAGS = -O3 -std=c++17 -fopenmp

# Main target
nbody: 
	g++ $(CXXFLAGS) -o nbody_sim \
		main.cpp \
		src/core.cpp \
		src/simplesimulation.cpp \
		src/barneshutt.cpp \
		&& ./nbody_sim \
		&& rm nbody_sim

.PHONY: clean

clean:
	rm -f nbody_sim