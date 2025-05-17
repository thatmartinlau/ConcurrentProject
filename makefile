CXX = g++
CXXFLAGS = -O3 -std=c++17 -fopenmp $(shell pkg-config --cflags Magick++)
LDFLAGS = $(shell pkg-config --libs Magick++)

# Main target
nbody: 
	$(CXX) $(CXXFLAGS) -o nbody_sim \
		main.cpp \
		src/core.cpp \
		src/simplesimulation.cpp \
		src/barneshutt.cpp \
		src/particlemesh.cpp \
		$(LDFLAGS) \
		&& ./nbody_sim \
		&& rm nbody_sim

.PHONY: clean

clean:
	rm -f nbody_sim *.gif
