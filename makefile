CXX = g++
CXXFLAGS = -O3 -std=c++17 -fopenmp $(shell pkg-config --cflags Magick++)
LDFLAGS = -fopenmp $(shell pkg-config --libs Magick++) -lfftw3

# Main target
nbody: 
	rm -f rockyplanets.mp4 \
	&& $(CXX) $(CXXFLAGS) -o nbody_sim \
		main.cpp \
		src/core.cpp \
		src/simplesimulation.cpp \
		src/barneshutt.cpp \
		$(LDFLAGS) \
		&& ./nbody_sim \
		&& rm -f nbody_sim \
		&& echo "Done!" 
		
nbody_particle: 
	$(CXX) $(CXXFLAGS) -o nbody_sim \
		main.cpp \
		mainparticlemesh.cpp \
		src/core.cpp \
		src/simplesimulation.cpp \
		src/barneshutt.cpp \
		src/particlemesh.cpp \
		$(LDFLAGS) \
		&& ./nbody_sim \
		&& rm nbody_sim

telemetry: 
	$(CXX) $(CXXFLAGS) -o telemetry \
		test_telemetry.cpp \
		src/core.cpp \
		src/simplesimulation.cpp \
		src/barneshutt.cpp \
		$(LDFLAGS) \
		&& ./telemetry \
		&& rm telemetry

.PHONY: clean

clean:
	rm -f nbody_sim *.gif telemetry nbody_sim_particle_mesh *.csv \
	&& rm -rf *frames*/

