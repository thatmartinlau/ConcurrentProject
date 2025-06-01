CXX = g++
CXXFLAGS = -O3 -std=c++17 -fopenmp $(shell pkg-config --cflags Magick++)
LDFLAGS = -fopenmp $(shell pkg-config --libs Magick++) -lfftw3

NTHREADS ?= 5  # Default value if not specified
VISUALIZE ?= false  # Default value
PRINT_TELEMETRY ?= false  # Default value
EXPORT_CSV ?= false  # Default value

nbody: 
	rm -f rockyplanets.mp4 \
	&& $(CXX) $(CXXFLAGS) \
		-DN_THREADS=$(NTHREADS) \
		-DVISUALIZE=$(VISUALIZE) \
		-DPRINT_TELEMETRY=$(PRINT_TELEMETRY) \
		-DEXPORT_CSV=$(EXPORT_CSV) \
		-o nbody_sim \
		main.cpp \
		src/core.cpp \
		src/simplesimulation.cpp \
		src/barneshutt.cpp \
		$(LDFLAGS) \
		&& ./nbody_sim \
		&& rm -f nbody_sim \
		&& echo "Done!"

nbody_particle: 
	$(CXX) $(CXXFLAGS) \ 
		-DN_THREADS=$(NTHREADS) \
		-DVISUALIZE=$(VISUALIZE) \
		-DPRINT_TELEMETRY=$(PRINT_TELEMETRY) \
		-DEXPORT_CSV=$(EXPORT_CSV) \
		-o nbody_sim \
		main.cpp \
		mainparticlemesh.cpp \
		src/core.cpp \
		src/simplesimulation.cpp \
		src/barneshutt.cpp \
		src/particlemesh.cpp \
		$(LDFLAGS) \
		&& ./nbody_sim \
		&& rm nbody_sim

bs:
	$(CXX) $(CXXFLAGS) \
		-DN_THREADS=$(NTHREADS) \
		-DVISUALIZE=$(VISUALIZE) \
		-DPRINT_TELEMETRY=$(PRINT_TELEMETRY) \
		-DEXPORT_CSV=$(EXPORT_CSV) \ 
		-o nbody_bs \
		main.cpp \
		src/core.cpp \
		src/simplesimulation.cpp \
		src/barneshutt.cpp \
		$(LDFLAGS)
	@echo "==== Running Barnes-Hut (serial) ===="
	OMP_NUM_THREADS=13 ./nbody_bs -method=barneshut
	rm -f nbody_bs

bsp:
	$(CXX) $(CXXFLAGS)\ 
		-DN_THREADS=$(NTHREADS) \
		-DVISUALIZE=$(VISUALIZE) \
		-DPRINT_TELEMETRY=$(PRINT_TELEMETRY) \
		-DEXPORT_CSV=$(EXPORT_CSV) \
		-o nbody_bsp \
		main.cpp \
		src/core.cpp \
		src/simplesimulation.cpp \
		src/barneshutt.cpp \
		$(LDFLAGS)
	@echo "==== Running Barnes-Hut (parallel, 13 threads) ===="
	OMP_NUM_THREADS=13 ./nbody_bsp -method=barneshut -parallel
	rm -f nbody_bsp


.PHONY: clean

clean:
	rm -f nbody_sim *.gif telemetry nbody_sim_particle_mesh *.csv \
	&& rm -rf *frames*/

