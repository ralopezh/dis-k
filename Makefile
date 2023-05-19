# Start of the makefile

SRC_DIR   = src
OBJ_DIR   = obj

FC         = gfortran
FFLAGS     = -O3 -J$(OBJ_DIR) #-fopenmp

#FC         = ifort
#FFLAGS     = -O3 -module $(OBJ_DIR) -qopenmp

SOURCES1 = params.f90 \
	main.f90 \
	load.f90 \
	plasmadf.f90 \
	disp.f90 \
	muller.f90 \
	solve.f90 \
	solve2d.f90 \
	polarization.f90 \
	extras.f90 \
	cubic_fit.f90 \
	zk12.f90 \
	integrals2.f90 \
	quadpack_double.f90
SOURCES2 = params.f90 \
	nyquist.f90 \
	load.f90 \
	plasmadf.f90 \
	disp.f90 \
	extras.f90 \
	zk12.f90 \
	integrals2.f90 \
	quadpack_double.f90
	
OBJECTS1 = $(SOURCES1:.f90=.o)	
OBJECTS2 = $(SOURCES2:.f90=.o)

SRC_LIST1 = $(addprefix $(SRC_DIR)/,$(SOURCES1))
OBJ_LIST1 = $(addprefix $(OBJ_DIR)/,$(OBJECTS1))
SRC_LIST2 = $(addprefix $(SRC_DIR)/,$(SOURCES2))
OBJ_LIST2 = $(addprefix $(OBJ_DIR)/,$(OBJECTS2))

all: dis-k nyquist
	
dis-k: $(OBJ_LIST1)
	@echo Compiling $@
	@$(FC) $(FFLAGS) $^ -o $@

nyquist: $(OBJ_LIST2)
	@echo Compiling $@
	@$(FC) $(FFLAGS) $^ -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	@mkdir -p $(OBJ_DIR)
	@echo Compiling and generating objects $@ ...
	@$(FC) -c $(FFLAGS) $< -o $@

clean:
	rm $(OBJ_DIR)/*.o
	rm $(OBJ_DIR)/*.mod
	rm dis-k nyquist
# End of the makefile
