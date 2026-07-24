# Makefile for 1-D SSW Model
# Requires: gfortran (or ifort) + NetCDF-Fortran

FC       ?= gfortran
FFLAGS   ?= -O3 -ffree-line-length-none
SRC       = src/fortran/1d_model.f90
EXE       = d.out

# NetCDF flags (auto-detect via nf-config or nc-config)
NFCONFIG := $(shell command -v nf-config 2>/dev/null)
NCCONFIG := $(shell command -v nc-config 2>/dev/null)

ifdef NFCONFIG
  NC_INC  := $(shell nf-config --fflags)
  NC_LIBS := $(shell nf-config --flibs)
else ifdef NCCONFIG
  NC_INC  := -I$(shell nc-config --includedir)
  NC_LIBS := -lnetcdff $(shell nc-config --libs)
else ifdef NETCDF_ROOT
  NC_INC  := -I$(NETCDF_ROOT)/include
  NC_LIBS := -L$(NETCDF_ROOT)/lib -lnetcdff -lnetcdf
else
  $(error Cannot find NetCDF. Set NETCDF_ROOT or ensure nf-config/nc-config is in PATH)
endif

.PHONY: all clean

all: $(EXE)

$(EXE): $(SRC)
	$(FC) $(FFLAGS) $(NC_INC) $< $(NC_LIBS) -o $@

clean:
	rm -f $(EXE) *.o *.mod
