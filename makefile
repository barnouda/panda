BIN=./bin
BUILD=./build
SRC=./src

CC=gcc
FC=gfortran
FFLAGS=-fopenmp -O3 -ffast-math -funroll-loops

export OMP_NUM_THREADS=4

OBJ = $(BUILD)/readparam.o $(BUILD)/readgridparam.o $(BUILD)/readdensitygrid.o $(BUILD)/numberofmuogdata.o $(BUILD)/readmuogdata.o $(BUILD)/readtopoparam.o $(BUILD)/readtopo.o $(BUILD)/intersfacex.o $(BUILD)/intersfacey.o $(BUILD)/intersfacez.o $(BUILD)/interscube.o $(BUILD)/bilinint.o $(BUILD)/trilinintcoeffs.o $(BUILD)/panda.o

.SUFFIXES: .o .f90 .mod

$(BUILD)/%.o: $(SRC)/%.f90
	$(FC) $(FLAGS) -c -o $@ $^

$(BUILD)/%.o: $(SRC)/%.c
	$(CC) $(FLAGS) -c -o $@ $^

panda: $(OBJ)
	$(FC) $(FFLAGS) -o $(BIN)/panda $(OBJ)

clean:
	/bin/rm -f $(BUILD)/*.o $(BIN)/*

