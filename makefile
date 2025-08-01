BIN=./bin
BUILD=./build
SRC=./src

CC=gcc
FC=gfortran
FFLAGS=-fopenmp -O3 -ffast-math -funroll-loops

export OMP_NUM_THREADS=4

obj = $(BUILD)/readparam.o $(BUILD)/readgridparam.o $(BUILD)/readdensitygrid.o $(BUILD)/numberofmuogdata.o $(BUILD)/readmuogdata.o $(BUILD)/readtopoparam.o $(BUILD)/readtopo.o $(BUILD)/bilinint.o $(BUILD)/trilinintcoeffs.o $(BUILD)/panda.o
obj2 = $(BUILD)/readparam.o $(BUILD)/readgridparam.o $(BUILD)/readdensitygrid.o $(BUILD)/numberofmuogdata.o $(BUILD)/readmuogdata2.o $(BUILD)/readtopoparam.o $(BUILD)/readtopo.o $(BUILD)/intersfacex.o $(BUILD)/intersfacey.o $(BUILD)/intersfacez.o $(BUILD)/interscube.o $(BUILD)/bilinint.o $(BUILD)/trilinintcoeffs2.o $(BUILD)/panda2.o

.SUFFIXES: .o .f90 .mod

$(BUILD)/%.o: $(SRC)/%.f90
	$(FC) $(FLAGS) -c -o $@ $^

$(BUILD)/%.o: $(SRC)/%.c
	$(CC) $(FLAGS) -c -o $@ $^

panda: $(obj)
	$(FC) $(FFLAGS) -o $(BIN)/panda $(obj)

panda2: $(obj2)
	$(FC) $(FFLAGS) -o $(BIN)/panda2 $(obj2)

clean:
	/bin/rm -f $(BUILD)/*.o $(BIN)/*

