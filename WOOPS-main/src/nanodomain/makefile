COMP=gfortran
OPT= -O3 -fPIC -fopenmp

#COMP= ifort
#OPT= -fPIC -shared-intel -mcmodel=large -fast -heap-arrays
##DEBUGGER FOR ifort INTEL COMPILER
##OPT= -fPIC -shared-intel -mcmodel=large -fast -O0 -g -traceback -fpe:0 -check all

OBJ = copyright.o modules.o main.o mainLoop.o aldeArrays.o dcdread.o psfread.o pdbread.o atominfo.o specArrays.o general.o restart.o writeStuff.o userDefinePlugins.o
all: $(OBJ)
	$(COMP) $(OPT) $(LINK) $(OPT1) -o map.e $(OBJ)
copyright.o: copyright.f90
	$(COMP) $(OPT) -c copyright.f90 -o copyright.o
modules.o: modules.f90
	$(COMP) $(OPT) -c modules.f90 -o modules.o
main.o: main.f90
	$(COMP) $(OPT) -c main.f90 -o main.o
mainLoop.o: mainLoop.f90
	$(COMP) $(OPT) -c mainLoop.f90 -o mainLoop.o
aldeArrays.o: aldeArrays.f90
	$(COMP) $(OPT) -c aldeArrays.f90 -o aldeArrays.o
dcdread.o: dcdread.f90
	$(COMP) $(OPT) -c dcdread.f90 -o dcdread.o
psfread.o: psfread.f90
	$(COMP) $(OPT) -c psfread.f90 -o psfread.o
pdbread.o: pdbread.f90
	$(COMP) $(OPT) -c pdbread.f90 -o pdbread.o
atominfo.o: atominfo.f90
	$(COMP) $(OPT) -c atominfo.f90 -o atominfo.o
specArrays.o: specArrays.f90
	$(COMP) $(OPT) -c specArrays.f90 -o specArrays.o
general.o: general.f90
	$(COMP) $(OPT) -c general.f90 -o general.o
restart.o: restart.f90
	$(COMP) $(OPT) -c restart.f90 -o restart.o
writeStuff.o: writeStuff.f90
	$(COMP) $(OPT) -c writeStuff.f90 -o writeStuff.o
userDefinePlugins.o: userDefinePlugins.f90
	$(COMP) $(OPT) $(OPT1) -c userDefinePlugins.f90 -o userDefinePlugins.o
#userDefinePluginsC.o: userDefinePlugins.cu
#	$(CCOMP) $(OPTC) -c userDefinePlugins.cu -o userDefinePluginsC.o
clean:
	rm -rf *.o *.mod
