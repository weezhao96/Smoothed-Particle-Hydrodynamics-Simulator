CXX=mpicxx
CXXFLAGS=-Wall
LIBARGS=-l boost_filesystem -l boost_system -l blas

default: SPH_sim

main.o: main.cpp SPH.h
	$(CXX) $(CXXFLAGS) -o main.o -c main.cpp

SPH.o: SPH.cpp SPH.h
	$(CXX) $(CXXFLAGS) -o SPH.o -c SPH.cpp

SPH_init.o: SPH_init.cpp SPH.h
	$(CXX) $(CXXFLAGS) -o SPH_init.o -c SPH_init.cpp

MPI_Manager.o: MPI_Manager.cpp MPI_Manager.h
	$(CXX) $(CXXFLAGS) -o MPI_Manager.o -c MPI_Manager.cpp

IO_Manager.o: IO_Manager.cpp IO_Manager.h
	$(CXX) $(CXXFLAGS) -o IO_Manager.o -c IO_Manager.cpp

SPH_sim: main.o SPH.o SPH_init.o MPI_Manager.o IO_Manager.o
	$(CXX) $(CXXFLAGS) -o SPH_main.out main.o SPH.o SPH_init.o MPI_Manager.o IO_Manager.o $(LIBARGS)

.PHONY: clean

clean:
	rm -f *.o SPH_sim