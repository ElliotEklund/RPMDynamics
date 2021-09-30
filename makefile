#build an executable named RPMDynamics from main.cpp
all: main.cpp
	mpic++ -O3 -o RPMDynamics  main.cpp functions.cpp
clean:
	$(RM) RPMDynamics
