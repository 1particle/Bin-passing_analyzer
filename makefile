 CC = g++
 CFLAGS  = -g -Wall
 main = src/Bin-passing_analyzer
 util = src/PMF_calculation
 bin = bin_Linux_x86_64/Bin-passing_analyzer__Linux_x86_64


$(bin):	Bin-passing_analyzer.o PMF_calculation.o
	$(CC) $(CFLAGS)	-o	$(bin)	Bin-passing_analyzer.o	PMF_calculation.o

Bin-passing_analyzer.o:  $(main).cpp
	$(CC) $(CFLAGS) -c $(main).cpp

PMF_calculation.o: $(util).cpp  $(util).h 
	$(CC) $(CFLAGS) -c $(util).cpp

clean: 
	$(RM) count *.o *~
