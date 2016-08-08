Alignment: main.o NW.o 	
	g++ main.o NW.o  -o NeedlemanWunch

main.o: main.cpp
	g++ -c main.cpp

NW.o: NW.cpp
	g++ -c NW.cpp


clean: 
	rm *.o  NeedlemanWunch
