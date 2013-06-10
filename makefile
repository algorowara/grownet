FOPENMP = -fopenmp

main/graphgen: main/graphgen.o growingnetwork2d/growingnetwork2d.o graph/growingnetwork.o graph/graph.o graph/vertex.o
	g++ main/graphgen.o growingnetwork2d/growingnetwork2d.o graph/growingnetwork.o graph/graph.o graph/vertex.o -o main/graphgen $(FOPENMP)
main/graphgen.o: main/graphgen.cpp growingnetwork2d/growingnetwork2d.cpp growingnetwork2d/growingnetwork2d.h
	g++ -c main/graphgen.cpp -o main/graphgen.o $(FOPENMP)
growingnetwork3d/growingnetwork3d.o: growingnetwork3d/growingnetwork3d.cpp growingnetwork3d/growingnetwork3d.h graph/*.cpp graph/*.h
	g++ -c growingnetwork3d/growingnetwork3d.cpp -o growingnetwork3d/growingnetwork3d.o
growingnetwork3d/spatialvertex.o: growingnetwork3d/spatialvertex.cpp growingnetwork3d/spatialvertex.h graph/vertex.cpp graph/vertex.h
	g++ -c growingnetwork3d/spatialvertex.cpp -o growingnetwork3d/spatialvertex.o
growingnetwork2d/growingnetwork2d.o: growingnetwork2d/growingnetwork2d.cpp growingnetwork2d/growingnetwork2d.h graph/*.cpp graph/*.h
	g++ -c growingnetwork2d/growingnetwork2d.cpp -o growingnetwork2d/growingnetwork2d.o
graph/growingnetwork.o: graph/growingnetwork.cpp graph/growingnetwork.h graph/graph.cpp graph/graph.h graph/vertex.cpp graph/vertex.h
	g++ -c graph/growingnetwork.cpp -o graph/growingnetwork.o
graph/graph.o: graph/graph.cpp graph/graph.h graph/vertex.cpp graph/vertex.h
	g++ -c graph/graph.cpp -o graph/graph.o
graph/vertex.o: graph/vertex.cpp graph/vertex.h
	g++ -c graph/vertex.cpp -o graph/vertex.o
backup: graph/*.cpp graph/*.h growingnetwork2d/*.cpp growingnetwork2d/*.h growingnetwork3d/*.cpp growingnetwork3d/*.h makefile README doc.odt
	zip ../backup.zip graph/*.cpp graph/*.h growingnetwork2d/*.cpp growingnetwork2d/*.h growingnetwork3d/*.cpp growingnetwork3d/*.h makefile README doc.odt
clean:
	rm -vf main/graphgen
	rm -vf */*.o
