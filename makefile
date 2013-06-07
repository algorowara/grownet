growingnetwork3d.o: growingnetwork3d/growingnetwork3d.cpp growingnetwork3d/growingnetwork3d.h graph/*.cpp graph/*.h
	g++ -c growingnetwork3d/growingnetwork3d.cpp -o growingnetwork3d/growingnetwork3d.o
spatialvertex.o: growingnetwork3d/spatialvertex.cpp growingnetwork3d/spatialvertex.h graph/vertex.cpp graph/vertex.h
	g++ -c growingnetwork3d/spatialvertex.cpp -o growingnetwork3d/spatialvertex.o
growingnetwork2d.o: growingnetwork2d/growingnetwork2d.cpp growingnetwork2d/growingnetwork2d.h graph/*.cpp graph/*.h
	g++ -c growingnetwork2d/growingnetwork2d.cpp -o growingnetwork2d/growingnetwork2d.o
growingnetwork.o: graph/*.cpp graph/*.h
	g++ -c graph/growingnetwork.cpp -o graph/growingnetwork.o
graph.o: graph/graph.cpp graph/graph.h graph/vertex.cpp graph/vertex.h
	g++ -c graph/graph.cpp -o graph/graph.o
vertex.o: graph/vertex.cpp graph/vertex.h
	g++ -c graph/vertex.cpp -o graph/vertex.o
backup: graph/*.cpp graph/*.h growingnetwork2d/*.cpp growingnetwork2d/*.h growingnetwork3d/*.cpp growingnetwork3d/*.h makefile
	zip ../backup.zip graph/*.cpp graph/*.h growingnetwork2d/*.cpp growingnetwork2d/*.h growingnetwork3d/*.cpp growingnetwork3d/*.h makefile
clean:
	rm -vf */*.o
