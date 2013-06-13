FLAGS = -fopenmp -O3

all: main/test main/graphgen main/parametertune
main/parametertune: main/parametertune.o growingnetwork3d/growingnetwork3d.o growingnetwork3d/spatialvertex.o graph/growingnetwork.o graph/graph.o graph/vertex.o
<<<<<<< HEAD
	g++ $(FLAGS)  main/parametertune.o growingnetwork3d/growingnetwork3d.o growingnetwork3d/spatialvertex.o graph/growingnetwork.o graph/graph.o graph/vertex.o -o main/parametertune 
=======
	g++ main/parametertune.o growingnetwork3d/growingnetwork3d.o growingnetwork3d/spatialvertex.o graph/growingnetwork.o graph/graph.o graph/vertex.o -o main/parametertune $(FOPENMP)
>>>>>>> 71f8691a972e8247abeb6f6b31b24d8d992dd458
main/test: main/test.o growingnetwork3d/growingnetwork3d.o growingnetwork3d/spatialvertex.o graph/growingnetwork.o graph/graph.o graph/vertex.o
	g++ $(FLAGS)  main/test.o growingnetwork3d/growingnetwork3d.o growingnetwork3d/spatialvertex.o graph/growingnetwork.o graph/graph.o graph/vertex.o -o main/test 
main/graphgen: main/graphgen.o growingnetwork3d/growingnetwork3d.o growingnetwork3d/spatialvertex.o growingnetwork2d/growingnetwork2d.o graph/growingnetwork.o graph/graph.o graph/vertex.o
	g++ $(FLAGS)  main/graphgen.o growingnetwork3d/growingnetwork3d.o growingnetwork3d/spatialvertex.o growingnetwork2d/growingnetwork2d.o graph/growingnetwork.o graph/graph.o graph/vertex.o -o main/graphgen 
main/pgrownet: main/pgrownet.o pgrownet2d/pgrownet2d.o graph/growingnetwork.o graph/graph.o graph/vertex.o growingnetwork3d/spatialvertex.o
<<<<<<< HEAD
	g++ $(FLAGS)  main/pgrownet.o pgrownet2d/pgrownet2d.o growingnetwork3d/spatialvertex.o graph/growingnetwork.o graph/graph.o graph/vertex.o -o main/pgrownet 
main/parametertune.o: main/parametertune.cpp growingnetwork3d/growingnetwork3d.cpp growingnetwork3d/growingnetwork3d.h
	g++ $(FLAGS)  -c main/parametertune.cpp -o main/parametertune.o 
main/test.o: main/test.cpp growingnetwork3d/growingnetwork3d.cpp growingnetwork3d/growingnetwork3d.h
	g++ $(FLAGS)  -c main/test.cpp -o main/test.o 
=======
	g++ main/pgrownet.o pgrownet2d/pgrownet2d.o growingnetwork3d/spatialvertex.o graph/growingnetwork.o graph/graph.o graph/vertex.o -o main/pgrownet $(FOPENMP)
main/parametertune.o: main/parametertune.cpp growingnetwork3d/growingnetwork3d.cpp growingnetwork3d/growingnetwork3d.h
	g++ -c main/parametertune.cpp -o main/parametertune.o $(FOPENMP)
main/test.o: main/test.cpp growingnetwork3d/growingnetwork3d.cpp growingnetwork3d/growingnetwork3d.h
	g++ -c main/test.cpp -o main/test.o $(FOPENMP)
>>>>>>> 71f8691a972e8247abeb6f6b31b24d8d992dd458
main/graphgen.o: main/graphgen.cpp growingnetwork2d/growingnetwork2d.cpp growingnetwork2d/growingnetwork2d.h
	g++ $(FLAGS)  -c main/graphgen.cpp -o main/graphgen.o 
main/pgrownet.o: main/pgrownet.cpp pgrownet2d/pgrownet2d.cpp pgrownet2d/pgrownet2d.h
	g++ $(FLAGS)  -c main/pgrownet.cpp -o main/pgrownet.o 
growingnetwork3d/growingnetwork3d.o: growingnetwork3d/growingnetwork3d.cpp growingnetwork3d/growingnetwork3d.h graph/*.cpp graph/*.h
	g++ $(FLAGS)  -c growingnetwork3d/growingnetwork3d.cpp -o growingnetwork3d/growingnetwork3d.o 
growingnetwork3d/spatialvertex.o: growingnetwork3d/spatialvertex.cpp growingnetwork3d/spatialvertex.h graph/vertex.cpp graph/vertex.h
	g++ $(FLAGS)  -c growingnetwork3d/spatialvertex.cpp -o growingnetwork3d/spatialvertex.o 
growingnetwork2d/growingnetwork2d.o: growingnetwork2d/growingnetwork2d.cpp growingnetwork2d/growingnetwork2d.h graph/*.cpp graph/*.h
	g++ $(FLAGS)  -c growingnetwork2d/growingnetwork2d.cpp -o growingnetwork2d/growingnetwork2d.o 
graph/growingnetwork.o: graph/growingnetwork.cpp graph/growingnetwork.h graph/graph.cpp graph/graph.h graph/vertex.cpp graph/vertex.h
	g++ $(FLAGS)  -c graph/growingnetwork.cpp -o graph/growingnetwork.o 
graph/graph.o: graph/graph.cpp graph/graph.h graph/vertex.cpp graph/vertex.h
	g++ $(FLAGS)  -c graph/graph.cpp -o graph/graph.o 
graph/vertex.o: graph/vertex.cpp graph/vertex.h
	g++ $(FLAGS)  -c graph/vertex.cpp -o graph/vertex.o 
pgrownet2d/pgrownet2d.o: pgrownet2d/pgrownet2d.cpp pgrownet2d/pgrownet2d.h growingnetwork3d/spatialvertex.cpp growingnetwork3d/spatialvertex.h graph/*.cpp graph/*.h
	g++ $(FLAGS)  -c pgrownet2d/pgrownet2d.cpp -o pgrownet2d/pgrownet2d.o
backup: graph/*.cpp graph/*.h growingnetwork2d/*.cpp growingnetwork2d/*.h growingnetwork3d/*.cpp growingnetwork3d/*.h main/*.cpp main/*.m makefile README doc.odt
	zip ../backup.zip graph/*.cpp graph/*.h growingnetwork2d/*.cpp growingnetwork2d/*.h growingnetwork3d/*.cpp growingnetwork3d/*.h main/*.cpp main/*.m makefile README doc.odt
clean:
	rm -vf main/test
	rm -vf main/graphgen
	rm -vf main/points.txt
	rm -vf */*.o
