FLAGS = -fopenmp -ffast-math -fsingle-precision-constant -O3


all: main/test main/graphgen main/parametertune main/pgrownet main/lineppnet main/figuredata
main/parametertune: main/parametertune.o ngraph/nball.o ngraph/nsphere.o ngraph/ngraph.o growingnetwork3d/growingnetwork3d.o growingnetwork3d/spatialvertex.o graph/growingnetwork.o graph/graph.o graph/vertex.o
	g++ $(FLAGS) main/parametertune.o ngraph/nball.o ngraph/nsphere.o ngraph/ngraph.o growingnetwork3d/growingnetwork3d.o growingnetwork3d/spatialvertex.o graph/growingnetwork.o graph/graph.o graph/vertex.o -o main/parametertune 
main/test: main/test.o ngraph/nball.o ngraph/nsphere.o ngraph/ngraph.o growingnetwork3d/growingnetwork3d.o growingnetwork3d/spatialvertex.o pgrownet2d/pgrownet2d.o growingnetwork2d/growingnetwork2d.o graph/growingnetwork.o graph/graph.o graph/vertex.o
	g++ $(FLAGS) main/test.o ngraph/nball.o ngraph/nsphere.o ngraph/ngraph.o growingnetwork3d/growingnetwork3d.o growingnetwork3d/spatialvertex.o pgrownet2d/pgrownet2d.o growingnetwork2d/growingnetwork2d.o graph/growingnetwork.o graph/graph.o graph/vertex.o -o main/test 
main/graphgen: main/graphgen.o ngraph/nball.o ngraph/nsphere.o ngraph/ngraph.o growingnetwork3d/growingnetwork3d.o growingnetwork3d/spatialvertex.o growingnetwork2d/growingnetwork2d.o graph/growingnetwork.o graph/graph.o graph/vertex.o
	g++ $(FLAGS) main/graphgen.o ngraph/nball.o ngraph/nsphere.o ngraph/ngraph.o growingnetwork3d/growingnetwork3d.o growingnetwork3d/spatialvertex.o growingnetwork2d/growingnetwork2d.o graph/growingnetwork.o graph/graph.o graph/vertex.o -o main/graphgen 
main/pgrownet: main/pgrownet.o ngraph/nball.o ngraph/ngraph.o pgrownet2d/pgrownet2d.o graph/growingnetwork.o graph/graph.o graph/vertex.o growingnetwork3d/spatialvertex.o growingnetwork3d/growingnetwork3d.o
	g++ $(FLAGS) main/pgrownet.o ngraph/nball.o ngraph/ngraph.o pgrownet2d/pgrownet2d.o growingnetwork3d/spatialvertex.o graph/growingnetwork.o graph/graph.o graph/vertex.o growingnetwork3d/growingnetwork3d.o -o main/pgrownet 
main/lineppnet: main/lineppnet.o pp1d/pp1d.o graph/growingnetwork.o graph/graph.o graph/vertex.o growingnetwork3d/spatialvertex.o growingnetwork2d/growingnetwork2d.o
	g++ $(FLAGS) main/lineppnet.o pp1d/pp1d.o growingnetwork3d/spatialvertex.o graph/growingnetwork.o graph/graph.o graph/vertex.o growingnetwork2d/growingnetwork2d.o -o main/lineppnet
main/figuredata: main/figuredata.o pp1d/pp1d.o graph/growingnetwork.o graph/graph.o graph/vertex.o growingnetwork3d/spatialvertex.o growingnetwork2d/growingnetwork2d.o ngraph/nball.o ngraph/ngraph.o ngraph/nsphere.o
	g++ $(FLAGS) main/figuredata.o pp1d/pp1d.o growingnetwork3d/spatialvertex.o graph/growingnetwork.o graph/graph.o graph/vertex.o growingnetwork2d/growingnetwork2d.o ngraph/nball.o ngraph/ngraph.o ngraph/nsphere.o -o main/figuredata
main/figuredata.o: main/figuredata.cpp ngraph/nball.cpp ngraph/nball.h ngraph/nsphere.cpp ngraph/nsphere.h pp1d/pp1d.cpp pp1d/pp1d.h growingnetwork2d/growingnetwork2d.cpp growingnetwork2d/growingnetwork2d.h
	g++ $(FLAGS) -c main/figuredata.cpp -o main/figuredata.o
main/parametertune.o: main/parametertune.cpp growingnetwork3d/growingnetwork3d.cpp growingnetwork3d/growingnetwork3d.h
	g++ $(FLAGS) -c main/parametertune.cpp -o main/parametertune.o 
main/test.o: main/test.cpp ngraph/nball.cpp ngraph/nball.h ngraph/ngraph.cpp ngraph/ngraph.h growingnetwork3d/growingnetwork3d.cpp growingnetwork3d/growingnetwork3d.h growingnetwork2d/growingnetwork2d.cpp growingnetwork2d/growingnetwork2d.h
	g++ $(FLAGS) -c main/test.cpp -o main/test.o 
main/graphgen.o: main/graphgen.cpp ngraph/nball.cpp ngraph/nball.h ngraph/nsphere.cpp ngraph/nsphere.h ngraph/ngraph.cpp ngraph/ngraph.h
	g++ $(FLAGS) -c main/graphgen.cpp -o main/graphgen.o 
main/pgrownet.o: main/pgrownet.cpp pgrownet2d/pgrownet2d.cpp pgrownet2d/pgrownet2d.h
	g++ $(FLAGS) -c main/pgrownet.cpp -o main/pgrownet.o
main/lineppnet.o: main/lineppnet.cpp pp1d/pp1d.cpp pp1d/pp1d.h
	g++ $(FLAGS) -c main/lineppnet.cpp -o main/lineppnet.o
ngraph/nball.o: ngraph/nball.cpp ngraph/nball.h growingnetwork3d/spatialvertex.cpp growingnetwork3d/spatialvertex.h graph/growingnetwork.cpp graph/growingnetwork.h graph/graph.cpp graph/graph.h graph/vertex.cpp graph/vertex.h
	g++ $(FLAGS) -c ngraph/nball.cpp -o ngraph/nball.o
ngraph/nsphere.o: ngraph/nsphere.cpp ngraph/nsphere.h ngraph/ngraph.cpp ngraph/ngraph.h growingnetwork3d/spatialvertex.cpp growingnetwork3d/spatialvertex.h graph/growingnetwork.cpp graph/growingnetwork.h graph/graph.cpp graph/graph.h graph/vertex.cpp graph/vertex.h
	g++ $(FLAGS) -c ngraph/nsphere.cpp -o ngraph/nsphere.o
ngraph/ngraph.o: ngraph/ngraph.cpp ngraph/ngraph.h growingnetwork3d/spatialvertex.h graph/growingnetwork.h graph/graph.h graph/vertex.h
	g++ $(FLAGS) -c ngraph/ngraph.cpp -o ngraph/ngraph.o
growingnetwork3d/growingnetwork3d.o: growingnetwork3d/growingnetwork3d.cpp growingnetwork3d/growingnetwork3d.h graph/*.cpp graph/*.h
	g++ $(FLAGS) -c growingnetwork3d/growingnetwork3d.cpp -o growingnetwork3d/growingnetwork3d.o 
growingnetwork3d/spatialvertex.o: growingnetwork3d/spatialvertex.cpp growingnetwork3d/spatialvertex.h graph/vertex.cpp graph/vertex.h
	g++ $(FLAGS) -c growingnetwork3d/spatialvertex.cpp -o growingnetwork3d/spatialvertex.o 
growingnetwork2d/growingnetwork2d.o: growingnetwork2d/growingnetwork2d.cpp growingnetwork2d/growingnetwork2d.h graph/*.cpp graph/*.h
	g++ $(FLAGS) -c growingnetwork2d/growingnetwork2d.cpp -o growingnetwork2d/growingnetwork2d.o 
graph/growingnetwork.o: graph/growingnetwork.cpp graph/growingnetwork.h graph/graph.cpp graph/graph.h graph/vertex.cpp graph/vertex.h
	g++ $(FLAGS) -c graph/growingnetwork.cpp -o graph/growingnetwork.o 
graph/graph.o: graph/graph.cpp graph/graph.h graph/vertex.cpp graph/vertex.h
	g++ $(FLAGS) -c graph/graph.cpp -o graph/graph.o 
graph/vertex.o: graph/vertex.cpp graph/vertex.h
	g++ $(FLAGS) -c graph/vertex.cpp -o graph/vertex.o 
pgrownet2d/pgrownet2d.o: pgrownet2d/pgrownet2d.cpp pgrownet2d/pgrownet2d.h growingnetwork3d/spatialvertex.cpp growingnetwork3d/spatialvertex.h graph/*.cpp graph/*.h
	g++ $(FLAGS) -c pgrownet2d/pgrownet2d.cpp -o pgrownet2d/pgrownet2d.o
pp1d/pp1d.o: pp1d/pp1d.cpp pp1d/pp1d.h growingnetwork3d/spatialvertex.cpp growingnetwork3d/spatialvertex.h graph/*.cpp graph/*.h
	g++ $(FLAGS) -c pp1d/pp1d.cpp -o pp1d/pp1d.o
backup: graph/*.cpp graph/*.h growingnetwork2d/*.cpp growingnetwork2d/*.h growingnetwork3d/*.cpp growingnetwork3d/*.h pgrownet2d/*.cpp pgrownet2d/*.h ngraph/*.cpp ngraph/*.h pp1d/*.cpp pp1d/*.h main/*.cpp main/*.m makefile README doc.odt
	zip ../backup.zip graph/*.cpp graph/*.h growingnetwork2d/*.cpp growingnetwork2d/*.h growingnetwork3d/*.cpp growingnetwork3d/*.h pgrownet2d/*.cpp pgrownet2d/*.h ngraph/*.cpp ngraph/*.h pp1d/*.cpp pp1d/*.h main/*.cpp main/*.m makefile README doc.odt design.odt
clean:
	rm -vf main/pgrownet
	rm -vf main/parametertune
	rm -vf main/test
	rm -vf main/graphgen
	rm -vf main/lineppnet
	rm -vf main/figuredata
	rm -vf spatialNetworksPaper/spatialNetworksPaper.aux
	rm -vf spatialNetworksPaper/spatialNetworksPaper.bbl
	rm -vf spatialNetworksPaper/spatialNetworksPaper.blg
	rm -vf spatialNetworksPaper/spatialNetworksPaper.dvi
	rm -vf spatialNetworksPaper/spatialNetworksPaper.log
	rm -vf spatialNetworksPaper/spatialNetworksPaper.out
	rm -vf spatialNetworksPaper/spatialNetworksPaper.pdf
	rm -vf spatialNetworksPaper/spatialNetworksPaperNotes.bib
	rm -vf */*.o
	
