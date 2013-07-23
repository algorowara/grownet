Flat profile:

Each sample counts as 0.01 seconds.
  %   cumulative   self              self     total           
 time   seconds   seconds    calls   s/call   s/call  name    
 92.31    386.62   386.62    18742     0.02     0.02  NBall::sumForces(SpatialVertex*)
  7.42    417.68    31.06 3656201218     0.00     0.00  NGraph::getNode(long) const
  0.22    418.59     0.91    20795     0.00     0.00  NBall::randomLocation()
  0.04    418.76     0.17  5786357     0.00     0.00  SpatialVertex::radialDistance()
  0.01    418.80     0.04      995     0.00     0.00  NBall::findMNearestNeighbors(SpatialVertex*)
  0.01    418.83     0.03      996     0.00     0.42  NBall::gradientDescent(float, float, long)
  0.00    418.83     0.00     4207     0.00     0.00  std::vector<Vertex*, std::allocator<Vertex*> >::_M_insert_aux(__gnu_cxx::__normal_iterator<Vertex**, std::vector<Vertex*, std::allocator<Vertex*> > >, Vertex* const&)
  0.00    418.83     0.00     3990     0.00     0.00  Vertex::addNeighbor(Vertex*)
  0.00    418.83     0.00     1000     0.00     0.00  SpatialVertex::SpatialVertex(long, float*, long)
  0.00    418.83     0.00     1000     0.00     0.00  GrowingNetwork::tick()
  0.00    418.83     0.00     1000     0.00     0.00  Graph::addNode(Vertex*)
  0.00    418.83     0.00     1000     0.00     0.00  Vertex::Vertex()
  0.00    418.83     0.00     1000     0.00     0.00  GrowingNetwork::getTime() const
  0.00    418.83     0.00      992     0.00     0.00  Vertex::clusteringCoefficient()
  0.00    418.83     0.00        1     0.00     0.00  _GLOBAL__sub_I__ZN16GrowingNetwork2DC2Ell
  0.00    418.83     0.00        1     0.00     0.00  _GLOBAL__sub_I__ZN16GrowingNetwork3DC2Ev
  0.00    418.83     0.00        1     0.00     0.00  _GLOBAL__sub_I__ZN30PositiveChargeGrowingNetwork2DC2Ellffl
  0.00    418.83     0.00        1     0.00     0.00  _GLOBAL__sub_I__ZN5Graph7addNodeEP6Vertex
  0.00    418.83     0.00        1     0.00     0.00  _GLOBAL__sub_I__ZN5NBallC2Elliffflll
  0.00    418.83     0.00        1     0.00     0.00  _GLOBAL__sub_I__ZN6VertexC2Ev
  0.00    418.83     0.00        1     0.00     0.00  _GLOBAL__sub_I__ZN7NSphereC2Elliffflll
  0.00    418.83     0.00        1     0.00     0.00  _GLOBAL__sub_I__ZNK14GrowingNetwork7getTimeEv
  0.00    418.83     0.00        1     0.00     0.00  _GLOBAL__sub_I_main
  0.00    418.83     0.00        1     0.00     0.00  GrowingNetwork::nodeAgeVsDegree()
  0.00    418.83     0.00        1     0.00   418.41  NBall::grow(long)
  0.00    418.83     0.00        1     0.00     0.00  NGraph::NGraph(long)

 %         the percentage of the total running time of the
time       program used by this function.

cumulative a running sum of the number of seconds accounted
 seconds   for by this function and those listed above it.

 self      the number of seconds accounted for by this
seconds    function alone.  This is the major sort for this
           listing.

calls      the number of times this function was invoked, if
           this function is profiled, else blank.
 
 self      the average number of milliseconds spent in this
ms/call    function per call, if this function is profiled,
	   else blank.

 total     the average number of milliseconds spent in this
ms/call    function and its descendents per call, if this 
	   function is profiled, else blank.

name       the name of the function.  This is the minor sort
           for this listing. The index shows the location of
	   the function in the gprof listing. If the index is
	   in parenthesis it shows where it would appear in
	   the gprof listing if it were to be printed.

		     Call graph (explanation follows)


granularity: each sample hit covers 4 byte(s) for 0.00% of 418.83 seconds

index % time    self  children    called     name
                                                 <spontaneous>
[1]    100.0    0.00  418.83                 NBall::NBall(long, long, int, float, float, float, long, long, long) [1]
                0.00  418.41       1/1           NBall::grow(long) [3]
                0.00    0.42       1/996         NBall::gradientDescent(float, float, long) [2]
                0.00    0.00       8/20795       NBall::randomLocation() [6]
                0.00    0.00      10/3656201218     NGraph::getNode(long) const [5]
                0.00    0.00      10/3990        Vertex::addNeighbor(Vertex*) [15]
                0.00    0.00       5/1000        GrowingNetwork::getTime() const [20]
                0.00    0.00       5/1000        SpatialVertex::SpatialVertex(long, float*, long) [16]
                0.00    0.00       5/1000        Graph::addNode(Vertex*) [18]
                0.00    0.00       5/1000        GrowingNetwork::tick() [17]
                0.00    0.00       1/1           NGraph::NGraph(long) [32]
-----------------------------------------------
                0.00    0.42       1/996         NBall::NBall(long, long, int, float, float, float, long, long, long) [1]
                0.03  418.26     995/996         NBall::grow(long) [3]
[2]    100.0    0.03  418.68     996         NBall::gradientDescent(float, float, long) [2]
              386.62   31.21   18742/18742       NBall::sumForces(SpatialVertex*) [4]
                0.85    0.00   19507/20795       NBall::randomLocation() [6]
-----------------------------------------------
                0.00  418.41       1/1           NBall::NBall(long, long, int, float, float, float, long, long, long) [1]
[3]     99.9    0.00  418.41       1         NBall::grow(long) [3]
                0.03  418.26     995/996         NBall::gradientDescent(float, float, long) [2]
                0.04    0.02     995/995         NBall::findMNearestNeighbors(SpatialVertex*) [8]
                0.06    0.00    1280/20795       NBall::randomLocation() [6]
                0.00    0.00    7960/5786357     SpatialVertex::radialDistance() [7]
                0.00    0.00    3980/3990        Vertex::addNeighbor(Vertex*) [15]
                0.00    0.00     995/1000        GrowingNetwork::getTime() const [20]
                0.00    0.00     995/1000        SpatialVertex::SpatialVertex(long, float*, long) [16]
                0.00    0.00     995/1000        Graph::addNode(Vertex*) [18]
                0.00    0.00     995/1000        GrowingNetwork::tick() [17]
-----------------------------------------------
                             6478471             NBall::sumForces(SpatialVertex*) [4]
              386.62   31.21   18742/18742       NBall::gradientDescent(float, float, long) [2]
[4]     99.8  386.62   31.21   18742+6478471 NBall::sumForces(SpatialVertex*) [4]
               31.04    0.00 3653682737/3656201218     NGraph::getNode(long) const [5]
                0.17    0.00 5778397/5786357     SpatialVertex::radialDistance() [7]
                             6478471             NBall::sumForces(SpatialVertex*) [4]
-----------------------------------------------
                0.00    0.00      10/3656201218     NBall::NBall(long, long, int, float, float, float, long, long, long) [1]
                0.02    0.00 2518471/3656201218     NBall::findMNearestNeighbors(SpatialVertex*) [8]
               31.04    0.00 3653682737/3656201218     NBall::sumForces(SpatialVertex*) [4]
[5]      7.4   31.06    0.00 3656201218         NGraph::getNode(long) const [5]
-----------------------------------------------
                0.00    0.00       8/20795       NBall::NBall(long, long, int, float, float, float, long, long, long) [1]
                0.06    0.00    1280/20795       NBall::grow(long) [3]
                0.85    0.00   19507/20795       NBall::gradientDescent(float, float, long) [2]
[6]      0.2    0.91    0.00   20795         NBall::randomLocation() [6]
-----------------------------------------------
                0.00    0.00    7960/5786357     NBall::grow(long) [3]
                0.17    0.00 5778397/5786357     NBall::sumForces(SpatialVertex*) [4]
[7]      0.0    0.17    0.00 5786357         SpatialVertex::radialDistance() [7]
-----------------------------------------------
                0.04    0.02     995/995         NBall::grow(long) [3]
[8]      0.0    0.04    0.02     995         NBall::findMNearestNeighbors(SpatialVertex*) [8]
                0.02    0.00 2518471/3656201218     NGraph::getNode(long) const [5]
-----------------------------------------------
                0.00    0.00      11/4207        Graph::addNode(Vertex*) [18]
                0.00    0.00    4196/4207        Vertex::addNeighbor(Vertex*) [15]
[14]     0.0    0.00    0.00    4207         std::vector<Vertex*, std::allocator<Vertex*> >::_M_insert_aux(__gnu_cxx::__normal_iterator<Vertex**, std::vector<Vertex*, std::allocator<Vertex*> > >, Vertex* const&) [14]
-----------------------------------------------
                0.00    0.00      10/3990        NBall::NBall(long, long, int, float, float, float, long, long, long) [1]
                0.00    0.00    3980/3990        NBall::grow(long) [3]
[15]     0.0    0.00    0.00    3990         Vertex::addNeighbor(Vertex*) [15]
                0.00    0.00    4196/4207        std::vector<Vertex*, std::allocator<Vertex*> >::_M_insert_aux(__gnu_cxx::__normal_iterator<Vertex**, std::vector<Vertex*, std::allocator<Vertex*> > >, Vertex* const&) [14]
-----------------------------------------------
                0.00    0.00       5/1000        NBall::NBall(long, long, int, float, float, float, long, long, long) [1]
                0.00    0.00     995/1000        NBall::grow(long) [3]
[16]     0.0    0.00    0.00    1000         SpatialVertex::SpatialVertex(long, float*, long) [16]
                0.00    0.00    1000/1000        Vertex::Vertex() [19]
-----------------------------------------------
                0.00    0.00       5/1000        NBall::NBall(long, long, int, float, float, float, long, long, long) [1]
                0.00    0.00     995/1000        NBall::grow(long) [3]
[17]     0.0    0.00    0.00    1000         GrowingNetwork::tick() [17]
-----------------------------------------------
                0.00    0.00       5/1000        NBall::NBall(long, long, int, float, float, float, long, long, long) [1]
                0.00    0.00     995/1000        NBall::grow(long) [3]
[18]     0.0    0.00    0.00    1000         Graph::addNode(Vertex*) [18]
                0.00    0.00      11/4207        std::vector<Vertex*, std::allocator<Vertex*> >::_M_insert_aux(__gnu_cxx::__normal_iterator<Vertex**, std::vector<Vertex*, std::allocator<Vertex*> > >, Vertex* const&) [14]
-----------------------------------------------
                0.00    0.00    1000/1000        SpatialVertex::SpatialVertex(long, float*, long) [16]
[19]     0.0    0.00    0.00    1000         Vertex::Vertex() [19]
-----------------------------------------------
                0.00    0.00       5/1000        NBall::NBall(long, long, int, float, float, float, long, long, long) [1]
                0.00    0.00     995/1000        NBall::grow(long) [3]
[20]     0.0    0.00    0.00    1000         GrowingNetwork::getTime() const [20]
-----------------------------------------------
                0.00    0.00     992/992         GrowingNetwork::nodeAgeVsDegree() [31]
[21]     0.0    0.00    0.00     992         Vertex::clusteringCoefficient() [21]
-----------------------------------------------
                0.00    0.00       1/1           __libc_csu_init [118]
[22]     0.0    0.00    0.00       1         _GLOBAL__sub_I__ZN16GrowingNetwork2DC2Ell [22]
-----------------------------------------------
                0.00    0.00       1/1           __libc_csu_init [118]
[23]     0.0    0.00    0.00       1         _GLOBAL__sub_I__ZN16GrowingNetwork3DC2Ev [23]
-----------------------------------------------
                0.00    0.00       1/1           __libc_csu_init [118]
[24]     0.0    0.00    0.00       1         _GLOBAL__sub_I__ZN30PositiveChargeGrowingNetwork2DC2Ellffl [24]
-----------------------------------------------
                0.00    0.00       1/1           __libc_csu_init [118]
[25]     0.0    0.00    0.00       1         _GLOBAL__sub_I__ZN5Graph7addNodeEP6Vertex [25]
-----------------------------------------------
                0.00    0.00       1/1           __libc_csu_init [118]
[26]     0.0    0.00    0.00       1         _GLOBAL__sub_I__ZN5NBallC2Elliffflll [26]
-----------------------------------------------
                0.00    0.00       1/1           __libc_csu_init [118]
[27]     0.0    0.00    0.00       1         _GLOBAL__sub_I__ZN6VertexC2Ev [27]
-----------------------------------------------
                0.00    0.00       1/1           __libc_csu_init [118]
[28]     0.0    0.00    0.00       1         _GLOBAL__sub_I__ZN7NSphereC2Elliffflll [28]
-----------------------------------------------
                0.00    0.00       1/1           __libc_csu_init [118]
[29]     0.0    0.00    0.00       1         _GLOBAL__sub_I__ZNK14GrowingNetwork7getTimeEv [29]
-----------------------------------------------
                0.00    0.00       1/1           __libc_csu_init [118]
[30]     0.0    0.00    0.00       1         _GLOBAL__sub_I_main [30]
-----------------------------------------------
                0.00    0.00       1/1           Graph::clustering() const [102]
[31]     0.0    0.00    0.00       1         GrowingNetwork::nodeAgeVsDegree() [31]
                0.00    0.00     992/992         Vertex::clusteringCoefficient() [21]
-----------------------------------------------
                0.00    0.00       1/1           NBall::NBall(long, long, int, float, float, float, long, long, long) [1]
[32]     0.0    0.00    0.00       1         NGraph::NGraph(long) [32]
-----------------------------------------------

 This table describes the call tree of the program, and was sorted by
 the total amount of time spent in each function and its children.

 Each entry in this table consists of several lines.  The line with the
 index number at the left hand margin lists the current function.
 The lines above it list the functions that called this function,
 and the lines below it list the functions this one called.
 This line lists:
     index	A unique number given to each element of the table.
		Index numbers are sorted numerically.
		The index number is printed next to every function name so
		it is easier to look up where the function in the table.

     % time	This is the percentage of the `total' time that was spent
		in this function and its children.  Note that due to
		different viewpoints, functions excluded by options, etc,
		these numbers will NOT add up to 100%.

     self	This is the total amount of time spent in this function.

     children	This is the total amount of time propagated into this
		function by its children.

     called	This is the number of times the function was called.
		If the function called itself recursively, the number
		only includes non-recursive calls, and is followed by
		a `+' and the number of recursive calls.

     name	The name of the current function.  The index number is
		printed after it.  If the function is a member of a
		cycle, the cycle number is printed between the
		function's name and the index number.


 For the function's parents, the fields have the following meanings:

     self	This is the amount of time that was propagated directly
		from the function into this parent.

     children	This is the amount of time that was propagated from
		the function's children into this parent.

     called	This is the number of times this parent called the
		function `/' the total number of times the function
		was called.  Recursive calls to the function are not
		included in the number after the `/'.

     name	This is the name of the parent.  The parent's index
		number is printed after it.  If the parent is a
		member of a cycle, the cycle number is printed between
		the name and the index number.

 If the parents of the function cannot be determined, the word
 `<spontaneous>' is printed in the `name' field, and all the other
 fields are blank.

 For the function's children, the fields have the following meanings:

     self	This is the amount of time that was propagated directly
		from the child into the function.

     children	This is the amount of time that was propagated from the
		child's children to the function.

     called	This is the number of times the function called
		this child `/' the total number of times the child
		was called.  Recursive calls by the child are not
		listed in the number after the `/'.

     name	This is the name of the child.  The child's index
		number is printed after it.  If the child is a
		member of a cycle, the cycle number is printed
		between the name and the index number.

 If there are any cycles (circles) in the call graph, there is an
 entry for the cycle-as-a-whole.  This entry shows who called the
 cycle (as parents) and the members of the cycle (as children.)
 The `+' recursive calls entry shows the number of function calls that
 were internal to the cycle, and the calls entry for each member shows,
 for that member, how many times it was called from other members of
 the cycle.


Index by function name

  [22] _GLOBAL__sub_I__ZN16GrowingNetwork2DC2Ell [7] SpatialVertex::radialDistance() [4] NBall::sumForces(SpatialVertex*)
  [23] _GLOBAL__sub_I__ZN16GrowingNetwork3DC2Ev [16] SpatialVertex::SpatialVertex(long, float*, long) [32] NGraph::NGraph(long)
  [24] _GLOBAL__sub_I__ZN30PositiveChargeGrowingNetwork2DC2Ellffl [31] GrowingNetwork::nodeAgeVsDegree() [15] Vertex::addNeighbor(Vertex*)
  [25] _GLOBAL__sub_I__ZN5Graph7addNodeEP6Vertex [17] GrowingNetwork::tick() [21] Vertex::clusteringCoefficient()
  [26] _GLOBAL__sub_I__ZN5NBallC2Elliffflll [18] Graph::addNode(Vertex*) [19] Vertex::Vertex()
  [27] _GLOBAL__sub_I__ZN6VertexC2Ev [6] NBall::randomLocation() [20] GrowingNetwork::getTime() const
  [28] _GLOBAL__sub_I__ZN7NSphereC2Elliffflll [2] NBall::gradientDescent(float, float, long) [5] NGraph::getNode(long) const
  [29] _GLOBAL__sub_I__ZNK14GrowingNetwork7getTimeEv [8] NBall::findMNearestNeighbors(SpatialVertex*) [14] std::vector<Vertex*, std::allocator<Vertex*> >::_M_insert_aux(__gnu_cxx::__normal_iterator<Vertex**, std::vector<Vertex*, std::allocator<Vertex*> > >, Vertex* const&)
  [30] _GLOBAL__sub_I_main     [3] NBall::grow(long)
