# include <cstddef>
# include <ctime>
# include <iostream> 

using namespace std;

# include <metis.h>

int main ( );
void partgraphrecursive_test ( );
void partgraphkway_test ( );
void timestamp ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    METIS_TEST tests the METIS library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 February 2017
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "METIS_TEST\n";
  cout << "  C++ version\n";
  cout << "  Test the METIS library for graph partitioning.\n";

  partgraphrecursive_test ( );
  partgraphkway_test ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "METIS_TEST\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void partgraphrecursive_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    PARTGRAPHRECURSIVE_TEST tests PARTGRAPHRECURSIVE.
//
//  Discussion:
//
//    The graph has the following form:
//
//      0 --- 1 --- 2
//      |     |     |
//      3 --- 4 --- 5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 February 2017
//
//  Author:
//
//    John Burkardt
//
{
//
//  The number of vertices.
//
  idx_t nvtxs = 6;
//
// Number of balancing constraints, which must be at least 1.
//
  idx_t ncon = 1;
//
//  Pointers to initial entries in adjacency array.
//
  idx_t xadj[nvtxs+1] = { 0, 2, 5, 7, 9, 12, 14 };
//
// Adjacent vertices in consecutive index order.
//
  idx_t nEdges = 7;
  idx_t adjncy[2 * nEdges] = {1,3,0,4,2,1,5,0,4,3,1,5,4,2};
//
//  The number of parts requested for the partition.
//
  idx_t nParts = 2;
//
//  On return, the edge cut volume of the partitioning solution.
//
  idx_t objval;
//
//  On return, the partition vector for the graph.
//
  idx_t part[nvtxs];

  cout << "\n";
  cout << "PARTGRAPHRECURSIVE_TEST:\n";
  cout << "  METIS_PartGraphRecursive partitions a graph into K parts\n";
  cout << "  using multilevel recursive bisection.\n";

  int ret = METIS_PartGraphRecursive ( &nvtxs, &ncon, xadj, adjncy, NULL, NULL, 
    NULL, &nParts, NULL, NULL, NULL, &objval, part );

  cout << "\n";
  cout << "  Return code = " << ret << "\n";
  cout << "  Edge cuts for partition = " << objval << "\n";
    
  cout << "\n";
  cout << "  Partition vector:\n";
  cout << "\n";
  cout << "  Node  Part\n";
  cout << "\n";
  for ( unsigned part_i = 0; part_i < nvtxs; part_i++ )
  {
	std::cout << "     " << part_i << "     " << part[part_i] << std::endl;
  }
  
  return;
}
//****************************************************************************80

void partgraphkway_test ( )

//****************************************************************************80
//  Purpose:
//
//    PARTGRAPHKWAY_TEST tests PARTGRAPHKWAY.
//
//  Discussion:
//
//    The graph has the following form:
//
//      0 --- 1 --- 2
//      |     |     |
//      3 --- 4 --- 5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 February 2017
//
//  Author:
//
//    John Burkardt
//
{
//
//  The number of vertices.
//
  idx_t nvtxs = 6;
//
// Number of balancing constraints, which must be at least 1.
//
  idx_t ncon = 1;
//
//  Pointers to initial entries in adjacency array.
//
  idx_t xadj[nvtxs+1] = { 0, 2, 5, 7, 9, 12, 14 };
//
// Adjacent vertices in consecutive index order.
//
  idx_t nEdges = 7;
  idx_t adjncy[2 * nEdges] = {1,3,0,4,2,1,5,0,4,3,1,5,4,2};
//
//  The number of parts requested for the partition.
//
  idx_t nParts = 2;
//
//  On return, the edge cut volume of the partitioning solution.
//
  idx_t objval;
//
//  On return, the partition vector for the graph.
//
  idx_t part[nvtxs];

  cout << "\n";
  cout << "PARTGRAPHKWAY_TEST:\n";
  cout << "  METIS_PartGraphKway partitions a graph into K parts\n";
  cout << "  using multilevel K-way partition.\n";

  int ret = METIS_PartGraphKway ( &nvtxs, &ncon, xadj, adjncy, NULL, NULL, 
    NULL, &nParts, NULL, NULL, NULL, &objval, part );

  cout << "\n";
  cout << "  Return code = " << ret << "\n";
  cout << "  Edge cuts for partition = " << objval << "\n";
    
  cout << "\n";
  cout << "  Partition vector:\n";
  cout << "\n";
  cout << "  Node  Part\n";
  cout << "\n";
  for ( unsigned part_i = 0; part_i < nvtxs; part_i++ )
  {
	std::cout << "     " << part_i << "     " << part[part_i] << std::endl;
  }
    
  return;
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
