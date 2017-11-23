#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Triangulation_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include <vector>
#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Triangulation_geomview_ostream_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//typedef CGAL::Triangulation_3<K>      Triangulation;
typedef CGAL::Delaunay_triangulation_3<K, CGAL::Fast_location> Triangulation;

typedef Triangulation::Cell_handle    Cell_handle;
typedef Triangulation::Vertex_handle  Vertex_handle;
typedef Triangulation::Locate_type    Locate_type;
typedef Triangulation::Point          Point;


std::list<Point> L;
using std::ifstream;
using std::cout;
using std::endl;

void read_input(char* fn)
{
  ifstream ifs(fn);
  double x,y,z;
  ifs>>x>>y>>z;
  while(ifs)
    {
      L.push_front(Point(x,y,z));
      cout<<x<<" " <<y<<" "<<z<<endl;
      ifs>>x>>y>>z;

    }
}

int main()
{
  CGAL::Geomview_stream gv(CGAL::Bbox_3(5, -60, 25, 15, -45, 50));
  gv.set_line_width(4);
  // gv.set_trace(true);
  gv.set_bg_color(CGAL::Color(0, 200, 200));
  
  read_input("expt_input.txt");
  // construction from a list of points :
   // L.push_front(Point(1,1,1));
   // L.push_front(Point(1,0,0));
   // L.push_front(Point(0,1,0));
   // L.push_front(Point(0,0,1));
  // std::vector<Point> V(1);
  // V[0] = Point(0,0,1);
  // //  V[1] = Point(1,1,1);
  // //  V[2] = Point(2,2,2);

  Triangulation T(L.begin(), L.end());
  cout<<T<<endl;
  //T.insert(V.begin(), V.end());
  Triangulation::size_type n = T.number_of_vertices();
  // Triangulation::All_facets_iterator bit = T.all_facets_begin();
  // Triangulation::All_facets_iterator eit = T.all_facets_end();
  // Triangulation::All_facets_iterator it;

  // for(it = bit ; it != eit; ++bit)
  //   {
      
  //   }

  Triangulation::All_cells_iterator bit = T.all_cells_begin();
  Triangulation::All_cells_iterator eit = T.all_cells_end();
  Triangulation::All_cells_iterator it;

  for(it = bit ; it != eit; ++it)
    {
      //it->vertex(0);
      if(!T.is_infinite(it))
	{
      std::cout<<"----------Cell----------------"<<std::endl;
      for(int i = 0; i < 4; ++i)
	{
	  Triangulation::Vertex_handle vh = it->vertex(i);
	  std::cout<<*vh<<std::endl;
	}
      std::cout<<"----------Cell End----------------"<<std::endl;
	}
    }
  gv.set_wired(true);
  gv << T;
  int a;
  std::cin>>a;
  sleep(30);
  return 0;
}
