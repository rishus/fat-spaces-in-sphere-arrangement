#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include <vector>
#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Triangulation_geomview_ostream_3.h>
#include <limits>
#include <unordered_set>
#include <set>
#include <map>
#include <string>
#include <sstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned int, K> Vb;
typedef CGAL::Triangulation_data_structure_3<Vb>                       Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location> Triangulation;
typedef Triangulation::Cell_handle    Cell_handle;
typedef Triangulation::Vertex_handle  Vertex_handle;
typedef Triangulation::Facet  Facet;
typedef Triangulation::Locate_type    Locate_type;
typedef Triangulation::Point          Point;
typedef Triangulation::Vertex         Vertex;
typedef Triangulation::Cell Cell;
typedef Triangulation::Facet_iterator Facet_iterator;

using std::list;
using std::ifstream;
using std::cout;
using std::endl;
using std::ofstream;
using std::numeric_limits;
using std::vector;
using std::pair;
using std::make_pair;
using std::unordered_set;
using std::set;
using std::map;
using std::string;
using std::ostringstream;

list<pair<Point, float> > L;

struct pixel_info
{
     int is_available;
     Cell_handle ch;
     
     pixel_info()
     {
         is_available = 0;
     }

};

typedef vector<vector<vector<pixel_info> > > pixel_T;

void read_input(char* fn);
void search(double vecx[], double vecy[], double vecz[], pixel_T pixels, Triangulation T);
void case_pixel_in_cell(Point p, Cell_handle ch, pixel_T pixels, int ix, int iy, int ik);
void case_pixel_on_face(Point p, Cell_handle ch, Triangulation T);
void case_pixel_on_edge(Point p, Cell_handle ch, Triangulation T);
void case_pixel_on_vertex(Point p, Cell_handle ch, Triangulation T);
void iterate_over_faces(Triangulation& T);
int get_inf_vertex_position_in_cell(Triangulation& T, Cell_handle ch);
void get_outside_faces(Triangulation& T, vector<Facet>& all_outside_facets, set<Facet>& set_of_outside_faces);
void get_adjacent_cells(set<Facet>& facets, unordered_set<string>& cells, map<string, Cell_handle>&,  Triangulation& T);
void get_adjacent_facets(unordered_set<string>& diffused_cells, map<string, Cell_handle>& name2ch, set<Facet>& diffused_facets, Triangulation& T);
string cell2string(Cell& c);

double minx, miny, minz, maxx, maxy, maxz;
int Nx, Ny, Nz;

int main()
{
	read_input("expt_input_point.txt");
	
	cout << minx << " " << maxx << " " << miny << " " << maxy << " " << minz << " " << maxz << endl;
	Triangulation T(L.begin(), L.end());
	cout<< "outputting T: " << endl;
	cout << T <<endl;
	Triangulation::size_type n = T.number_of_vertices();
	//iterate_over_faces(T);
	 
	Locate_type lt;
	int li, lj;
	Cell_handle target_ch = T.locate(Point(0,0,0), lt, li, lj); //TODO: make it a Fe location
	Cell& target_cell = *target_ch;
	string target_cell_name = cell2string(target_cell);
	
	vector<Facet> all_outside_facets;
	set<Facet>  diffused_facets;
  	get_outside_faces(T, all_outside_facets, diffused_facets);
	for (int i=0; i < all_outside_facets.size(); i++)
	{
		Facet f = all_outside_facets[i];
		Cell c = *(f.first);
		cout << "pivot ["<<*(c.vertex(f.second))<<"] : ["<<*(c.vertex(0))<<"]  ["<<*(c.vertex(1))<<"] ["<<*(c.vertex(2))<<"]  ["<<*(c.vertex(3)) << "]" << endl;
	} 	
	unordered_set<string> diffused_cells;
	map<string, Cell_handle> name2ch;
	while(1)
	{
		if(diffused_cells.find(target_cell_name) != diffused_cells.end()) break;
		get_adjacent_cells(diffused_facets, diffused_cells, name2ch, T);
		get_adjacent_facets(diffused_cells, name2ch, diffused_facets, T);
	}
	cout<<"num diffused cells "<<diffused_cells.size()<<endl;
	////// possibly useful functions ////////
	//is_edge(cell_handle, int i, int j)/////
	//is_facet(u,v,w,&cell_handle, &i,&j,&k) ////
	//facet_circulator : all_facets_begin()///
	// build the search box

	return 0;

}

void read_input(char* fn)
{
	ifstream ifs(fn);
	double x,y,z;
	minx=numeric_limits<double>::max() , miny=numeric_limits<double>::max(), minz=numeric_limits<double>::max();
	maxx=numeric_limits<double>::min(), maxy=numeric_limits<double>::min(),maxz=numeric_limits<double>::min();
	ifs >> x >> y >> z;
	while(ifs)
	{
		L.push_front(make_pair(Point(x,y,z), 1.0));
		cout << x << " " << y << " " << z << endl;
		ifs >> x >> y >> z;
		if (x - 1.0 < minx)
		{
			minx = x - 1.0;
		}
		if (x + 1.0 > maxx)
		{
			maxx = x + 1.0;
		}
		if (y - 1.0 < miny)
		{
			miny = y - 1.0;
		}
		if (y + 1.0 > maxy)
		{
			maxy = y + 1.0;
		}
		if (z - 1.0 < minz)
		{
			minz = z - 1.0;
		}
		if (z + 1.0 > maxz)
		{
			maxz = z + 1.0;
		}

	}

}

int get_inf_vertex_position_in_cell(Triangulation& T, Cell_handle ch)
{

	for(int i = 0; i < 4; ++i)
		if(T.is_infinite(ch->vertex(i))) return i;
	assert(0);
	return -1;

}

void get_outside_faces(Triangulation& T, vector<Facet>& all_outside_facets, set<Facet>& set_of_outside_faces)
{
	Vertex_handle infv = T.infinite_vertex();
	Triangulation::All_cells_iterator bit = T.all_cells_begin();
	Triangulation::All_cells_iterator eit = T.all_cells_end();
	Triangulation::All_cells_iterator it;
	int count = 1;
	for(it = bit ; it != eit; ++it)
	{
		if(T.is_infinite(it))
		{
			Cell_handle ch = it;
			//Facet outside_face(make_pair<Cell_handle, int>(ch, get_inf_vertex_position_in_cell(T,ch)));
			Facet outside_face(make_pair(ch, get_inf_vertex_position_in_cell(T,ch)));
			all_outside_facets.push_back(outside_face);
			set_of_outside_faces.insert(outside_face);
			//Facet f(it,0);
			//Facet mf = T.mirror_facet(f);
			//Cell_handle adj_cell_handle = mf.first;
			//Vertex_handle u = it->vertex(0);
			//Face_handle f = it->face(0);

			//for(int i = 0; i < 4; ++i)
			//{
			//	Vertex_handle vh = it->vertex(i);
			//	cout<<*vh<<endl;
			//}

		}
		count++;
	}

}

void iterate_over_faces(Triangulation& T)
{
	Facet_iterator fit = T.all_facets_begin();
	Facet_iterator eit = T.all_facets_end();
	for (; fit!=eit; ++fit)
	{
		Facet f = *fit;
		Cell_handle ch = f.first;
		Cell c = *ch;
		Facet mf = T.mirror_facet(f);
		Cell_handle mch = mf.first;
		Cell mc = *mch;
		cout << "pivot ["<<*(c.vertex(f.second))<<"] : ["<<*(c.vertex(0))<<"]  ["<<*(c.vertex(1))<<"] ["<<*(c.vertex(2))<<"]  ["<<*(c.vertex(3)) << endl;
		cout << "pivot ["<<*(mc.vertex(mf.second))<<"] :["<<*(mc.vertex(0))<<"] ["<<*(mc.vertex(1))<<"] ["<<*(mc.vertex(2))<<"] ["<<*(mc.vertex(3)) << endl;
		cout<<"--------"<<endl;
	}

}

void explore_path()
{
	//for (Delaunay::Finite_faces_iterator fit = triangulation.finite_faces_begin(); 
	//		fit != triangulation.finite_faces_end(); ++fit)
	//{
		

	//}

}

void get_adjacent_cells(set<Facet>& facets, unordered_set<string>& cells, map<string,Cell_handle>& name2ch, Triangulation& T)
{
	set<Facet>::iterator setIter;
	for (setIter = facets.begin(); setIter!=facets.end(); ++setIter)
	{
		Facet f = *setIter;
		Facet mf = T.mirror_facet(f);
		Cell_handle adj_cell_handle = mf.first;
		Cell adj_cell = *adj_cell_handle;
		string adj_cell_name = cell2string(adj_cell);
		cells.insert(adj_cell_name);
		if(name2ch.find(adj_cell_name) == name2ch.end())
			name2ch[adj_cell_name] = adj_cell_handle;

	}
	
}

void get_bounding_faces(unordered_set<string>& cells, set<Facet>& facets, Triangulation& T)
{
	unordered_set<string>::iterator sit;
	for(sit = cells.begin(); sit != cells.end(); ++sit)
	{
		Cell_handle ch; //get the corresponding cell
		for(int i = 0; i < 4; ++i)
		{
				Facet f(make_pair(ch, i));
				if(!T.is_infinite(f))
					facets.insert(f);
		}
	}	

}

void get_adjacent_facets(unordered_set<string>& diffused_cells, map<string, Cell_handle>& name2ch, set<Facet>& diffused_facets, Triangulation& T)
{
	unordered_set<string>::iterator sit;
	for(sit = diffused_cells.begin(); sit != diffused_cells.end(); ++sit)
	{
		Cell_handle ch = name2ch[*sit];
		for(int i = 0; i < 4; ++i)
		{
				Facet f(make_pair(ch, i));
				if(!T.is_infinite(f))
					diffused_facets.insert(f);
		}

	}
		

}
string cell2string(Cell& c)
{
	ostringstream cstr;
	for(int i = 0; i < 4; ++i)
	{	
		cstr<<*(c.vertex(i))<<":";
	}
	
	return cstr.str();

}

