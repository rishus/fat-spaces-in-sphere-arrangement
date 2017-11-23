#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include <vector>
#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Triangulation_geomview_ostream_3.h>
#include <limits>
#include <map>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_3<K, CGAL::Fast_location> Triangulation;
typedef Triangulation::Cell_handle    Cell_handle;
typedef Triangulation::Vertex_handle  Vertex_handle;
typedef Triangulation::Locate_type    Locate_type;
typedef Triangulation::Point          Point;
typedef Triangulation::Vertex         Vertex;

std::list<Point> L;
using std::ifstream;
using std::cout;
using std::endl;
using std::ofstream;
using std::numeric_limits;
using std::vector;
using std::map;

struct pixel_info
{
     int is_cavity;
     Cell_handle ch;
     
     pixel_info()
     {
         is_cavity = 0;
     }

};

typedef vector<vector<vector<pixel_info> > > pixel_T;
typedef map<Vertex,float> atom_radius_T;
void read_input(char* fn);
void search(double vecx[], double vecy[], double vecz[], pixel_T pixels, Triangulation T);
void case_pixel_in_cell(Point p, Cell_handle ch, Triangulation T);
void case_pixel_on_face(Point p, Cell_handle ch, Triangulation T);
void case_pixel_on_edge(Point p, Cell_handle ch, Triangulation T);
void case_pixel_on_vertex(Point p, Cell_handle ch, Triangulation T);

double minx, miny, minz, maxx, maxy, maxz;
int Nx, Ny, Nz;

int main()
{
	read_input("expt_input_point.txt");

	Triangulation T(L.begin(), L.end());
	cout<< "outputting T: " << endl;
	cout << T <<endl;
	Triangulation::size_type n = T.number_of_vertices();


	Triangulation::All_cells_iterator bit = T.all_cells_begin();
	Triangulation::All_cells_iterator eit = T.all_cells_end();
	Triangulation::All_cells_iterator it;
	int count = 1;
	////// possibly useful functions ////////
	//is_edge(cell_handle, int i, int j)/////
	//is_facet(u,v,w,&cell_handle, &i,&j,&k) ////
	//facet_circulator : all_facets_begin()///
	for(it = bit ; it != eit; ++it)
	{
		if(!T.is_infinite(it))
		{
			//Cell_handle ch = *it;
			Vertex_handle u = it->vertex(0);

			cout<<"----------Cell--" << count << "--------------"<<endl;
			for(int i = 0; i < 4; ++i)
			{
				Vertex_handle vh = it->vertex(i);
				cout<<*vh<<endl;
			}
			cout<<"----------Cell End----------------"<<endl;
		}
		count++;
	}

	// build the search box
	double dx = 0.1;
	Nx = (maxx-minx)/dx; 
	Ny = (maxy-miny)/dx; 
	Nz=(maxz - minz)/dx;
	cout << Nx << Ny << Nz << endl;
	double vecx[Nx], vecy[Ny], vecz[Nz];
	pixel_T  pixels;
	pixels.resize(Nx);
	for (int i = 0; i < Nx; i++)
	{
		pixels[i].resize(Ny);
		for (int j=0; j<Ny; j++)
		{
		    pixels[i][j].resize(Nz);
		}
	}

	for (int ix = 0; ix < Nx; ix++){
		vecx[ix] = minx + ix * dx;
	}
	for (int iy = 0; iy < Ny; iy++){
		vecy[iy] = miny + iy * dx;
	}
	for (int iz = 0; iz < Nz; iz++){
		vecz[iz] = (minz + iz*dx);
	}
			
	// now do the search
	cout << "Now going to search " << endl;
	search(vecx, vecy, vecz, pixels, T);

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
		L.push_front(Point(x,y,z));
		cout << x << " " << y << " " << z << endl;
		ifs >> x >> y >> z;
		if (x < minx)
		{
			minx = x;
		}
		if (x > maxx)
		{
			maxx = x;
		}
		if (y < miny)
		{
			miny = y;
		}
		if (y > maxy)
		{
			maxy = y;
		}
		if (z < minz)
		{
			minz = z;
		}
		if (z > maxz)
		{
			maxz = z;
		}

	}

}

void search(double vecx[], double vecy[], double vecz[], pixel_T pixels, Triangulation T)
{

	Locate_type lt;
	int li, lj;
	for (int ix=0; ix<Nx; ix++)
	{
		for (int iy=0; iy<Ny; iy++)
		{
			for (int iz=0; iz<Nz; iz++)
			{
				Point p(vecx[ix], vecy[iy], vecz[iz]);  //p(0.9, 0.9, 0.5);  //
				// locating the cell this point belongs to
				Cell_handle ch = T.locate(p, lt, li, lj);
				if (lt == Triangulation::CELL)
				{
					cout << "Point (" <<  p << ") belongs to Cell " << 
						*(ch->vertex(0)) << ", " << *(ch->vertex(1)) << ", " << *(ch->vertex(2)) << ", " << *(ch->vertex(3)) 
						<<  endl;
					cout << "Point individual coordinates are " << p.x() << ", " << p.y() << ", " << p.z() << endl;
					case_pixel_in_cell(p, ch, T);
				}
				else if (lt == Triangulation::FACET)
				{
					cout << "Point " << p << "is on Facet " << li <<  "of cell " << 
						*(ch->vertex(0)) << ", " << *(ch->vertex(1)) << ", " << *(ch->vertex(2)) << ", " << *(ch->vertex(3)) 
						<< endl;
				}
				else if (lt == Triangulation::EDGE)
				{
					cout << "This point is on Edge " << li << " , " << lj << " of cell "  <<
						*(ch->vertex(0)) << ", " << *(ch->vertex(1)) << ", " << *(ch->vertex(2)) << ", " << *(ch->vertex(3)) 
						<< endl;
				}
				else if (lt == Triangulation::VERTEX)
				{
					cout << "This point is Vertex " << li << "of cell " << 
						*(ch->vertex(0)) << ", " << *(ch->vertex(1)) << ", " << *(ch->vertex(2)) << ", " << *(ch->vertex(3))
						<< endl;
					assert(ch->vertex(li)->point() == p);
				}
				else if (lt == Triangulation::OUTSIDE_CONVEX_HULL)
				{
					cout << "Point lies outside convex hull " << endl;
				}
			}
		}
	}

}

void case_pixel_in_cell(Point p, Cell_handle ch, Triangulation T)
{
	
	for (int i=0; i<4; i++)
	{
		float dist0 = (p.x() - (ch->vertex(i)->point().x())) * (p.x() - (ch->vertex(i)->point().x())) 
			+ (p.y() - (ch->vertex(i)->point().y())) * (p.y() - (ch->vertex(i)->point().y()))  
			+ (p.z() - (ch->vertex(i)->point().z())) * (p.z() - (ch->vertex(i)->point().z()));
		dist0 = sqrt(dist0);
		cout << "distance from vertex " << i << " = " << dist0 << endl;
		/**float dist1 = (p.x() - (ch->vertex(1)->point().x())) * (p.x() - (ch->vertex(1)->point().x())) 
			+ (p.y() - (ch->vertex(1)->point().y())) * (p.y() - (ch->vertex(1)->point().y()))  
			+ (p.z() - (ch->vertex(1)->point().z())) * (p.z() - (ch->vertex(1)->point().z()));
		dist1 = sqrt(dist1);
		cout << "distance from vertex 1 = " << dist1 << endl;
		float dist2 = (p.x() - (ch->vertex(2)->point().x())) * (p.x() - (ch->vertex(2)->point().x())) 
			+ (p.y() - (ch->vertex(2)->point().y())) * (p.y() - (ch->vertex(2)->point().y()))  
			+ (p.z() - (ch->vertex(2)->point().z())) * (p.z() - (ch->vertex(2)->point().z()));
		dist2 = sqrt(dist1);
		cout << "distance from vertex 2 = " << dist2 << endl;
		float dist3 = (p.x() - (ch->vertex(3)->point().x())) * (p.x() - (ch->vertex(3)->point().x())) 
			+ (p.y() - (ch->vertex(3)->point().y())) * (p.y() - (ch->vertex(3)->point().y()))  
			+ (p.z() - (ch->vertex(3)->point().z())) * (p.z() - (ch->vertex(3)->point().z()));
		dist3 = sqrt(dist3);
		cout << "distance from vertex 3 = " << dist3 << endl;*/

	}


}

void case_pixel_on_face(Point p, Cell_handle ch, Triangulation T)
{

}

void case_pixel_on_edge(Point p, Cell_handle ch, Triangulation T)
{

}

void case_pixel_on_vertex(Point p, Cell_handle ch, Triangulation T)
{

}

