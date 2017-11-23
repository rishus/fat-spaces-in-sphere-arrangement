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

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned int, K> Vb;
typedef CGAL::Triangulation_data_structure_3<Vb>                       Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location> Triangulation;
typedef Triangulation::Cell_handle    Cell_handle;
typedef Triangulation::Vertex_handle  Vertex_handle;
typedef Triangulation::Locate_type    Locate_type;
typedef Triangulation::Point          Point;
typedef Triangulation::Vertex         Vertex;

using std::list;
using std::ifstream;
using std::cout;
using std::endl;
using std::ofstream;
using std::numeric_limits;
using std::vector;
using std::pair;
using std::make_pair;


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
	double dx = 1.4;
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
			/*		cout << "Point (" <<  p << ") belongs to Cell " << 
						*(ch->vertex(0)) << ", " << *(ch->vertex(1)) << ", " << *(ch->vertex(2)) << ", " << *(ch->vertex(3)) 
						<<  endl; */
					case_pixel_in_cell(p, ch, pixels, ix, iy, iz);
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

void case_pixel_in_cell(Point p, Cell_handle ch, pixel_T pixels, int ix, int iy, int iz)
{
	for (int i=0; i<4; i++)
	{
		float radiusi = ch->vertex(i)->info();
		float disti = (p.x() - (ch->vertex(i)->point().x())) * (p.x() - (ch->vertex(i)->point().x())) 
			+ (p.y() - (ch->vertex(i)->point().y())) * (p.y() - (ch->vertex(i)->point().y()))  
			+ (p.z() - (ch->vertex(i)->point().z())) * (p.z() - (ch->vertex(i)->point().z()));
		disti = sqrt(disti);
		cout << "distance from vertex " << i << " = " << disti << endl;
		if (disti < radiusi)
		{
			return;
		}
	}
	pixels[ix][iy][iz].is_available = 1;

}

void case_pixel_on_face(Point p, Cell_handle ch, Triangulation T)
{
	for (int i=0; i<3; i++)
	{
		float radiusi = ch->vertex(i)->info();
		float disti = (p.x() - (ch->vertex(i)->point().x())) * (p.x() - (ch->vertex(i)->point().x())) 
			+ (p.y() - (ch->vertex(i)->point().y())) * (p.y() - (ch->vertex(i)->point().y()))  
			+ (p.z() - (ch->vertex(i)->point().z())) * (p.z() - (ch->vertex(i)->point().z()));
		disti = sqrt(disti);
		cout << "distance from vertex " << i << " = " << disti << endl;
		if (disti < radiusi)
		{
			return;
		}
	}
	pixels[ix][iy][iz].is_available = 1;

}

void case_pixel_on_edge(Point p, Cell_handle ch, Triangulation T)
{
	for (int i=0; i<2; i++)
	{
		float radiusi = ch->vertex(i)->info();
		float disti = (p.x() - (ch->vertex(i)->point().x())) * (p.x() - (ch->vertex(i)->point().x())) 
			+ (p.y() - (ch->vertex(i)->point().y())) * (p.y() - (ch->vertex(i)->point().y()))  
			+ (p.z() - (ch->vertex(i)->point().z())) * (p.z() - (ch->vertex(i)->point().z()));
		disti = sqrt(disti);
		cout << "distance from vertex " << i << " = " << disti << endl;
		if (disti < radiusi)
		{
			return;
		}
	}
	pixels[ix][iy][iz].is_available = 1;

}

void case_pixel_on_vertex(Point p, Cell_handle ch, Triangulation T)
{

}

void explore_path()
{
	for (Delaunay::Finite_faces_iterator fit = triangulation.finite_faces_begin(); 
			fit != triangulation.finite_faces_end(); ++fit)
	{


	}

}

