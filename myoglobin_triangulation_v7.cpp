#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Origin.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/cartesian_homogeneous_conversion.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include <vector>
#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Triangulation_geomview_ostream_3.h>
#include <CGAL/intersections.h>
#include <CGAL/squared_distance_3.h> 
#include <limits>
#include <CGAL/Line_3.h>
#include <CGAL/utils.h>
#include <unordered_set>
#include <set>
#include <map>
#include <algorithm>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <iterator>
#include <CGAL/Vector_3.h>
#include <CGAL/enum.h>

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
typedef K::Line_3 Line_3;
typedef K::Plane_3 Plane_3;
typedef K::Vector_3 Vector_3;
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
using std::ofstream;
using CGAL::intersection;
using CGAL::Object;
using CGAL::ORIGIN;
using std::insert_iterator;


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

map<Facet,  Point> cavity_location;
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
void get_adjacent_cells(set<Facet>& facets, set<string>& cells, map<string, Cell_handle>&,  Triangulation& T);
void get_adjacent_facets(set<string>& diffused_cells, map<string, Cell_handle>& name2ch, set<Facet>& diffused_facets, Triangulation& T);
string cell2string(Cell& c);
void draw_faces(set<Facet>& facets, ofstream& ofcav, ofstream& ofedge);
void draw_face(const Facet& facet, ofstream& ofcav, ofstream& ofedge);
void draw_atoms(set<string>& diffused_cells, map<string, Cell_handle>& name2ch, ofstream& ofs,Triangulation& T);
void draw_atoms(set<Facet>& facets, ofstream& ofs);
void generate_grid_input();
void generate_packed_sphere();
Point read_modified_pdb(char* fn);

double minx, miny, minz, maxx, maxy, maxz;
int Nx, Ny, Nz;
double probe_rad;

int main()
{
	//read_input("expt_input_point.txt");
	//generate_grid_input();
//	generate_packed_sphere();
	
	Point target_loc = read_modified_pdb("1A6G.pqr");
	
	//probe_rad = 0.4;   	
	double epsilon = 0.105;
	probe_rad = (2.0/sqrt(3.0)) - 1.0 - epsilon;
	probe_rad = 1.4;
	
	//cout << minx << " " << maxx << " " << miny << " " << maxy << " " << minz << " " << maxz << endl;
	Triangulation T(L.begin(), L.end());
	cout<< "outputting T: " << endl;
	cout << T <<endl;
	Triangulation::size_type nverts = T.number_of_vertices();
	Triangulation::size_type ncells = T.number_of_finite_cells();
	Triangulation::size_type nfaces = T.number_of_facets();
	//iterate_over_faces(T);
 	cout << "Total number of tetrahedrons (or, cells)" << ncells << " "<<T.number_of_cells()<<endl;
	cout << "Outputting the  triangulation with coordinates " << endl;
	Triangulation::All_cells_iterator bit = T.all_cells_begin();
	Triangulation::All_cells_iterator eit = T.all_cells_end();
	Triangulation::All_cells_iterator it;
	int count = 1;
	/*for(it = bit ; it != eit; ++it)
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
	}*/

	Locate_type lt;
	int li, lj;
	Vertex_handle target_vertex;
	Cell_handle target_ch  = T.locate(target_loc, lt, li, lj); //TODO: make it a Fe location
	set<string> target_cell_names;
	vector<Cell_handle> all_target_cells;
	switch(lt)
	{
		case Triangulation::VERTEX:
			cout<<"target is a vertex "<<endl;
			target_vertex = target_ch->vertex(li);
			T.incident_cells(target_vertex, std::back_inserter(all_target_cells));

			for(int i = 0;  i < all_target_cells.size(); ++i)
			{
				target_cell_names.insert(cell2string(*(all_target_cells[i])));

			}
			break;
		
		case Triangulation::CELL:
		{
			cout<<"target is a cell"<<endl;
			target_cell_names.insert(cell2string(*(target_ch)));

		}

	}
	/*if(assign(target_vertex, o))
	{
		cout<<"target is a vertex "<<endl;

		vector<Cell_handle> all_target_cells;
		T.incident_cells(target_vertex, std::back_inserter(all_target_cells));

		for(int i = 0;  i < all_target_cells.size(); ++i)
		{
			target_cell_names.insert(cell2string(*(all_target_cells[i])));

		}

	}
	if(assign(target_ch, o))
	{
		cout<<"target is a cell"<<endl;
		target_cell_names.insert(cell2string(*(target_ch)));

	}*/
	cout<<"Number of target cells "<<target_cell_names.size()<<endl;
	
	vector<Facet> all_outside_facets;
	set<Facet>  diffused_facets;
  	get_outside_faces(T, all_outside_facets, diffused_facets);
	/*cout << "Listing all outside faces ... " << endl;
	for (int i=0; i < all_outside_facets.size(); i++)
	{
		Facet f = all_outside_facets[i];
		Cell c = *(f.first);
		cout << "pivot ["<<*(c.vertex(f.second))<<"] : ["<<*(c.vertex(0))<<"]  ["<<*(c.vertex(1))<<"] ["<<*(c.vertex(2))<<"]  ["<<*(c.vertex(3)) << "]" << endl;
	} */	
	ofstream edge_ofs("edge_3D_iter_0.dat");
	ofstream cavity_ofs("cavity_3D_iter_0.dat");
	ofstream oa_ofs1("atom_3D_iter_0.dat");
	draw_faces(diffused_facets, cavity_ofs, edge_ofs);
	draw_atoms(diffused_facets, oa_ofs1);

	set<string> diffused_cells;
	map<string,  Cell_handle> name2ch;
	count = 1;
	int prev_cell_count  = 0;
	while(1)
	{
		get_adjacent_cells(diffused_facets, diffused_cells, name2ch, T);
		if(prev_cell_count == diffused_cells.size()) break;
		prev_cell_count = diffused_cells.size();
		set<string>::iterator sit;
		cout << "Num Diffused cells inserted in iteration " << count << ": "<<diffused_cells.size()<<endl;
		

		get_adjacent_facets(diffused_cells, name2ch, diffused_facets, T);
		//cout << "Diffused faces inserted in iteration " << count << endl;
		cout << " Num Diffused faces inserted in iteration " << count <<": "<<diffused_facets.size()<< endl;
		set<Facet>::iterator fit;
	  	/*cout << "Current list of facets" << endl;
		for (fit = diffused_facets.begin(); fit != diffused_facets.end(); fit++)
		{
			const Facet& f = *fit;
			Cell& c = *f.first;
			cout << "pivot ["<<*(c.vertex(f.second))<<"] : ["<<*(c.vertex(0))<<"]  ["<<*(c.vertex(1))<<"] ["<<*(c.vertex(2))<<"]  ["<<*(c.vertex(3)) << "]" << endl;
		}*/
		ostringstream oefn_cav, oefn_edge; 
		oefn_edge<<"edge_3D_iter_"<<count<<".dat";
		ofstream oe_ofs_edge(oefn_edge.str());
		oefn_cav<<"cavity_3D_iter_"<<count<<".dat";
		ofstream oe_ofs_cav(oefn_cav.str());
		ostringstream oatm_str;
		oatm_str<<"atom_3D_iter_"<<count<<".dat";
		ofstream oa_ofs(oatm_str.str());
		draw_faces(diffused_facets, oe_ofs_cav, oe_ofs_edge);
		//draw_atoms(diffused_cells, name2ch, oa_ofs, T);
		draw_atoms(diffused_facets, oa_ofs);
		cout << "All Diffused cells in iteration " << count << endl;
		for (sit = diffused_cells.begin(); sit != diffused_cells.end(); sit++)
		{
			//cout << *sit << endl;
			/*if(*sit == target_cell_name){
				
				set<string>::iterator bit = diffused_cells.find(target_cell_name);
				if(bit == diffused_cells.end()){cout<<"target cell not found "<<endl;}
				if(bit != diffused_cells.end()){cout<<"target cell found "<<endl;}
				std::cout<<"target cell "<< *sit<<std::endl;
				break;
			}*/
		}
		set<string> found_targets; 	
		std::set_intersection(diffused_cells.begin(), diffused_cells.end(), target_cell_names.begin(), target_cell_names.end(), std::inserter(found_targets, found_targets.begin()));
		if(found_targets.size() > 0)
		{
			cout<<"found num targets "<<found_targets.size()<<endl;
			break;
		}
		count++;

	}
	cout<<"num diffused cells "<<diffused_cells.size()<<endl;
	////// possibly useful functions ////////
	//is_edge(cell_handle, int i, int j)/////
	//is_facet(u,v,w,&cell_handle, &i,&j,&k) ////
	//facet_circulator : all_facets_begin()///
	// build the search box

	return 0;

}

bool explore_face_edges(Facet& facet, Triangulation& P)
{
	//return true;
	Cell_handle ch = facet.first;
	int  opp_vertex = facet.second;
	float tx, ty, tz;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if ((i != opp_vertex) && (j != opp_vertex) && (i != j))
			{
				int k = 6 - opp_vertex - i - j;
				float rad_i = ch->vertex(i)->info();
				float rad_j = ch->vertex(j)->info();
				float rad_k = ch->vertex(k)->info();
				Point pi = ch->vertex(i)->point();
				Point pj = ch->vertex(j)->point();
				Point pk = ch->vertex(k)->point();
				//cout<<"Check "<<rad_i<<" "<<rad_j<<" "<<rad_k<<endl;
				tx = (ch->vertex(i)->point().x() - ch->vertex(j)->point().x())*(ch->vertex(i)->point().x() - ch->vertex(j)->point().x());
				ty = (ch->vertex(i)->point().y() - ch->vertex(j)->point().y())*(ch->vertex(i)->point().y() - ch->vertex(j)->point().y());
				tz = (ch->vertex(i)->point().z() - ch->vertex(j)->point().z())*(ch->vertex(i)->point().z() - ch->vertex(j)->point().z());

				float len_ij = sqrt(tx + ty + tz); 
				if (len_ij - (rad_i + rad_j) >=  2.0 * probe_rad)
				{
					// check whether the opposite vertex atom is blocking
					// * Is the oppsite vertex even inside the base?:
					//CGAL::angle(const CGAL::Point_3< Kernel > &pi, const CGAL::Point_3< Kernel > &pj, const CGAL::Point_3< Kernel > &pk);
					CGAL::Angle ang = CGAL::angle(pi, pj, pk);

					// IF it is inside, then DOES it block?:
					if ((ang == CGAL::ACUTE) || (ang == CGAL::ACUTE))
					{
						Line_3 l(pi, pj);
						double dist = sqrt(CGAL::squared_distance(pk, l));
						if (rad_k >= dist)
						{
							// find the intersection and edge availability
						}
						else if ((rad_k < dist) && (dist - rad_k) >= (2.0*1.4) )
						{
							return true;
						}
						/*else if ((rad_k < dist) && (dist - rad_k) < (2.0*1.4) )
						  {
						// check the possibility of squeezing in from either sides of this edge
						// 1. Find the intersection of the sphere and the segment.
						// 2. Check the residual area between the i-sphere and the k-sphere
						// 3. Check the residual area between the j-sphere and the k-sphere

						}
						else
						{

						}*/
					}
					else
					{
						// (iii) if it is outside the base, then?? 

					}

					return true;
				}
			}
		}
	}
	int i,j,k;
	if(opp_vertex == 0) { i = 1; j = 2; k =3;}
	if(opp_vertex == 1) { i = 0; j = 2; k =3;}
	if(opp_vertex == 2) { i = 0; j = 1; k =3;}
	if(opp_vertex == 3) { i = 0; j = 1; k =2;}

					// side (i, j)
	Point p0(0,0,0);
	Point pi = ch->vertex(i)->point();
	Point pj = ch->vertex(j)->point();
	Point pk = ch->vertex(k)->point();
	float rad_i = ch->vertex(i)->info();
	float rad_j = ch->vertex(j)->info();
	float rad_k = ch->vertex(k)->info();
	float len_ij = sqrt(CGAL::squared_distance(pi, pj));
	float len_ik = sqrt(CGAL::squared_distance(pi, pk));
	float len_jk = sqrt(CGAL::squared_distance(pj, pk));
	Vector_3 vi(p0, pi);
	Vector_3 vj(p0, pj);
	Vector_3 vk(p0, pk);
	Point seg_ij_l = ORIGIN + ((rad_i/len_ij) * (vj-vi) + vi) ;
	Point seg_ik_l = ORIGIN + ((rad_i/len_ik) * (vk-vi) + vi) ;
	Point seg_jk_l = ORIGIN + ((rad_j/len_jk) * (vk-vj) + vj) ;

	Point seg_ij_r = ORIGIN + (((len_ij-rad_j)/len_ij) * (vj-vi) + vi) ;
	Point seg_ik_r = ORIGIN + (((len_ik- rad_k)/len_ik) * (vk-vi) + vi) ;
	Point seg_jk_r = ORIGIN + (((len_jk-rad_k)/len_jk) * (vk-vj) + vj) ;

	Point mid_ij = ORIGIN + ((Vector_3(p0, seg_ij_l) + Vector_3(p0, seg_ij_r))/2.0);
	Point mid_ik = ORIGIN + ((Vector_3(p0, seg_ik_l) + Vector_3(p0, seg_ik_r))/2.0);
	Point mid_jk = ORIGIN + ((Vector_3(p0, seg_jk_l) + Vector_3(p0, seg_jk_r))/2.0);

	Line_3 l_ij(pi, pj);
	Plane_3 p_ij = l_ij.perpendicular_plane(mid_ij);

	Line_3 l_ik(pi, pk);
	Plane_3 p_ik = l_ik.perpendicular_plane(mid_ik);

	Line_3 l_jk(pj, pk);
	Plane_3 p_jk = l_jk.perpendicular_plane(mid_jk);

	Line_3 plane_intersect_ij_jk ; 
	Object o = intersection(p_ij, p_jk);
	if(CGAL::assign(plane_intersect_ij_jk, o))
	{
		//Object p = intersection(plane_intersect_ij_jk, p_ik);
		Plane_3 face_plane(pi,pj,pk);
		Object p = intersection(face_plane, plane_intersect_ij_jk);
		Point cav_center;
		if(CGAL::assign(cav_center, p))
		{

			double sq_distance_i = sqrt(CGAL::squared_distance(cav_center, pi)) - rad_i;
			double sq_distance_j = sqrt(CGAL::squared_distance(cav_center, pj)) - rad_j;
			double sq_distance_k = sqrt(CGAL::squared_distance(cav_center, pk)) - rad_k;
			double max_cav_rad = CGAL::min(sq_distance_i, CGAL::min(sq_distance_j, sq_distance_k));
			cavity_location[facet] = cav_center;
			if(max_cav_rad > probe_rad) return true;
			return false; 

		}
		else
		{
			cout<<plane_intersect_ij_jk<< " "<<plane_intersect_ij_jk.direction()<<endl;
			cout<<p_ik<<" "<<endl;
			Line_3 ltmp;
			bool is_line = CGAL::assign(ltmp, p);
			Plane_3 lplane;	
			bool is_plane = CGAL::assign(lplane, p);
			cout<<"is line "<<is_line<<endl;
			cout<<"is plane"<<is_plane<<endl;
			assert(0);
		}
	}
	else
	{
		assert(0);
	}
	return false;
}

void get_adjacent_facets(set<string>& diffused_cells, map<string, Cell_handle>& name2ch, set<Facet>& diffused_facets, Triangulation& T)
{
	set<string>::iterator sit;
	for(sit = diffused_cells.begin(); sit != diffused_cells.end(); ++sit)
	{
		Cell_handle ch = name2ch[*sit];
		for(int i = 0; i < 4; ++i)
		{
				Facet f(make_pair(ch, i));
				if ((!T.is_infinite(f)) && explore_face_edges(f, T))
					diffused_facets.insert(f);
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
			if (explore_face_edges(outside_face, T))
			{
				all_outside_facets.push_back(outside_face);
				set_of_outside_faces.insert(outside_face);
			}
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

void draw_atoms(set<string>& diffused_cells, map<string, Cell_handle>& name2ch, ofstream& ofs,Triangulation& T)
{
	set<string>::iterator cit = diffused_cells.begin();
	for(; cit != diffused_cells.end(); ++cit)
	{
		Cell_handle ch = name2ch[*cit];
		for(int i = 0; i < 4; ++i)
		{
			Point p = ch->vertex(i)->point();
			ofs<<p.x()<<" "<<p.y()<<" "<<p.z()<<" "<<ch->vertex(i)->info()<<endl;

		}

	}
}
void draw_atoms(set<Facet>& facets, ofstream& ofs)
{
	set<Facet>::iterator fit = facets.begin();
	set<Facet>::iterator eit = facets.end();
	for (; fit!=eit; ++fit)
	{
		Facet f = *fit;
		Cell_handle ch = f.first;
		for(int i = 0; i < 4; ++i)
		{
			if(i != f.second)
			{
				Point p = ch->vertex(i)->point();
				ofs<<p.x()<<" "<<p.y()<<" "<<p.z()<<" "<<ch->vertex(i)->info()<<endl;
				
			}

		}
	}
}
void draw_faces(set<Facet>& facets, ofstream& ofcav, ofstream& ofedge)
{
	set<Facet>::iterator fit = facets.begin();
	set<Facet>::iterator eit = facets.end();
	for (; fit!=eit; ++fit)
	{
		Facet f = *fit;
		Cell_handle ch = f.first;
		draw_face(f, ofcav, ofedge);
	}

}

void draw_face(const Facet& facet, ofstream& ofcav, ofstream& ofedge)
{
	Cell_handle ch = facet.first;
	int  opp_vertex = facet.second;
	float tx, ty, tz;
	float xitil, yitil, zitil, xjtil, yjtil, zjtil;
	float len_ij, pi, pj, rad_i=0.2, rad_j=0.2;
	int cavity_type = 0;
	if(cavity_location.find(facet) != cavity_location.end()) cavity_type = 1;
	if(cavity_type == 1)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				if ((i != opp_vertex)&&(j != opp_vertex) && (i != j))
				{
					// write edges to edge file
					ofedge << ch->vertex(i)->point().x() << "  " << ch->vertex(i)->point().y() << "  " << ch->vertex(i)->point().z() << "  " <<
						ch->vertex(j)->point().x() << "  " << ch->vertex(j)->point().y() << "  " << ch->vertex(j)->point().z() << endl;
				}
			}
		}
				
		Point cc = cavity_location[facet];
		ofcav << cc.x() << "  " << cc.y() << "  " << cc.z() << "  " << 0 << "  " << 0 << "  " << 0 << endl;
		return ;
	}
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if ((i != opp_vertex)&&(j != opp_vertex) && (i != j))
			{
				// write edges to edge file
				
				float rad_i = ch->vertex(i)->info();
				float rad_j = ch->vertex(j)->info();
				
				Point p0(0,0,0);
				Point pi = ch->vertex(i)->point();
				Point pj = ch->vertex(j)->point();
				float len_ij = sqrt(CGAL::squared_distance(pi, pj));
				Vector_3 vi(p0, pi);
				Vector_3 vj(p0, pj);
				Point l = ORIGIN + ((rad_i/len_ij) * (vj-vi) + vi) ;
				Point r = ORIGIN + (((len_ij-rad_j)/len_ij) * (vj-vi) + vi) ;

				Point mid = ORIGIN + ((Vector_3(p0, l) + Vector_3(p0, r))/2.0);
				if (len_ij - (rad_i + rad_j) >=  2.0 * probe_rad)
				{
					ofedge << ch->vertex(i)->point().x() << "  " << ch->vertex(i)->point().y() << "  " << ch->vertex(i)->point().z() << "  " <<
						ch->vertex(j)->point().x() << "  " << ch->vertex(j)->point().y() << "  " << ch->vertex(j)->point().z() << endl;
					ofcav << mid.x() << "  " << mid.y()<< "  " << mid.z() << "  " << rad_i << "  " << rad_j << "  " << len_ij << endl;
				}
				/*tx = (ch->vertex(i)->point().x() - ch->vertex(j)->point().x())*(ch->vertex(i)->point().x() - ch->vertex(j)->point().x());
				ty = (ch->vertex(i)->point().y() - ch->vertex(j)->point().y())*(ch->vertex(i)->point().y() - ch->vertex(j)->point().y());
				tz = (ch->vertex(i)->point().z() - ch->vertex(j)->point().z())*(ch->vertex(i)->point().z() - ch->vertex(j)->point().z());
				len_ij = sqrt(tx + ty + tz); 
				
				pi = rad_i/len_ij;	
				xitil = (ch->vertex(i)->point().x()) + pi * (ch->vertex(j)->point().x() - ch->vertex(i)->point().x());
				yitil = (ch->vertex(i)->point().y()) + pi * (ch->vertex(j)->point().y() - ch->vertex(i)->point().y());
				zitil = (ch->vertex(i)->point().z()) + pi * (ch->vertex(j)->point().z() - ch->vertex(i)->point().z());

				pj = 1.0 - (rad_j/len_ij);	
				xjtil = (ch->vertex(i)->point().x()) + pj * (ch->vertex(j)->point().x() - ch->vertex(i)->point().x());
				yjtil = (ch->vertex(i)->point().y()) + pj * (ch->vertex(j)->point().y() - ch->vertex(i)->point().y());
				zjtil = (ch->vertex(i)->point().z()) + pj * (ch->vertex(j)->point().z() - ch->vertex(i)->point().z());

				tx = (xitil + xjtil)/2.0;
				ty = (yitil + yjtil)/2.0;
				tz = (zitil + zjtil)/2.0;	*/
				/*tx =  (ch->vertex(i)->point().x() + ch->vertex(j)->point().x())/2.0;
				ty =  (ch->vertex(i)->point().y() + ch->vertex(j)->point().y())/2.0;
				tz =  (ch->vertex(i)->point().z() +  ch->vertex(j)->point().z())/2.0; */
		
			}
		}

	}
}

void explore_path()
{
	//for (Delaunay::Finite_faces_iterator fit = triangulation.finite_faces_begin(); 
	//		fit != triangulation.finite_faces_end(); ++fit)
	//{
		

	//}

}

void get_adjacent_cells(set<Facet>& facets, set<string>& cells, map<string,Cell_handle>& name2ch, Triangulation& T)
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


string cell2string(Cell& c)
{
	ostringstream cstr;
	for(int i = 0; i < 4; ++i)
	{	
		cstr<<*(c.vertex(i))<<":";
	}
	
	return cstr.str();

}

Point read_modified_pdb(char* fn)
{
	ifstream ifs(fn);
	string c1_atm;
	int c2_dummy;
	string c3_atmname;
	string c4_atmtype;
	int c5_dummy;
	
	ofstream ofs("atoms_data.dat");
	Point fe_loc;
	double x, y, z, q, r, rmax = 0.0;
	minx=numeric_limits<double>::max() , miny=numeric_limits<double>::max(), minz=numeric_limits<double>::max();
	maxx=numeric_limits<double>::min(), maxy=numeric_limits<double>::min(),maxz=numeric_limits<double>::min();
	char aline[5000];
	ifs.getline(aline, 5000);
	ifs>>c1_atm>>c2_dummy>>c3_atmname>>c4_atmtype>>c5_dummy>>x>>y>>z>>q>>r;
	while(ifs)
	{
		ofs << x << " " << y << " " << z << " " << r << endl;
		if(c1_atm == "ATOM")
		{
			L.push_front(make_pair(Point(x,y,z), r));
			if (c3_atmname == "FE")
			{
				// set target
				cout<<c1_atm<<" "<<c3_atmname<<" "<<x<<" "<<y<<" "<<z<<" "<<r<<endl;
				fe_loc = Point(x,y,z);
			}
		}
		
		//cout<<c1_atm<<" "<<c3_atmname<<" "<<x<<" "<<y<<" "<<z<<" "<<r<<endl;
		
		ifs>>c1_atm>>c2_dummy>>c3_atmname>>c4_atmtype>>c5_dummy>>x>>y>>z>>q>>r;


	}
	return fe_loc;

}
void read_input(char* fn)
{
	ifstream ifs(fn);
	double x, y, z, r, rmax = 0.0;
	minx=numeric_limits<double>::max() , miny=numeric_limits<double>::max(), minz=numeric_limits<double>::max();
	maxx=numeric_limits<double>::min(), maxy=numeric_limits<double>::min(),maxz=numeric_limits<double>::min();
	ifs >> x >> y >> z >> r;
	while(ifs)
	{
		L.push_front(make_pair(Point(x,y,z), r));
		cout << x << " " << y << " " << z << endl;
		ifs >> x >> y >> z >> r;
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
		if (r > rmax)
		{
			r = rmax;
		}

	}
	minx = minx - rmax;
	maxx = maxx + rmax;
	miny = miny - rmax;
	maxy = maxy + rmax;
	minz = minz - rmax;
	maxz = maxz + rmax;

}

void generate_grid_input()
{
	int Nz = 5, Nx = 5, Ny = 5;
	float xmin = -2.5, ymin = -2.5, zmin = -2.5;
	float x, y, z, dx = 1.1, r = 0.2;
	ofstream ofs("atoms_data.dat");

	for (int i = 0; i < Nx; i++)
	{
	    x = xmin + i * dx; 
	    for (int j = 0; j < Ny; j++)
	    {
		y = ymin + j * dx; 
		for (int k = 0; k < Nz; k++)
		{
			z = zmin + k * dx;
			float toss = drand48();
			if (toss < 0.3)
			{
				L.push_front(make_pair(Point(x,y,z), r));
				ofs << x << " " << y << " " << z << " " << r << endl;
			}
			
		}
	    }
	}


}

void generate_packed_sphere()
{

	ofstream ofs("atoms_data.dat");
	int side = 16;
	double r = 1.0;
	for(int i = 0; i < side; ++i)
	{
		for(int j = 0; j < side; ++j)
		{
			for(int k = 0; k < side; ++k)
			{
				double x = 2 * i + j%2 + k%2;
				double y = sqrt(3.0) * (j + ((1.0/3.0) * (k%2)));
				double z = (2 * sqrt(6))/3.0 * k;
				float toss = drand48();
				if (toss < 0.3)
				{
					L.push_front(make_pair(Point(x,y,z), r));
					ofs << x << " " << y << " " << z << " " << r << endl;
				}
				//cout<<x<<" "<<y<<" "<<z<<endl;

			}
		}
	}


}
/*void get_bounding_faces(unordered_set<string>& cells, set<Facet>& facets, Triangulation& T)
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

}*/

