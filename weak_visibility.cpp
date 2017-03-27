#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Triangular_expansion_visibility_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Line_2.h>
#include <iostream>
#include <vector>
#include <deque>
#include <map>
// Define the used kernel and arrangement  
typedef CGAL::Exact_predicates_exact_constructions_kernel       Kernel;
typedef Kernel::Point_2                                         Point_2;
typedef CGAL::Polygon_2<Kernel>                                 Polygon_2;
typedef Kernel::Segment_2                                       Segment_2;
typedef CGAL::Line_2<Kernel>                                    Line_2;
typedef CGAL::Arr_segment_traits_2<Kernel>                      Traits_2;
typedef CGAL::Arrangement_2<Traits_2>                           Arrangement_2;
typedef CGAL::Polygon_with_holes_2<Kernel>                      Polygon_with_holes_2;
typedef Arrangement_2::Halfedge_const_handle                    Halfedge_const_handle;
typedef Arrangement_2::Face_handle                              Face_handle;
typedef std::list<Polygon_with_holes_2>                         Pwh_list_2;
typedef Polygon_2::Vertex_iterator                              Vertex_iterator;
typedef Polygon_2::Edge_const_iterator                          Edge_iterator;

#include "print_utils.h"
// Define the used visibility class 
typedef CGAL::Triangular_expansion_visibility_2<Arrangement_2>  TEV;

int main() {
  // Defining the input geometry 
  std::vector<Segment_2> segments;
  int j,a,b,c,d,i,n,qn,qj,holes,points;
  Point_2 x,y;
  Polygon_2 outer;
  std::cout<<"give the polygon size: ";
  std::cin>>n;
  std::cout<<"Enter the each segment of outer polygon:"<<std::endl;
  Point_2 p[1000];
  for(i=0;i<n;++i){
    std::cin>>x;
    outer.push_back (x);
    std::cin>>y;
    Point_2 p1=x,p2=y;
    p[i] = p1;
    p[(i+1)%n] = p2;
    segments.push_back(Segment_2(p[i],p[(i+1)%n]));
  }
  int edg,edge;
  std::cin >> edge;
  edg = edge;
  Polygon_2 outers(outer),poly[2*edg];
  Polygon_with_holes_2 unionR;
  while(edg--){
  int quer,sd=2;
  while(sd--){
  std::map<Point_2,Point_2> Maps;
  std::multimap< Point_2, Point_2 > Maps2;
  std::map<Point_2,int> Maps3;
  std::map<Point_2,int> Maps4;
  std::vector<Polygon_2> visiPolygons;
  std::vector<Point_2> myqueue;
  std::deque<Polygon_2> set_polygons;
  set_polygons.push_back(outers);
  std::cout<<"Shortest path tree has to be caluculated for Point: ";
  std::cin>>quer;
  Point_2 query_point = p[quer];
  myqueue.push_back(query_point);

 while(!set_polygons.empty()) {
  outer = set_polygons.front();
  Arrangement_2 env;
  segments.clear();
  Edge_iterator vit;
  for(vit = outer.edges_begin(); vit != outer.edges_end(); vit++) {
      segments.push_back(*vit);
  }
  CGAL::insert_non_intersecting_curves(env,segments.begin(),segments.end());
  Halfedge_const_handle he = env.halfedges_begin();
  Vertex_iterator vin;
  for(std::vector<Point_2>::iterator vec = myqueue.begin(); vec != myqueue.end(); vec++) {
    int idj = 0;
    for(vin = outer.vertices_begin(); vin != outer.vertices_end(); vin++) {
      idj = 0;
      if(*vin == *vec) {query_point = *vin; idj = 1; break;}
    }
    if(idj) break;
  }
  
  while (he->target()->point() != query_point || he->face()->is_unbounded() )
      he++;
  Arrangement_2 output_arr;
  TEV tev(env);
  Face_handle fh = tev.compute_visibility(query_point, he, output_arr);
  Polygon_2 polygon;
  Arrangement_2::Ccb_halfedge_circulator curr = fh->outer_ccb();
  std::cout << "Visibility region from vertex: " << query_point << std::endl;
  polygon.push_back(Point_2(curr->source()->point()));
  std::cout << "[" << curr->source()->point() << " -> " << curr->target()->point() << "]" << std::endl;
  while (++curr != fh->outer_ccb()){
    polygon.push_back(Point_2(curr->source()->point()));
    std::cout << "[" << curr->source()->point() << " -> " << curr->target()->point() << "]"<< std::endl;
  }
  Vertex_iterator fi, afi, qfi;
  int flg = 0, mi = 0; std::size_t ct = -1;
  for(fi = polygon.vertices_begin(); fi != polygon.vertices_end(); fi++) {
    flg = 0; ++ct; mi = 0;
    for(afi = outers.vertices_begin(); afi != outers.vertices_end(); afi++) {
      if(*fi == *afi) {
        int dil = 0;
        for (std::map<Point_2,Point_2>::iterator its=Maps.begin(); its!=Maps.end(); ++its)
           {
              if( its->first == *fi) { 
                 dil = 1;
                 break;
              } 
           }
        if(!dil){
           std::cout << *fi << " => " << query_point << std::endl;
           Maps[*fi] = query_point;
           Maps2.insert(std::pair<Point_2,Point_2>(query_point,*fi));
           Maps3[query_point] = mi;
           if(!Maps4.count(query_point)) Maps4[query_point] = mi;
        }
        flg = 1; break;
      }
      ++mi;
    }
    if(!flg) {
      Point_2 p1,p2;
      if(ct != polygon.size()-1) p1 = polygon.vertex(ct+1);
      else p1 = *(polygon.vertices_begin());
      if(ct != 0) p2 = polygon.vertex(ct-1);
      else p2 = polygon.vertex(polygon.size()-1);
        if(CGAL::collinear(*fi, query_point, p1) && p1 != *fi && p1 != query_point){
          std::vector<Point_2>::iterator itr = myqueue.begin();
          int fl2 = 0;
          while(itr != myqueue.end()) {
            if(*itr == p1) fl2 = 1;
            ++itr;
          }
          if(!fl2) myqueue.push_back(p1);
        }
       if(CGAL::collinear(*fi, query_point, p2) && p2 != *fi && p2 != query_point){
          std::vector<Point_2>::iterator itr = myqueue.begin();
          int fl2 = 0;
          while(itr != myqueue.end()) {
            if(*itr == p2) fl2 = 1;
            ++itr;
          }
          if(!fl2) myqueue.push_back(p2);
        }
    }
  //for (std::map<Point_2,Point_2>::iterator itj = Maps.begin(); itj!= Maps.end(); ++itj)
    //std::cout << itj->first << " =><== " << itj->second << '\n'; 
  }
  Pwh_list_2 symmR;
  Pwh_list_2::const_iterator iter;
  CGAL::symmetric_difference (outer, polygon, std::back_inserter(symmR));
  for (iter = symmR.begin(); iter != symmR.end(); ++iter) {
    Polygon_2 polys = iter->outer_boundary();
    set_polygons.push_back(polys);
  }
  set_polygons.pop_front();
  }
  std::cout << "The required SPT is as follows:" << std::endl;
  for (std::map<Point_2,Point_2>::iterator it = Maps.begin(); it!= Maps.end(); ++it)
    std::cout << it->first << " => " << it->second << '\n';
  std::cout << "Done." <<std::endl;
  for (std::multimap<Point_2,int>::iterator it = Maps3.begin(); it!= Maps3.end(); ++it)
    std::cout << it->first << " <= " << p[it->second] << '\n';
  std::cout << "Done2." <<std::endl;
  std::multimap<Point_2,Point_2>::iterator it,itr;
  Point_2 r,q,m,q1,q2;
  if(sd) {
     r = p[(quer+1)%n];
     q = p[quer];
     m = p[(quer+2)%n]; 
     q2 = p[(quer+3)%n];
    }
  else {
     int quer2 = quer + n;
     r = p[(quer2-1)%n];
     q = p[quer];
     m = p[(quer2-2)%n]; 
     q2 = p[(quer2-3)%n];
  }
  Line_2 lin1(r,m);
  //std::cout << r << " " << q << " " << m << std::endl;
  int side = lin1.has_on_positive_side(q2);
  std::cout << side << std::endl ;
  it = Maps2.find(q);
  std::cout << it->first << "=>" << it->second << std::endl; 
  while(it != Maps2.end()) {
    if(it->first == it->second){
      Maps2.erase(it);
      it = Maps2.find(q);
    } else{
    itr = Maps2.find(it->second);
    if(itr == Maps2.end()){
      Maps2.erase(it);
      it = Maps2.find(q);
    }
    else{
      //std::cout << it->first << " wins  " << it->second << std::endl;
      Line_2 lin2(it->first, it->second);
      q1 = it->second;
      if(side){
        if(lin2.has_on_positive_side(itr->second)) it = Maps2.find(q1);
        else break;
      } else {
          if(!lin2.has_on_positive_side(itr->second)) it = Maps2.find(q1);
          else break;
      }
      if(it == Maps2.end()) it = Maps2.find(q);
    }
  }}
  std::cout << "extend this " << it->first << " =>> " << it->second << std::endl;
  Line_2 lin2(it->first,it->second);
  Point_2 tmpt = itr->second;
  Point_2 creo1 = itr->first, creo2 = itr->second;
  Maps2.erase(itr);
  itr = Maps2.find(it->second);
  while(itr != Maps2.end()) {
    Line_2 lins(creo1,creo2);
    if(lins.has_on_positive_side(itr->second) == side){
        tmpt = itr->second;
        creo1 = itr->first, creo2 = itr->second;
    }
      Maps2.erase(itr);
      itr = Maps2.find(it->second);     
  }
  std::cout << "new line pt: "<< tmpt << std::endl;
  Point_2 tmpt2;
  for(Vertex_iterator vir = outers.vertices_begin(); vir != outers.vertices_end(); ++vir) {
    if(*vir == tmpt){
     if(side) tmpt2 = *(++vir);
     else tmpt2 = *(--vir);
     break;
    }
  }
  std::cout << "other pt: " << tmpt2 <<std::endl;
  Line_2 new_lin(tmpt,tmpt2);
  CGAL::cpp11::result_of<Kernel::Intersect_2(Line_2, Line_2)>::type itrsc;
  itrsc = CGAL::intersection(lin2,new_lin);
  const Point_2* intrsctn = boost::get<Point_2>(&*itrsc);
  std::cout << "req new point of intersection: " << *intrsctn << std::endl;
  int qwer1 = 0, sd2 = edg*2+sd, quer1 = quer;
  std::cout << tmpt << " " << it->second << std::endl;
  Point_2 pnt, pnt2;
  for(pnt = p[quer1]; pnt != tmpt; ++quer1){
    pnt = p[quer1%n];
    std::cout << "add " << pnt << std::endl;
    poly[sd2].push_back(pnt);
  }
  poly[sd2].push_back(*intrsctn);
  std::cout << "viw " << pnt << std::endl;
  if(pnt == tmpt) pnt2 = it->second;
  else pnt2 = tmpt;
  quer1 = quer;
  for(pnt = p[(quer1-1)%n]; pnt != pnt2; --quer1){
    pnt = p[quer1%n];
  }
  for(pnt = p[(quer1+1)%n]; pnt != p[(quer-1)%n]; ++quer1){
    pnt = p[quer1%n];
    poly[sd2].push_back(pnt);
  }
  //for(Vertex_iterator vits = outers.vertices_begin(); vits != outers.vertices_end(); ++vits){
  //  if(*vits == tmpt2 || *vits == it->second){
  //    poly[sd2].push_back(*vits);
  //    if(!qwer1) {
  //      poly[sd2].push_back(*intrsctn);
  //      qwer1 = 1;
  //    } else qwer1 = 0;
  //  }
  //  else{
  //    if(!qwer1) poly[sd2].push_back(*vits);
  //  }
  //}
  for(Vertex_iterator vits = poly[sd2].vertices_begin(); vits != poly[sd2].vertices_end(); ++vits){
    std::cout << *vits << " plo " << edg*2+sd << std::endl;
  }
  }
  Polygon_with_holes_2 wk_vis_poly;
  Pwh_list_2                  intR;
  Pwh_list_2::const_iterator  it5;
  CGAL::intersection (poly[edg*2], poly[edg*2+1], std::back_inserter(intR));
  std::cout << "The intersection:" << std::endl;
  for (it5 = intR.begin(); it5 != intR.end(); ++it5) {
    std::cout << "--> ";
    print_polygon_with_holes (*it5);
    std::cout << "union poly: " << std::endl;
    if(edg == edge-1){ if(CGAL::join(*it5,*it5,unionR)){
      print_polygon_with_holes(unionR);
    }}else{if(CGAL::join(unionR,*it5,unionR)){
      print_polygon_with_holes(unionR);
    }}
  }
  }
  Pwh_list_2 symmR;
  Pwh_list_2::const_iterator it9;
  CGAL::symmetric_difference (outers, unionR, std::back_inserter(symmR));
  if(symmR.size() == 0){
    std::cout<<"Edge Guard set holds"<<std::endl;
  }
  else{
    std::cout << "\n\nThe symmetric difference:" << std::endl;
    for (it9 = symmR.begin(); it9 != symmR.end(); ++it9) {
      std::cout << "--> ";
      print_polygon_with_holes (*it9);
    }
    std::cout<<"Edge Guard set is invalid"<<std::endl;
  }
  return 0;
}
