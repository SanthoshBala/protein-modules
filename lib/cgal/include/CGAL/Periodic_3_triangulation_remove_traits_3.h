// Copyright (c) 2009   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/releases/CGAL-4.1-branch/Periodic_3_triangulation_3/include/CGAL/Periodic_3_triangulation_remove_traits_3.h $
// $Id: Periodic_3_triangulation_remove_traits_3.h 67117 2012-01-13 18:14:48Z lrineau $
// 
//
// Author(s)     : Manuel Caroli <Manuel.Caroli@sophia.inria.fr>

#ifndef CGAL_PERIODIC_3_TRIANGULATION_REMOVE_TRAITS_3_H
#define CGAL_PERIODIC_3_TRIANGULATION_REMOVE_TRAITS_3_H

#include <CGAL/basic.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Periodic_3_offset_3.h>

namespace CGAL { 

template < class Traits_, class Predicate_ >
class Point_offset_adaptor {
  typedef Traits_ Traits;
  typedef Predicate_ Predicate;

  typedef typename Traits::Point_3       Point;

public:
  typedef typename Predicate::result_type result_type;

 Point_offset_adaptor(const Predicate & pred) : _pred(pred) {}

  result_type operator()(const Point& p0, const Point& p1) const {
    return _pred(p0.first, p1.first,
	p0.second, p1.second);
  }
  result_type operator()(const Point& p0, const Point& p1,
      const Point& p2) const {
    return _pred(p0.first, p1.first, p2.first,
	p0.second, p1.second, p2.second);
  }
  result_type operator()(const Point& p0, const Point& p1,
      const Point& p2, const Point& p3) const {
    return _pred(p0.first, p1.first, p2.first, p3.first,
	p0.second, p1.second, p2.second, p3.second);
  }
  result_type operator()(const Point& p0, const Point& p1,
      const Point& p2, const Point& p3, const Point& p4) const {
    return _pred(p0.first, p1.first, p2.first, p3.first, p4.first,
	p0.second, p1.second, p2.second, p3.second, p4.second);
  }

private:
  Predicate _pred;
};

template < class P3DTTraits, class Off = typename CGAL::Periodic_3_offset_3 >
class Periodic_3_triangulation_remove_traits_3 : public P3DTTraits::K
{
public:
  typedef P3DTTraits                                            PT;
  typedef typename P3DTTraits::K                                Base;
  typedef Off                                                   Offset;
  typedef Periodic_3_triangulation_remove_traits_3< PT,Offset > Self;  

  typedef typename PT::RT                RT;
  typedef typename PT::FT                FT;
  typedef typename PT::Point_3           Bare_point;
  typedef std::pair<Bare_point,Offset>       Point_3;
  typedef typename PT::Iso_cuboid_3      Iso_cuboid_3;

  Periodic_3_triangulation_remove_traits_3(const Iso_cuboid_3& domain)
    : _pt() {
    _pt.set_domain(domain);
  }

  // Triangulation traits
  typedef Point_offset_adaptor<Self, typename PT::Compare_xyz_3>
      Compare_xyz_3;
  typedef Point_offset_adaptor<Self, typename PT::Coplanar_orientation_3>
      Coplanar_orientation_3;
  typedef Point_offset_adaptor<Self, typename PT::Orientation_3>
      Orientation_3;
  
  // Delaunay Triangulation traits
  typedef Point_offset_adaptor<Self,
	  typename PT::Coplanar_side_of_bounded_circle_3>
      Coplanar_side_of_bounded_circle_3;
  typedef Point_offset_adaptor<Self, typename PT::Side_of_oriented_sphere_3>
      Side_of_oriented_sphere_3;
  typedef Point_offset_adaptor<Self, typename PT::Compare_distance_3>
      Compare_distance_3;

  // Operations
  Compare_xyz_3
  compare_xyz_3_object() const {
    return Compare_xyz_3(_pt.compare_xyz_3_object());
  }
  Coplanar_orientation_3
  coplanar_orientation_3_object() const {
    return Coplanar_orientation_3(_pt.coplanar_orientation_3_object());
  }
  Orientation_3
  orientation_3_object() const {
    return Orientation_3(_pt.orientation_3_object());
  }
  Coplanar_side_of_bounded_circle_3
  coplanar_side_of_bounded_circle_3_object() const {
    return Coplanar_side_of_bounded_circle_3(_pt.coplanar_side_of_bounded_circle_3_object());
  }
  Side_of_oriented_sphere_3
  side_of_oriented_sphere_3_object() const {
    return Side_of_oriented_sphere_3(_pt.side_of_oriented_sphere_3_object());
  }
  Compare_distance_3
  compare_distance_3_object() const {
    return Compare_distance_3(_pt.compare_distance_3_object());
  }
public:
  PT _pt;
};

} //namespace CGAL

#endif // CGAL_PERIODIC_3_TRIANGULATION_REMOVE_TRAITS_3_H
