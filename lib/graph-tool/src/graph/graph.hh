// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006-2013 Tiago de Paula Peixoto <tiago@skewed.de>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#ifndef GRAPH_HH
#define GRAPH_HH
#include "config.h"

#include <Python.h>
#include <boost/python/object.hpp>
#include <boost/python/dict.hpp>

#include <deque>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include "fast_vector_property_map.hh"
#include <boost/variant.hpp>
#include <boost/mpl/vector.hpp>
#include "graph_properties.hh"
#include "graph_exceptions.hh"

namespace graph_tool
{
using namespace std;
using namespace boost;

// GraphInterface
// this class is the main interface to the internally kept graph. This is how
// the external world will manipulate the graph. All the algorithms should be
// registered here. This class will be exported to python in graph_bind.hh

namespace detail
{
// Generic graph_action functor. See graph_filtering.hh for details.
template <class Action, class GraphViews, class Wrap = mpl::false_,
          class TR1=boost::mpl::vector<>, class TR2=boost::mpl::vector<>,
          class TR3=boost::mpl::vector<>, class TR4=boost::mpl::vector<> >
struct graph_action;
}

class GraphInterface
{
public:
    GraphInterface();
    GraphInterface(const GraphInterface& g, bool keep_ref,
                   python::object ovprops, python::object oeprops);
    ~GraphInterface();

    // useful enums

    typedef enum
    {
        IN_DEGREE,
        OUT_DEGREE,
        TOTAL_DEGREE
    } degree_t;

    // general "degree" type, i.e., either a degree_t above or a string
    // representing a scalar vertex property
    typedef boost::variant<degree_t, boost::any> deg_t;

    //
    // Basic manipulation
    //

    size_t GetNumberOfVertices();
    size_t GetNumberOfEdges();
    void SetDirected(bool directed) {_directed = directed;}
    bool GetDirected() {return _directed;}
    void SetReversed(bool reversed) {_reversed = reversed;}
    bool GetReversed() {return _reversed;}

    // graph filtering
    void SetVertexFilterProperty(boost::any prop, bool invert);
    bool IsVertexFilterActive() const;
    void SetEdgeFilterProperty(boost::any prop, bool invert);
    bool IsEdgeFilterActive() const;

    // graph modification
    void InsertPropertyMap(string name, boost::any map);
    void ReIndexEdges();
    void PurgeVertices(boost::any old_index); // removes filtered vertices
    void PurgeEdges();    // removes filtered edges
    void Clear();
    void ClearEdges();
    void ShiftVertexProperty(boost::any map, size_t index) const;
    void ReIndexVertexProperty(boost::any map, boost::any old_index) const;
    void CopyVertexProperty(const GraphInterface& src, boost::any prop_src,
                            boost::any prop_tgt);
    void CopyEdgeProperty(const GraphInterface& src, boost::any prop_src,
                          boost::any prop_tgt);

    //
    // python interface
    //
    python::object DegreeMap(string deg) const;

    // used for graph properties
    graph_property_tag GetDescriptor() const { return graph_property_tag(); }
    bool CheckValid() const {return true;}

    // I/O
    void WriteToFile(string s, python::object pf, string format,
                     python::list properties);
    python::tuple ReadFromFile(string s, python::object pf, string format);

    //
    // Internal types
    //

    // the following defines the edges' internal properties
    typedef property<edge_index_t, size_t> EdgeProperty;

    // this is the main graph type
    typedef adjacency_list <vecS, // edges
                            vecS, // vertices
                            bidirectionalS,
                            no_property,
                            EdgeProperty>  multigraph_t;
    typedef graph_traits<multigraph_t>::vertex_descriptor vertex_t;
    typedef graph_traits<multigraph_t>::edge_descriptor edge_t;

    typedef property_map<multigraph_t,vertex_index_t>::type vertex_index_map_t;
    typedef property_map<multigraph_t,edge_index_t>::type edge_index_map_t;
    typedef ConstantPropertyMap<size_t,graph_property_tag> graph_index_map_t;

    // internal access

    multigraph_t& GetGraph() {return _state->_mg;}
    vertex_index_map_t GetVertexIndex() {return _vertex_index;}
    edge_index_map_t   GetEdgeIndex()   {return _edge_index;}
    size_t             GetMaxEdgeIndex(){return _state->_max_edge_index;}

    graph_index_map_t  GetGraphIndex()  {return graph_index_map_t(0);}

    void           AddEdgeIndex(const edge_t& e);
    void           RemoveEdgeIndex(const edge_t& e);

    // Gets the encapsulated graph view. See graph_filtering.cc for details
    boost::any GetGraphView() const;

private:

    // Generic graph_action functor. See graph_filtering.hh for details.
    template <class Action, class GraphViews, class Wrap, class TR1, class TR2,
              class TR3, class TR4>
    friend struct detail::graph_action;

    // python interface
    friend class PythonVertex;
    template <class Graph>
    friend class PythonEdge;

    struct state_t
    {
        // this is the main graph
        multigraph_t _mg;

        // keep track of the number of edges, since num_edges() is O(V) in
        // adjacency_list... :-(
        size_t _nedges;

        deque<size_t> _free_indexes; // indexes of deleted edges to be used up
                                     // for new edges to avoid very large
                                     // indexes, and property map memory usage
        size_t _max_edge_index;
    };

    shared_ptr<state_t> _state;

    // vertex index map
    vertex_index_map_t _vertex_index;

    // edge index map
    edge_index_map_t _edge_index;

    // this will hold an instance of the graph views at run time
    vector<boost::any> _graph_views;

    // reverse and directed states
    bool _reversed;
    bool _directed;

    // graph index map
    graph_index_map_t _graph_index;

    // vertex filter
    typedef unchecked_vector_property_map<uint8_t,vertex_index_map_t>
        vertex_filter_t;
    vertex_filter_t _vertex_filter_map;
    bool _vertex_filter_invert;
    bool _vertex_filter_active;

    // edge filter
    typedef unchecked_vector_property_map<uint8_t,edge_index_map_t>
        edge_filter_t;
    edge_filter_t _edge_filter_map;
    bool _edge_filter_invert;
    bool _edge_filter_active;
};

} //namespace graph_tool

#endif
