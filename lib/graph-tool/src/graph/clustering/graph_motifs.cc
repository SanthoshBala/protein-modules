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

#include "graph_filtering.hh"

#include "graph.hh"
#include "graph_selectors.hh"
#include "graph_properties.hh"
#include "graph_util.hh"

#include "graph_motifs.hh"

#include <boost/python.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;


struct append_to_list
{

    template <class Graph>
    void operator()(Graph& g, boost::any& list) const
    {
        typedef typename mpl::if_<typename is_directed::apply<Graph>::type,
                                  d_graph_t,
                                  u_graph_t>::type graph_t;
        vector<graph_t>& glist = any_cast<vector<graph_t>&>(list);
        glist.push_back(graph_t());
        graph_copy(g, glist.back());
    }
};

struct retrieve_from_list
{
    template <class Graph>
    void operator()(Graph& g, boost::any& list, bool& done) const
    {
        typedef typename mpl::if_<typename is_directed::apply<Graph>::type,
                                  d_graph_t,
                                  u_graph_t>::type graph_t;
        vector<graph_t>& glist = any_cast<vector<graph_t>&>(list);
        if (glist.empty())
        {
            done = true;
            return;
        }
        graph_copy(glist.back(), g);
        glist.pop_back();
    }
};

void get_motifs(GraphInterface& g, size_t k, python::list subgraph_list,
                python::list hist, python::list p, bool comp_iso,
                bool fill_list, rng_t& rng)
{
    boost::any list;
    if (g.GetDirected())
        list = vector<d_graph_t>();
    else
        list = vector<u_graph_t>();
    try
    {
        for (int i = 0; i < python::len(subgraph_list); ++i)
        {
            GraphInterface& sub =
                python::extract<GraphInterface&>(subgraph_list[i]);
            run_action<>()(sub, bind<void>(append_to_list(), _1,
                                           boost::ref(list)))();
        }
    }
    catch (bad_any_cast&)
    {
        throw ValueException("All motif graphs must be either directed or "
                             "undirected!");
    }

    vector<size_t> phist;
    vector<double> plist;
    double total = 1;
    for (int i = 0; i < python::len(p); ++i)
    {
        plist.push_back(python::extract<double>(p[i]));
        total *= plist[i];
    }

    boost::any sampler;
    if (total == 1.0)
        sampler = sample_all();
    else
        sampler = sample_some(plist, rng);

    run_action<>()
        (g, boost::bind<void>(get_all_motifs(), _1, k, boost::ref(list),
                              boost::ref(phist), _2,
                              plist[0], comp_iso, fill_list, boost::ref(rng)),
         mpl::vector<sample_all,sample_some>())(sampler);

    for (size_t i = 0; i < phist.size(); ++i)
        hist.append(phist[i]);

    if (fill_list)
    {
        for (int i = 0; i < python::len(subgraph_list); ++i)
            subgraph_list.pop();

        bool done = false;
        while (!done)
        {

            GraphInterface sub;
            sub.SetDirected(g.GetDirected());
            typedef graph_tool::detail::get_all_graph_views::apply
                <graph_tool::detail::scalar_pairs,
                mpl::bool_<false>,mpl::bool_<false>,
                mpl::bool_<false>,mpl::bool_<true>,
                mpl::bool_<true> >::type gviews;
            run_action<gviews>()
                (sub, boost::bind<void>(retrieve_from_list(), _1,
                                        boost::ref(list), boost::ref(done)))();
            if (!done)
            {
                sub.ReIndexEdges();
                subgraph_list.append(sub);
            }
        }
        subgraph_list.reverse();
    }
}
