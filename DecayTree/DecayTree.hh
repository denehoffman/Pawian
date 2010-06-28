#ifndef DECAYTREE_HH
#define DECAYTREE_HH

#include <iostream>
#include <cstdlib>
#include <iterator>
#include <algorithm>
#include <string>
#include <map>

#include <boost/utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/visitors.hpp>

class ParticleTable;
class Particle;

namespace decayGraph
{
  struct Decay {
    int mother;
    std::vector<int> daughters;
    Decay(int mom) {mother = mom;};
    void print(std::ostream& o) const {
      std::vector<int>::const_iterator iter;
      o << mother << " -> ";
      for (iter = daughters.begin(); iter != daughters.end(); iter++)
	o << *iter << " ";
    o << std::endl;
    }
  };

  struct EdgeList {
    std::vector<Decay*> decays;
    std::vector<Decay*>::reverse_iterator lastVertex;
    std::map<int, std::string> particleNames;
    int lastVertexNumber;
    void print(std::ostream& o) const {
      std::vector<Decay*>::const_iterator iter;
      for (iter = decays.begin(); iter != decays.end(); iter++)
	(*iter)->print(o);
      std::map<int, std::string>::const_iterator nameIter;
      for (nameIter = particleNames.begin(); nameIter != particleNames.end(); nameIter++)
	o << nameIter->first << " " << nameIter->second << std::endl;
    }
  };

  using namespace boost;
  typedef std::pair<int, int> Edge;
  typedef adjacency_list<vecS, vecS, bidirectionalS> Graph;
  typedef graph_traits<Graph>::vertex_descriptor Vertex;
}

class DecayTree
{
public:
  DecayTree();
  ~DecayTree();
  DecayTree(const decayGraph::EdgeList* edgeList);
  
  bool fillParticleRefs(ParticleTable& ptable);
  void print(std::ostream& o) const;

private:
  decayGraph::Graph *theDecayTree;
  std::map<int, std::string> *nameMap;
  std::map<int, const Particle*> *particleMap;
};

#endif
