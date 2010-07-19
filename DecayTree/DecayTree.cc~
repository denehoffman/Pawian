#include "DecayTree/DecayTree.hh"
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/iteration_macros.hpp>

#include <map> 

using namespace decayGraph;

DecayTree::DecayTree()
{
}

DecayTree::~DecayTree()
{
}

DecayTree::DecayTree(const EdgeList* edgeList)
{
  Edge edges[edgeList->decays.size() * 2];
  int counter;
  for (counter = 0; counter < edgeList->decays.size(); counter++) {
    edges[2*counter]   = std::pair<int,int>(edgeList->decays[counter]->mother,
					    edgeList->decays[counter]->daughters[0]);
    edges[2*counter+1] = std::pair<int,int>(edgeList->decays[counter]->mother,
					    edgeList->decays[counter]->daughters[1]);
  }

  const std::size_t nedges = sizeof(edges)/sizeof(Edge);
  theDecayTree = new Graph(edges, edges+nedges, edgeList->lastVertexNumber);

  nameMap = &(const_cast<EdgeList*>(edgeList)->particleNames);
  for (counter = 1; counter <= nameMap->size(); ++counter)
    std::cout << 
      " name " << counter << " " << nameMap->find(counter)->second << std::endl;
}

bool DecayTree::fillParticleRefs(ParticleTable& ptable)
{
  particleMap = new std::map<int, const Particle*>;
  std::map<int, std::string>::const_iterator nameIter;
  for (nameIter = nameMap->begin(); nameIter != nameMap->end(); ++nameIter) {
    if (nameIter->second.size() > 0)
      particleMap->insert(std::pair<int, const Particle*>(nameIter->first,
							  ptable.particle(nameIter->second)));
  }
  return true;
}


void DecayTree::print(std::ostream& o) const
{
  o << "a decay" << std::endl;
  return;
}
