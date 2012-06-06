#pragma once

#include <vector>
#include <string>

#include <boost/config.hpp> // put this first to suppress some VC++ warnings
#include <boost/variant/recursive_variant.hpp>

class ParticleTable;

namespace decayGraph {
  struct EdgeList;
}

namespace setupGrammar {

  struct decay_tree;

  typedef
  boost::variant<
    boost::recursive_wrapper<decay_tree>
    , std::string
    >
  decay_node;

  struct decay_tree {
    std::string name;                        // mother particle name
    std::vector<decay_node> children;        // children
    std::vector<std::string> addParticle;
    std::vector<std::string> cloneParticle;
    std::vector<std::string> modParticle;
    std::vector<std::string> beamInput;
    std::vector<std::string> mcInput;
    std::vector<std::string> defineTuple;
    std::vector<std::string> fitVars;
    std::vector<std::string> initialProps;
  };

}

class SetupParser
{
public:

  SetupParser();
  ~SetupParser();

  bool parse(std::string& fileName, ParticleTable* pdtTable = 0);
  const setupGrammar::decay_tree* setup() const;
  const decayGraph::EdgeList* edgeList() const;

private:
  
  setupGrammar::decay_tree* thisDecay;
};
