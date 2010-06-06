#ifndef SETUPGRAMMAR_HH
#define SETUPGRAMMAR_HH

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/variant/recursive_variant.hpp>
#include <boost/foreach.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>


namespace setupGrammar
{
  namespace fusion = boost::fusion;
  namespace phoenix = boost::phoenix;
  namespace qi = boost::spirit::qi;
  namespace ascii = boost::spirit::ascii;
  
  ///////////////////////////////////////////////////////////////////////////
  //  Our mini DECAY tree representation
  ///////////////////////////////////////////////////////////////////////////

  struct decay_tree;
  
  
}

// We need to tell fusion about our setup struct
// to make it a first-class fusion citizen

BOOST_FUSION_ADAPT_STRUCT(
    setupGrammar::decay_tree,
    (std::string, name)
    (std::vector<setupGrammar::decay_node>, children)
    (std::vector<std::string>, addParticle)
    (std::vector<std::string>, cloneParticle)
    (std::vector<std::string>, modParticle)
    (std::vector<std::string>, beamInput)
    (std::vector<std::string>, mcInput)
    (std::vector<std::string>, defineTuple)
    (std::vector<std::string>, fitVars)
    (std::vector<std::string>, initialProps)
)

namespace setupGrammar
{
    ///////////////////////////////////////////////////////////////////////////
    //  Print out the mini decay tree
    ///////////////////////////////////////////////////////////////////////////
    int const tabsize = 4;

    void tab(int indent)
    {
        for (int i = 0; i < indent; ++i)
            std::cout << ' ';
    }

    struct decay_tree_printer
    {
        decay_tree_printer(int indent = 0)
          : indent(indent)
        {
        }

        void operator()(decay_tree const& decay) const;

        int indent;
    };

    struct decay_node_printer : boost::static_visitor<>
    {
        decay_node_printer(int indent = 0)
          : indent(indent)
        {
        }

        void operator()(decay_tree const& decay) const
        {
            decay_tree_printer(indent+tabsize)(decay);
        }

        void operator()(std::string const& text) const
        {
            tab(indent+tabsize);
            std::cout << "text: \"" << text << '"' << std::endl;
        }

        int indent;
    };

    void decay_tree_printer::operator()(decay_tree const& decay) const
    {
        tab(indent);
        std::cout << "tag: " << decay.name << std::endl;
        tab(indent);
        std::cout << '{' << std::endl;

        BOOST_FOREACH(decay_node const& node, decay.children)
        {
            boost::apply_visitor(decay_node_printer(indent), node);
        }

        tab(indent);

        std::vector<std::string>::const_iterator iter;
        
        iter = decay.addParticle.begin();
        while (iter != decay.addParticle.end()) {
          std::cout << "new particle: " << *iter << std::endl;
          ++iter;
        }

        iter = decay.cloneParticle.begin();
        while (iter != decay.cloneParticle.end()) {
          std::cout << "cloned particle: " << *iter << std::endl;
          ++iter;
        }

        iter = decay.modParticle.begin();
        while (iter != decay.modParticle.end()) {
          std::cout << "modified particle: " << *iter << std::endl;
          ++iter;
        }

        std::cout << '}' << std::endl;
    }

    ///////////////////////////////////////////////////////////////////////////
    //  Our mini DECAY grammar definition
    ///////////////////////////////////////////////////////////////////////////
    //[tutorial_decay1_grammar
    template <typename Iterator>
    struct setup_file_grammar : qi::grammar<Iterator, decay_tree(), ascii::space_type>
    {
      setup_file_grammar() : setup_file_grammar::base_type(decay)
        {
          using qi::lit;
          using qi::lexeme;
          using ascii::char_;
          using ascii::string;
          using namespace qi::labels;
          
          using phoenix::at_c;
          using phoenix::push_back;
          using boost::phoenix::ref;

          restOfLine = lexeme[*(char_ -'\n') [_val += _1]];
          comment = (char_('#') >> lexeme[*(char_ -'\n')]);

          particleName = lexeme[+(char_('!','z') -',') [_val += _1]];
          decNode = (decay | particleName);
          goesTo = char_('-') >> char_('>');
          
          decay = +( comment ||
                     "addParticle"        >> restOfLine[push_back(at_c<2>(_val), _1)] || 
                     "cloneParticle"      >> restOfLine[push_back(at_c<3>(_val), _1)] || 
                     "modParticle"        >> restOfLine[push_back(at_c<4>(_val), _1)] || 
                     "beamInput"          >> restOfLine[push_back(at_c<5>(_val), _1)] || 
                     "mcInput"            >> restOfLine[push_back(at_c<6>(_val), _1)] || 
                     "defineTuple"        >> restOfLine[push_back(at_c<7>(_val), _1)] || 
                     "fitVariables"       >> restOfLine[push_back(at_c<8>(_val), _1)] || 
                     "initialProperties"  >> restOfLine[push_back(at_c<9>(_val), _1)] || 
                    ( char_('{')
                      >> (decay[push_back(at_c<1>(_val), _1)] || particleName[at_c<0>(_val) = _1])
                      >> goesTo
                      >> decNode   [push_back(at_c<1>(_val), _1)]
                      >> ','
                      >> decNode   [push_back(at_c<1>(_val), _1)]
                      >> char_('}') )
                    );
          
        }

        
        qi::rule<Iterator, decay_tree(), ascii::space_type> decay;
        qi::rule<Iterator, decay_node(), ascii::space_type> decNode;
        qi::rule<Iterator, std::string(), ascii::space_type> goesTo;
        qi::rule<Iterator, std::string(), ascii::space_type> comment;
        qi::rule<Iterator, std::string(), ascii::space_type> restOfLine;
        qi::rule<Iterator, std::string(), ascii::space_type> particleName;
    };

}

#endif
