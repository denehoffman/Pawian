#pragma once

#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"
#include "ErrLogger/ErrLogger.hh"

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>
#include <boost/spirit/include/qi_char_.hpp>
#include <boost/spirit/include/qi_lexeme.hpp>

#include <iostream>
#include <fstream>
#include <string>

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

namespace particleGrammar
{
  struct types_ : qi::symbols<char, ParticleType> {
    types_() {
      add
	("Parton", parton)
	("Lepton", lepton)
	("Boson", boson)
	("Meson", meson)
	("Baryon", baryon)
	("Nucleus", nucleus)
	("Special", special)
	;
    }
  } types;
 
  struct functions_ : qi::symbols<char, DynFunctionType> {
    functions_() {
      add
	("BW", BW)
	("relBW", relBW)
	("relBWBlattWK", relBWBlattWK)
	("undefined", undefined)
	;
    }  
  } functions;
  
  template <typename Iterator>
  bool parse_particle(Iterator first, Iterator last, ParticleData& pData) {
    using boost::spirit::qi::double_;
    using boost::spirit::qi::_1;
    using boost::spirit::qi::phrase_parse;
    using boost::spirit::ascii::space;
    using boost::phoenix::ref;
    using qi::char_;
    using qi::int_;
    using qi::lexeme;
    
    bool wasComment = false;
    double mass=0.0, plusDm=0.0, minusDm=0.0;
    double width=0.0, plusDw=0.0, minusDw=0.0;
    
    bool r = phrase_parse(first, last,
			  
	  //  Begin grammar
	  (
	     '*' >> +char_[ref(wasComment)=true] || char_('*')[ref(wasComment)=true]
	     || '#' >> +char_[ref(wasComment)=true] || char_('#')[ref(wasComment)=true]
	     ||
	     types[ref(pData.type) = _1] >>
	     lexeme[+char_('!', 'z')[boost::phoenix::ref(pData.name) += _1]] >>
	     lexeme[+char_('!', '}')[boost::phoenix::ref(pData.texName) += _1]] >>
	     int_[ref(pData.charge) = _1] >>
	     int_[ref(pData.twoJ) = _1] >>
 	     int_[ref(pData.theParity) = _1] >> // P
 	     int_[ref(pData.theCParity) = _1] >> // C
 	     int_    // G
           ),
           //  End grammar

	   space);


    if(!wasComment && first != last) {
      r = phrase_parse(first, last,
        //  Begin grammar
        (
          int_[ref(pData.iso) = _1] >>
	  int_[ref(pData.isoThree) = _1] >>
	  int_[ref(pData.strange) = _1] >>
	  int_[ref(pData.charm) = _1] >>
	  double_[ref(mass) = _1] >> -(char_('(') >> double_[ref(plusDm) = _1] >> double_[ref(minusDm) = _1] >> char_(')')) >>
	  double_[ref(width) = _1] >> -(char_('(') >> double_[ref(plusDw) = _1] >> double_[ref(minusDw) = _1] >> char_(')')) >>
	  //	     double_[ref(pData.mass) = _1] >>
	  //	     double_[ref(pData.width).mean(_1)] >>
	  functions[ref(pData.dynamicFunction) = _1] >>
	  lexeme[+char_]
	 ),
         //  End grammar
         space);
    }

    if (wasComment) return false;
    if (!r || first != last) { // fail if we did not get a full match
      std::string rest(first, last);
      ErrMsg << "failed parsing at: " << rest << endmsg;
      return false;
    }

    // looks like success
    pData.mass.mean(mass); pData.mass.plusErr(plusDm); 
    pData.mass.minusErr(fabs(minusDm));
    pData.width.mean(width); pData.width.plusErr(plusDw); 
    pData.width.minusErr(fabs(minusDw));

    return r;
  }
}
