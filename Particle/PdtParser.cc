#include "Particle/PdtParser.hh"
// Particle.hh defines ParticleData struct
#include "Particle/Particle.hh"
#include "Particle/ParticleTable.hh"

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


namespace client
{

    struct types_ : qi::symbols<char, ParticleType>
    {
        types_()
        {
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
 

    template <typename Iterator>
    bool parse_particle(Iterator first, Iterator last, ParticleData& pData)
    {
        using boost::spirit::qi::double_;
        using boost::spirit::qi::_1;
        using boost::spirit::qi::phrase_parse;
        using boost::spirit::ascii::space;
        using boost::phoenix::ref;
	using qi::char_;
	using qi::int_;
	using qi::lexeme;

	bool wasComment = false;

        bool r = phrase_parse(first, last,

            //  Begin grammar
            (
	     '*' >> +char_[ref(wasComment)=true] || char_('*')[ref(wasComment)=true]
	     || '#' >> +char_[ref(wasComment)=true] || char_('#')[ref(wasComment)=true]
	     ||
	     types[ref(pData.type) = _1] >>
	     lexeme[+char_('!', 'z')[ref(pData.name) += _1]] >>
	     lexeme[+char_('!', '}')[ref(pData.texName) += _1]] >>
	     int_[ref(pData.charge) = _1] >>
	     int_[ref(pData.twoJ) = _1] >>
 	     int_ >> // P
 	     int_ >> // C
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
	     double_[ref(pData.mass) = _1] >>
	     double_[ref(pData.width) = _1] >>
	     lexeme[+char_]
            ),
            //  End grammar
            space);
	}

	if (wasComment) return false;
        if (!r || first != last) // fail if we did not get a full match
            return false;

        return r;
    }
}



PdtParser::PdtParser()
{
}

PdtParser::~PdtParser()
{
}

bool PdtParser::parse(std::string& fileName, ParticleTable& table)
{
  ParticleData* pData;

  std::ifstream file;
  file.open(fileName.c_str());

  if (file.is_open()) {
    std::string str;
    while (!file.eof()) {
      std::getline(file, str);
      if (!file.eof()) {
	pData = new ParticleData;
	if (client::parse_particle(str.begin(), str.end(), *pData)) { // success
	  Particle* newParticle = new Particle(*pData);
	  table.addParticle(newParticle);
	}
// 	else {
// 	  std::cout << "failed: " << str << std::endl;
// 	}
      }
    }
    file.close();
    return true; // success
  }
  else
    return false;
}

