#ifndef SETUPPARSER_HH
#define SETUPPARSER_HH

#include <string>

class SetupParser
{
public:

  SetupParser();
  ~SetupParser();

  bool parse(std::string& fileName);

};

#endif
