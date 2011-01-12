#ifndef STREAMPARSER_H
#define STREAMPARSER_H

#include <iostream>
using namespace std;

class streamParser
{
  public:
    virtual void ParseStream(istream &theistream) = 0;
    virtual ~streamParser(){};
};
#endif // STREAMPARSER_H
