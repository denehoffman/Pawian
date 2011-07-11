#ifndef STREAMPARSERLS_H
#define STREAMPARSERLS_H

#include <iostream>
using namespace std;

class StreamParserLS
{
  public:
    virtual void ParseStream(istream &theistream) = 0;
    virtual ~StreamParserLS(){};
};
#endif // STREAMPARSERLS_H
