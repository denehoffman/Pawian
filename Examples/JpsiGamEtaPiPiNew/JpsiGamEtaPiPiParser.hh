#ifndef JpsiGamEtaPiPiParser_HH
#define JpsiGamEtaPiPiParser_HH


#include "PwaUtils/ParserBase.hh"

// Boost headers go here



class JpsiGamEtaPiPiParser : public ParserBase 
{

  public:

  JpsiGamEtaPiPiParser(int argc,char **argv);
  virtual ~JpsiGamEtaPiPiParser(){;}
  
  const bool massIndependentFit() const {return _massIndependentFit; }
  const bool useCommonProductionPhases() const {return _useCommonProductionPhases; }
  const std::pair<double, double> massRange() const { return std::make_pair( _massMin, _massMax  ) ; }
  
  
protected:
  bool _massIndependentFit;
  bool _useCommonProductionPhases;
  
  double _massMin;
  double _massMax;

};


#endif /* JpsiGamEtaPiPiParser_HH */
