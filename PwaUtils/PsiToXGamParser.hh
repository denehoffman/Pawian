#ifndef PsiToXGamParser_HH
#define PsiToXGamParser_HH


#include "PwaUtils/ParserBase.hh"

// Boost headers go here



class PsiToXGamParser : public ParserBase 
{

  public:

  PsiToXGamParser(int argc,char **argv);
  virtual ~PsiToXGamParser(){;}
  
  const bool massIndependentFit() const {return _massIndependentFit; }
  const bool useCommonProductionPhases() const {return _useCommonProductionPhases; }
  const bool useEvtWeight() const {return _useEvtWeight; }
  const std::pair<double, double> massRange() const { return std::make_pair( _massMin, _massMax  ) ; }
  const std::string getjobOption() const {return _jobOption;}
  
protected:
  bool _massIndependentFit;
  bool _useCommonProductionPhases;
  bool _useEvtWeight;
  
  double _massMin;
  double _massMax;

  std::string _jobOption;
};


#endif /*PsiToXGamParser_HH */
