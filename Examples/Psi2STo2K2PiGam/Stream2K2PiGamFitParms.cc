#include <fstream>

#include "Examples/Psi2STo2K2PiGam/Stream2K2PiGamFitParms.hh"
#include "ErrLogger/ErrLogger.hh"

Stream2K2PiGamFitParms::Stream2K2PiGamFitParms(std::string& filePath, boost::shared_ptr<const Psi2STo2K2PiGamStates> theStates) :
  AbsFitParamStreamer(filePath),
  _states(theStates)
{
  fillParamMap();

}

Stream2K2PiGamFitParms::~Stream2K2PiGamFitParms(){;}

void Stream2K2PiGamFitParms::fillParamMap(){
  
  std::vector< boost::shared_ptr<const JPCLS> >::const_iterator itJPCLS;
  StringPairMap::const_iterator stringPairIter;

  std::vector< boost::shared_ptr<const JPCLS> > JPCLSPsiToChiGam=_states->PsiToChiGamStates();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTo2K892States=_states->ChiTo2K892States();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTo2K_2_1430States=_states->ChiTo2K_2_1430States();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTo2K_0_1430States=_states->ChiTo2K_0_States();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSK1400ToKst1PiStates=_states->K1400ToKst1PiStates();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiTof0f0States=_states->ChiTof0f0States();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSChiToK0K0States=_states->ChiToK0K0States();
  std::vector< boost::shared_ptr<const JPCLS> > JPCLSK1ToK0PiStates=_states->K1ToK0PiStates();
 
  std::string theSuffix("ChiGamVal");
  fillParamMapAmps(JPCLSPsiToChiGam, theSuffix, _paramVal.PsiToChiGam);
  theSuffix="ChiGamErr";
  fillParamMapAmps(JPCLSPsiToChiGam, theSuffix, _paramErr.PsiToChiGam);

  theSuffix="K892K892Val";
  fillParamMapAmps(JPCLSChiTo2K892States, theSuffix, _paramVal.ChiTo2K892);

  theSuffix="K892K892Err";
  fillParamMapAmps(JPCLSChiTo2K892States, theSuffix, _paramErr.ChiTo2K892);

  theSuffix="K_0_1430K_0_1430Val";
  fillParamMapAmps(JPCLSChiTo2K_0_1430States, theSuffix, _paramVal.ChiTo2K_0_1430);

  theSuffix="K_0_1430K_0_1430Err";
  fillParamMapAmps(JPCLSChiTo2K_0_1430States, theSuffix, _paramErr.ChiTo2K_0_1430);

  theSuffix="K_2_1430K_2_1430Val";
  fillParamMapAmps(JPCLSChiTo2K_2_1430States, theSuffix, _paramVal.ChiTo2K_2_1430);

  theSuffix="K_2_1430K_2_1430Err";
  fillParamMapAmps(JPCLSChiTo2K_2_1430States, theSuffix, _paramErr.ChiTo2K_2_1430);

  theSuffix="K_1_1400KVal";
  fillParamMapAmps(JPCLSK1400ToKst1PiStates, theSuffix, _paramVal.K1400ToK892Pi);

  theSuffix="K_1_1400KErr";
  fillParamMapAmps(JPCLSK1400ToKst1PiStates, theSuffix, _paramErr.K1400ToK892Pi);


  theSuffix="K_1_1270ToK892PiVal";
  fillParamMapAmps(JPCLSK1400ToKst1PiStates, theSuffix, _paramVal.K_1_1270ToK892Pi);

  theSuffix="K_1_1270ToK892PiErr";
  fillParamMapAmps(JPCLSK1400ToKst1PiStates, theSuffix, _paramErr.K_1_1270ToK892Pi);


  theSuffix="K_1_1270ToK_0_1430PiVal";
  fillParamMapAmps(JPCLSK1ToK0PiStates, theSuffix, _paramVal.K_1_1270ToK_0_1430Pi);

  theSuffix="K_1_1270ToK_0_1430PiErr";
  fillParamMapAmps(JPCLSK1ToK0PiStates, theSuffix, _paramErr.K_1_1270ToK_0_1430Pi);


  theSuffix="K_0_1430K_0_1950Val";
  fillParamMapAmps(JPCLSChiTo2K_0_1430States, theSuffix, _paramVal.ChiToK_0_1430K_0_1950);

  theSuffix="K_0_1430K_0_1950Err";
  fillParamMapAmps(JPCLSChiTo2K_0_1430States, theSuffix, _paramErr.ChiToK_0_1430K_0_1950);


  theSuffix="f980_pif1710_kVal";
  fillParamMapAmps(JPCLSChiTof0f0States, theSuffix, _paramVal.ChiTof980_pif1710_k);

  theSuffix="f980_pif1710_kErr";
  fillParamMapAmps(JPCLSChiTof0f0States, theSuffix, _paramErr.ChiTof980_pif1710_k);

  theSuffix="f980_kf1710_piVal";
  fillParamMapAmps(JPCLSChiTof0f0States, theSuffix, _paramVal.ChiTof980_kf1710_pi);

  theSuffix="f980_kf1710_piErr";
  fillParamMapAmps(JPCLSChiTof0f0States, theSuffix, _paramErr.ChiTof980_kf1710_pi);

  theSuffix="f980f980Val";
  fillParamMapAmps(JPCLSChiTof0f0States, theSuffix, _paramVal.ChiTof980f980);

  theSuffix="f980f980Err";
  fillParamMapAmps(JPCLSChiTof0f0States, theSuffix, _paramErr.ChiTof980f980);

  theSuffix="f980f2200Val";
  fillParamMapAmps(JPCLSChiTof0f0States, theSuffix, _paramVal.ChiTof980f2200);

  theSuffix="f980f2200Err";
  fillParamMapAmps(JPCLSChiTof0f0States, theSuffix, _paramErr.ChiTof980f2200);

  theSuffix="f980_pif1370_kVal";
  fillParamMapAmps(JPCLSChiTof0f0States, theSuffix, _paramVal.ChiTof980_pif1370_k);

  theSuffix="f980_pif1370_kErr";
  fillParamMapAmps(JPCLSChiTof0f0States, theSuffix, _paramErr.ChiTof980_pif1370_k);


  theSuffix="f980_kf1370_piVal";
  fillParamMapAmps(JPCLSChiTof0f0States, theSuffix, _paramVal.ChiTof980_kf1370_pi);

  theSuffix="f980_kf1370_piErr";
  fillParamMapAmps(JPCLSChiTof0f0States, theSuffix, _paramErr.ChiTof980_kf1370_pi);


  theSuffix="K_0_2400KToKf980KVal";
  fillParamMapAmps(JPCLSChiToK0K0States, theSuffix, _paramVal.ChiToK_0_2400ToKf980);

  theSuffix="K_0_2400KToKf980KErr";
  fillParamMapAmps(JPCLSChiToK0K0States, theSuffix, _paramErr.ChiToK_0_2400ToKf980);


  theSuffix="K892K_1_1680Val";
  fillParamMapAmps(JPCLSChiTo2K892States, theSuffix, _paramVal.ChiToK892K1680);

  theSuffix="K892K_1_1680Err";
  fillParamMapAmps(JPCLSChiTo2K892States, theSuffix, _paramErr.ChiToK892K1680);


  theSuffix="K892K_1_2300Val";
  fillParamMapAmps(JPCLSChiTo2K892States, theSuffix, _paramVal.ChiToK892K2300);

  theSuffix="K892K_1_2300Err";
  fillParamMapAmps(JPCLSChiTo2K892States, theSuffix, _paramErr.ChiToK892K2300);


  theSuffix="sigmaf1710Val";
  fillParamMapAmps(JPCLSChiTof0f0States, theSuffix, _paramVal.ChiToSigmaf1710);

  theSuffix="sigmaf1710Err";
  fillParamMapAmps(JPCLSChiTof0f0States, theSuffix, _paramErr.ChiToSigmaf1710);


  theSuffix="K_0_1460ToK892PiVal";
  fillParamMapAmps(JPCLSChiToK0K0States, theSuffix, _paramVal.K_0_1460ToK892Pi);

  theSuffix="K_0_1460ToK892PiErr";
  fillParamMapAmps(JPCLSChiToK0K0States, theSuffix, _paramErr.K_0_1460ToK892Pi);



  theSuffix="K_0_1460ToK_0_1430PiVal";
  fillParamMapAmps(JPCLSChiToK0K0States, theSuffix, _paramVal.K_0_1460ToK_0_1430Pi);

  theSuffix="K_0_1460ToK_0_1430PiErr";
  fillParamMapAmps(JPCLSChiToK0K0States, theSuffix, _paramErr.K_0_1460ToK_0_1430Pi);


  std::string theKey("K892massVal");
  fillParamMapMass(theKey, _paramVal.BwK892);
  theKey="K892massErr";
  fillParamMapMass(theKey, _paramErr.BwK892);

  theKey="K_0_1430massVal";
  fillParamMapMass(theKey, _paramVal.BwK_0_1430);
  theKey="K_0_1430massErr";
  fillParamMapMass(theKey, _paramErr.BwK_0_1430);

  theKey="K_2_1430massVal";
  fillParamMapMass(theKey, _paramVal.BwK_2_1430);
  theKey="K_2_1430massErr";
  fillParamMapMass(theKey, _paramErr.BwK_2_1430);

  theKey="K_1_1400massVal";
  fillParamMapMass(theKey, _paramVal.BwK1400);
  theKey="K_1_1400massErr";
  fillParamMapMass(theKey, _paramErr.BwK1400);

  theKey="K_1_1270massVal";
  fillParamMapMass(theKey, _paramVal.BwK_1_1270);
  theKey="K_1_1270massErr";
  fillParamMapMass(theKey, _paramErr.BwK_1_1270);

  theKey="f980_pimassVal";
  fillParamMapMass(theKey, _paramVal.Bwf980_pi);
  theKey="f980_pimassErr";
  fillParamMapMass(theKey, _paramErr.Bwf980_pi);

  theKey="FlatMf980Val";
  fillParamFlatte(theKey, _paramVal.Flatf980);
  theKey="FlatMf980Err";
  fillParamFlatte(theKey, _paramErr.Flatf980);

  theKey="Flatf980gPiPiVal";
  fillParamFlatte(theKey, _paramVal.Flatf980gPiPi);
  theKey="Flatf980gKKErr";
  fillParamFlatte(theKey, _paramErr.Flatf980gPiPi);

  theKey="Flatf980gKKVal";
  fillParamFlatte(theKey, _paramVal.Flatf980gKK);
  theKey="Flatf980gKKErr";
  fillParamFlatte(theKey, _paramErr.Flatf980gKK);

  theKey="f1710_kmassVal";
  fillParamMapMass(theKey, _paramVal.Bwf1710_k);
  theKey="f1710_kmassErr";
  fillParamMapMass(theKey, _paramErr.Bwf1710_k);

  theKey="f2200massVal";
  fillParamMapMass(theKey, _paramVal.Bwf2200);
  theKey="f2200massErr";
  fillParamMapMass(theKey, _paramErr.Bwf2200);


  theKey="sigmamassVal";
  fillParamMapMass(theKey, _paramVal.BwSigma);
  theKey="sigmamassErr";
  fillParamMapMass(theKey, _paramErr.BwSigma);

  theKey="f1370_kmassVal";
  fillParamMapMass(theKey, _paramVal.Bwf1370_k);
  theKey="f1370_kmassErr";
  fillParamMapMass(theKey, _paramErr.Bwf1370_k);

  theKey="K_0_2400Val";
  fillParamMapMass(theKey, _paramVal.BwK_0_2400);
  theKey="K_0_2400Err";
  fillParamMapMass(theKey, _paramErr.BwK_0_2400);

  theKey="K_0_1950Val";
  fillParamMapMass(theKey, _paramVal.BwK_0_1950);
  theKey="K_0_1950Err";
  fillParamMapMass(theKey, _paramErr.BwK_0_1950);

  theKey="K_1_1680Val";
  fillParamMapMass(theKey, _paramVal.BwK_1_1680);
  theKey="K_1_1680Err";
  fillParamMapMass(theKey, _paramErr.BwK_1_1680);

  theKey="K_1_2300Val";
  fillParamMapMass(theKey, _paramVal.BwK_1_2300);
  theKey="K_1_2300Err";
  fillParamMapMass(theKey, _paramErr.BwK_1_2300);

  theKey="K_0_1460Val";
  fillParamMapMass(theKey, _paramVal.BwK_0_1460);
  theKey="K_0_1460Err";
  fillParamMapMass(theKey, _paramErr.BwK_0_1460);


  theKey="PhaseSpaceVal";
  fillParamFlatte(theKey, _paramVal.phaseSpace);
  theKey="PhaseSpaceErr";
  fillParamFlatte(theKey, _paramErr.phaseSpace);

}
