// In this example we serialize an object and construct another from
// the XML representation. We do this twice, once directly through the
// derived class and once through an abstract base class.
//
// Compile with a command similar to
// g++ -o serializationDemo -I /opt/boost143/include/ -L /opt/boost143/lib/ -lboost_serialization -lboost_system serializationDemo.C
//
// Author: Ruediger Berlich

// Standard headers
#include <iostream>
#include <sstream>
#include <string>

// Boost headers required for serialization
#include <boost/shared_ptr.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/tracking.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/export.hpp>

// PWA test includes
#include "Examples/MATpbarpToOmegaPi/serSpin.hh"
#include "Examples/MATpbarpToOmegaPi/pbarpStates.hh"
#include "Examples/MATpbarpToOmegaPi/pbarpToOmegaPi0States.hh"

#include "Examples/MATpbarpToOmegaPi/GOmegaPiIndividual.hh"
#include "Examples/MATpbarpToOmegaPi/OmegaPiEventList.hh"
#include "Examples/EtacToapi0Fit/EtacToapi0Lh.hh"
//#include "Examples/MATpbarpToOmegaPi/OmegaPiHist.hh"

#include "Examples/MATpbarpToOmegaPi/OmegaPiData.hh"

//#include "Setup/PwaEnv.hh"
#include "Particle/ParticleTable.hh"
#include "Particle/Particle.hh"
#include "Event/EventList.hh"
#include "Event/Event.hh"
#include "Event/CBElsaReader.hh"
#include "Particle/PdtParser.hh"

#include "ErrLogger/ErrLogger.hh"

using namespace Gem::GenEvA;

/****************************************************************/

class serializationBaseClass 
{
  ////////////////////////////////////////////////////////////////
  friend class boost::serialization::access;

  template<typename Archive>
  void serialize(Archive & ar, const unsigned int)  {
    using boost::serialization::make_nvp;

    // Trigger serialization of the local data
    ar & BOOST_SERIALIZATION_NVP(testVar_);
  }
  ////////////////////////////////////////////////////////////////

public:
  serializationBaseClass() 
  { /* nothing */ }

  serializationBaseClass(int testVar): testVar_(testVar) 
  { /* nothing */ }

  virtual ~serializationBaseClass()  
  { /* nothing */ }

  int getSecret() const { return testVar_; }

  // The base class is purely virtual so it cannot be instantiated
  virtual void ImPurelyVirtual() = 0;

private:
  int testVar_;
};

// Let the Boost serialization library know that this class is abstract / purely virtual
BOOST_SERIALIZATION_ASSUME_ABSTRACT(serializationBaseClass)

/****************************************************************/

class serializationClass 
  :public serializationBaseClass
{
  ////////////////////////////////////////////////////////////////
  friend class boost::serialization::access;

  template<typename Archive>
  void serialize(Archive & ar, const unsigned int) {
    using boost::serialization::make_nvp;
    
    // Make sure the base class (and its data) gets serialized
    ar & make_nvp("serializationBaseClass", boost::serialization::base_object<serializationBaseClass>(*this));
  }
  ////////////////////////////////////////////////////////////////

public:
  serializationClass() 
  { /* nothing */ }
  
  serializationClass(int testVar) : serializationBaseClass(testVar)  
  { /* nothing */ }
  
  virtual ~serializationClass() 
  { /* nothing */ }

  // Make sure we can instantiate this class -- need to re-implement the
  // purely virtual function of the base class.
  virtual void ImPurelyVirtual() { /* nothing */ }
};

// Register this class with the boost serialization library. 
// Note: This macro and the associated header would usually reside in 
// the implementation file (.cpp), not in the header file.
#include <boost/serialization/export.hpp>
BOOST_CLASS_EXPORT(serializationClass)

/****************************************************************/

int main() {
  // The objects that should be (de-)serialized
  boost::shared_ptr<serializationClass> p_derived(new serializationClass(123));
  boost::shared_ptr<serializationBaseClass> p_base(new serializationClass(456));

  // Will hold reconstructed objects. Here these pointers are empty
  boost::shared_ptr<serializationClass> p_derived_reconstructed;
  boost::shared_ptr<serializationBaseClass> p_base_reconstructed;

  // The PWA-objects that should be (de-)serialized and changed inbetween
  //Spin
  boost::shared_ptr<serSpin> mySpin_init(new serSpin());
  boost::shared_ptr<serSpin> mySpin_reco1;
  boost::shared_ptr<serSpin> mySpin_reco2;
  //pbarpStates
  boost::shared_ptr<pbarpStates> myStates_init(new pbarpStates(2));
  boost::shared_ptr<pbarpStates> myStates_reco1;
  boost::shared_ptr<pbarpStates> myStates_reco2;
  //pbarpToOmegaStates
  boost::shared_ptr<pbarpToOmegaPi0States> myOmStates_init(new pbarpToOmegaPi0States(myStates_init));
  boost::shared_ptr<pbarpToOmegaPi0States> myOmStates_reco1;
  boost::shared_ptr<pbarpToOmegaPi0States> myOmStates_reco2;
  //OmegaPiEvent List / GOmegaPiIndividual
  std::string theSourcePath=getenv("CMAKE_SOURCE_DIR");
  
  std::string piomegaDatFile;
  std::string piomegaMcFile; 
  piomegaDatFile=theSourcePath+"/Examples/MATpbarpToOmegaPi/data/510_0600.dat";
  piomegaMcFile=theSourcePath+"/Examples/MATpbarpToOmegaPi/data/mc510_0600.dat";
  ParticleTable pTable;
  PdtParser parser;
  std::string pdtFile(theSourcePath+"/Particle/pdt.table");
  if (!parser.parse(pdtFile, pTable)) {
    Alert << "Error: could not parse " << pdtFile << endmsg;
    exit(1);
  }
  std::vector<std::string> fileNames;

  fileNames.push_back(piomegaDatFile);
  CBElsaReader eventReader(fileNames, 3, 0); 
  EventList piOmegaEventsData;
  eventReader.fillAll(piOmegaEventsData);

  if (!piOmegaEventsData.findParticleTypes(pTable))
    Warning << "could not find all particles" << endmsg;

  Info << "\nFile has " << piOmegaEventsData.size() << " events. Each event has "
                  <<  piOmegaEventsData.nextEvent()->size() << " final state particles.\n" << endmsg;
  piOmegaEventsData.rewind();


  Event* anEvent;
  int evtCount = 0;
  while ((anEvent = piOmegaEventsData.nextEvent()) != 0 && evtCount < 20) {
    Info << "\n" 
                    << *(anEvent->p4(0)) << "\tm = " << anEvent->p4(0)->Mass() << "\n"
                    << *(anEvent->p4(1)) << "\tm = " << anEvent->p4(1)->Mass() << "\n"
                    << *(anEvent->p4(2)) << "\tm = " << anEvent->p4(2)->Mass() << "\n"
                    << endmsg;
    ++evtCount;
  }
  piOmegaEventsData.rewind();

  std::vector<std::string> fileNamesMc;

  fileNamesMc.push_back(piomegaMcFile);
  CBElsaReader eventReaderMc(fileNamesMc, 3, 0); 
  EventList piOmegaEventsMc;
  eventReaderMc.fillAll(piOmegaEventsMc);
  piOmegaEventsMc.rewind();


  boost::shared_ptr<pbarpStates> pbarpStatesPtr(new pbarpStates(1));
  boost::shared_ptr<const pbarpToOmegaPi0States> pbarpToOmegaPi0StatesPtr(new pbarpToOmegaPi0States(pbarpStatesPtr));
  //boost::shared_ptr<const OmegaPiEventList> theOmegaPiEventPtr(new OmegaPiEventList(piOmegaEventsData, piOmegaEventsMc, 1,  600));
  boost::shared_ptr<OmegaPiEventList> theOmegaPiEventPtr = OmegaPiEventList::getList();
  theOmegaPiEventPtr->initList(piOmegaEventsData, piOmegaEventsMc, 1,  600);
  //boost::shared_ptr<const OmegaPiEventList> theOmegaPiEventPtr_reco1;
 // boost::shared_ptr<const OmegaPiEventList> theOmegaPiEventPtr_reco2;
  boost::shared_ptr<GOmegaPiIndividual> gdii_ptr_init(new GOmegaPiIndividual(pbarpToOmegaPi0StatesPtr));
  boost::shared_ptr<GOmegaPiIndividual> gdii_ptr_reco1;
  boost::shared_ptr<GOmegaPiIndividual> gdii_ptr_reco2;

  //////////////////////////////////////////////////////////////////////////////////
  std::ostringstream oarchive_stream_pwa1;
  // Serialization (through boost::shared_ptr)
  {
    boost::archive::xml_oarchive oa(oarchive_stream_pwa1);
    oa << boost::serialization::make_nvp("mySpin", gdii_ptr_init);
  } // Explicit scope necessary -- makes sure oa's destructor is called
  
  // Here is the serialized data
  std::string spin_init_ser = oarchive_stream_pwa1.str();

  // We now have an XML representation of the object. Print it out just to look at it
  std::cout << "=============================================================" << std::endl
            << "Class representation, when serializing initial OmegaPiEventList class:" << std::endl
            << std::endl
            << spin_init_ser << std::endl;

  // Feed the data into an istringstream object
  std::istringstream iarchive_stream_pwa1(spin_init_ser);

  // De-Serialization
  {
    boost::archive::xml_iarchive ia(iarchive_stream_pwa1);
    ia >> boost::serialization::make_nvp("mySpin", gdii_ptr_reco1);
  } // note: explicit scope here is essential so the ia-destructor gets called

  //make changes (pretend being network client)
  boost::shared_ptr<pbarpStates> pbarpStatesPtr_new(new pbarpStates(1));
  boost::shared_ptr<const pbarpToOmegaPi0States> pbarpToOmegaPi0StatesPtr_new(new pbarpToOmegaPi0States(pbarpStatesPtr_new));
  //boost::shared_ptr<const OmegaPiEventList> theOmegaPiEventPtr_new(new OmegaPiEventList(piOmegaEventsData, piOmegaEventsMc, 2,  600));
  boost::shared_ptr<OmegaPiEventList> theOmegaPiEventPtr_new = OmegaPiEventList::getList();
  //theOmegaPiEventPtr_new->initList(piOmegaEventsData, piOmegaEventsMc, 2,  600);
  boost::shared_ptr<GOmegaPiIndividual> gdii_ptr_new(new GOmegaPiIndividual(pbarpToOmegaPi0StatesPtr_new));
  gdii_ptr_reco1=gdii_ptr_new;

  std::ostringstream oarchive_stream_pwa2;
  // Serialization (through boost::shared_ptr)
  {
    boost::archive::xml_oarchive oa(oarchive_stream_pwa2);
    oa << boost::serialization::make_nvp("mySpin", gdii_ptr_reco1);
  } // Explicit scope necessary -- makes sure oa's destructor is called
  
  // Here is the serialized data
  std::string spin_changed_ser = oarchive_stream_pwa2.str();

  // We now have an XML representation of the object. Print it out just to look at it
  std::cout << "=============================================================" << std::endl
            << "Class representation, when serializing changed OmegaPiEventList class:" << std::endl
            << std::endl
            << spin_changed_ser << std::endl;

  // Feed the data into an istringstream object
  std::istringstream iarchive_stream_pwa2(spin_changed_ser);

  // De-Serialization through derived object
  {
    boost::archive::xml_iarchive ia(iarchive_stream_pwa2);
    ia >> boost::serialization::make_nvp("mySpin", gdii_ptr_reco2);
  } // note: explicit scope here is essential so the ia-destructor gets called

  // Print testVar_ of original and reconstructed class
  OmegaPiData::fitParamVal myFitParmVal;
  std::cout << std::endl;
  std::cout << " Init: " << gdii_ptr_init->getFitParams(myFitParmVal) << std::endl; 
  //theOmegaPiEventPtr_init->print(cout);
  std::cout << std::endl << " Reco1: " << gdii_ptr_reco1->getFitParams(myFitParmVal) << std::endl; 
  //theOmegaPiEventPtr_reco1->print(cout);
  std::cout << std::endl << " Reco2: " << gdii_ptr_reco2->getFitParams(myFitParmVal) << std::endl; 
  //theOmegaPiEventPtr_reco2->print(cout);
  std::cout << std::endl;


  //////////////////////////////////////////////////////////////////////////////////
  std::ostringstream oarchive_stream_derived;
  // Serialization through derived class (through boost::shared_ptr)
  {
    boost::archive::xml_oarchive oa(oarchive_stream_derived);
    oa << boost::serialization::make_nvp("classHierarchyFromDerived", p_derived);
  } // Explicit scope necessary -- makes sure oa's destructor is called
  
  // Here is the serialized data
  std::string derived_class_serialized = oarchive_stream_derived.str();

  // We now have an XML representation of the object. Print it out just to look at it
  std::cout << "=============================================================" << std::endl
            << "Class representation, when serializing the derived class:" << std::endl
            << std::endl
            << derived_class_serialized << std::endl;

  // Feed the data into an istringstream object
  std::istringstream iarchive_stream_derived(derived_class_serialized);

  // De-Serialization through derived object
  {
    boost::archive::xml_iarchive ia(iarchive_stream_derived);
    ia >> boost::serialization::make_nvp("classHierarchyFromDerived", p_derived_reconstructed);
  } // note: explicit scope here is essential so the ia-destructor gets called

  // Print testVar_ of original and reconstructed class
  std::cout << std::endl 
	    << "Before: " << p_derived->getSecret() << " After: " << p_derived_reconstructed->getSecret() << std::endl;


  //////////////////////////////////////////////////////////////////////////////////
  std::ostringstream oarchive_stream_base;
  // Serialization through a base pointer to purely abstract class
  {
    boost::archive::xml_oarchive oa(oarchive_stream_base);
    oa << boost::serialization::make_nvp("classHierarchyFromBasePointer", p_base);
  } // Explicit scope necessary -- makes sure oa's destructor is called

  // Here is the serialized data
  std::string base_class_serialized = oarchive_stream_base.str();

  // We now have an XML representation of the object. Print it out just to look at it
  std::cout << "=============================================================" << std::endl
            << "Class representation, when serialized through a base pointer" << std::endl
            << std::endl
            << base_class_serialized << std::endl;

  // Feed the data in the istringstream
  std::istringstream iarchive_stream_base(base_class_serialized);

  // De-Serialization through a base pointer
  {
    boost::archive::xml_iarchive ia(iarchive_stream_base);
    ia >> boost::serialization::make_nvp("classHierarchyFromBasePointer", p_base_reconstructed);
  } // note: explicit scope here is essential so the ia-destructor gets called

  // Print testVar_ of original and reconstructed class
  std::cout << std::endl 
	    << "Before: " << p_base->getSecret() << " After: " << p_base_reconstructed->getSecret() << std::endl;


//  Destructor's correct?
 // while(true) {
 //  boost::shared_ptr<pbarpStates> pbarpStatesPtr_test(new pbarpStates(2));
 //   boost::shared_ptr<const pbarpToOmegaPi0States> pbarpToOmegaPi0StatesPtr_test(new pbarpToOmegaPi0States(pbarpStatesPtr_test));
 //   boost::shared_ptr<const OmegaPiEventList> theOmegaPiEventPtr_test(new OmegaPiEventList(piOmegaEventsData, piOmegaEventsMc, 2,  600));
 //   boost::shared_ptr<GOmegaPiIndividual> gdii_ptr_test(new GOmegaPiIndividual(theOmegaPiEventPtr_test, pbarpToOmegaPi0StatesPtr_test));
 // }


}

/****************************************************************/
