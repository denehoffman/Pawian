/*
 * This file is part of Pawian, which is distributed under
 * GNU General Public License v3 (https://www.gnu.org/licenses/gpl-3.0.html)
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE EventReader

#include "Event/Event.hh"
#include "Event/EventList.hh"
#include "Event/EventReaderDefault.hh"
#include "qft++/topincludes/tensor.hh"

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(TestEventReaderDefault);

BOOST_AUTO_TEST_CASE(test_fill) {
  std::vector<std::string> inputFiles = {"Event/test/SampleFileNoHeader.dat"};
  int numberOfParticles = 3;
  int skipFirstLines = 0;
  bool hasWeights = true;
  EventReaderDefault gevEventReader(inputFiles, numberOfParticles,
                                    skipFirstLines, hasWeights);
  EventReaderDefault mevEventReader(inputFiles, numberOfParticles,
                                    skipFirstLines, hasWeights);
  EventList gevEvents, mevEvents;

  gevEventReader.setUnit("GEV");
  gevEventReader.setOrder("Px Py Pz E");
  gevEventReader.fill(gevEvents);

  mevEventReader.setUnit("MEV");
  mevEventReader.setOrder("E Px Py Pz"); // ! note mock order
  mevEventReader.fill(mevEvents);

  BOOST_CHECK_EQUAL(gevEvents.size(), 2);
  BOOST_CHECK_EQUAL(mevEvents.size(), 2);
  if (gevEvents.size() != mevEvents.size())
    throw std::runtime_error("Same number of events required to continue");

  std::vector<double> expectedWeights = {.3, .5};
  std::vector<float> gevParticle1Energy = {1.41796, 1.38767};
  std::vector<float> gevParticle2Pz = {0.00401993, -0.416389};
  std::vector<float> mevParticle2Energy = {-0.185839, -0.127614};
  std::vector<float> mevParticle3Pz = {0.790078, 1.19683};
  for (int i = 0; i < gevEvents.size(); ++i) {
    auto gevEvent = gevEvents.nextEvent();
    auto mevEvent = mevEvents.nextEvent();
    if (!gevEvent || !mevEvent)
      throw std::runtime_error("Nullpointer event!");
    BOOST_CHECK_EQUAL(gevEvent->size(), 3);
    BOOST_CHECK_EQUAL(mevEvent->size(), 3);
    BOOST_CHECK_EQUAL(mevEvent->Weight(), expectedWeights.at(i));
    BOOST_CHECK_EQUAL(gevEvent->Weight(), expectedWeights.at(i));
    auto gevParticle1 = gevEvent->p4(0);
    auto gevParticle2 = gevEvent->p4(1);
    auto mevParticle2 = mevEvent->p4(1);
    auto mevParticle3 = mevEvent->p4(2);
    BOOST_CHECK_CLOSE(gevParticle1->E(), gevParticle1Energy.at(i), 1e-5);
    BOOST_CHECK_CLOSE(gevParticle2->Pz(), gevParticle2Pz.at(i), 1e-5);
    BOOST_CHECK_CLOSE(mevParticle2->E(), mevParticle2Energy.at(i) / 1e3, 1e-5);
    BOOST_CHECK_CLOSE(mevParticle3->Pz(), mevParticle3Pz.at(i) / 1e3, 1e-5);
  }
}

BOOST_AUTO_TEST_SUITE_END();
