#include "ErrLogger.hh"

int main(int argc, char* argv[]) {

  ErrLogger::instance().setThreshold(logging::log_level::TRACE);

  InfoMsg << "An Info" << endmsg;
  DebugMsg << "Debug Message" << endmsg;
  WarningMsg << "A Warning" << endmsg;
  ErrMsg << "Error Message" << endmsg;
  Alert << "A Fatal Error" << endmsg;

  ErrLogger::instance().setThreshold(logging::log_level::ERROR);

  DebugMsg << "2nd Debug Message" << endmsg;
  InfoMsg << "2nd Info" << endmsg;
  WarningMsg << "2nd Warning" << endmsg;
  ErrMsg << "2nd Error Message" << endmsg;
  Alert << "2nd Fatal Error" << endmsg;

  return 0;
}
