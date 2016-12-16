#include "ErrLogger/ErrLogger.hh"

int main(int argc, char* argv[])
{

  ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);

  DebugMsg << "Debug Message";
  InfoMsg << "An Info";
  WarningMsg << "A Warning";
  ErrMsg << "Error Message";
  Alert << "A Fatal Error";

  ErrLogger::instance()->setLevel(log4cpp::Priority::ERROR);

  DebugMsg << "2nd Debug Message";
  InfoMsg << "2nd Info";
  WarningMsg << "2nd Warning";
  ErrMsg << "2nd Error Message";
  Alert << "2nd Fatal Error";

  return 0;
}
