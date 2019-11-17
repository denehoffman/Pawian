#include "ErrLogger.hh"

ErrLogger* ErrLogger::theErrLogger = NULL;
std::ostringstream ErrLogger::buffer;

ErrLogger& ErrLogger::instance() {
  if (0 == theErrLogger)
    theErrLogger = new ErrLogger();
  return *theErrLogger;
}

ErrLogger::ErrLogger() {
  setLevel(logging::log_level::TRACE);
  setThreshold(logging::log_level::ERROR);
}

ErrLogger::~ErrLogger() {
}

void ErrLogger::setLevel(logging::log_level newLevel) {
  logLevel = newLevel;
}

void ErrLogger::setThreshold(logging::log_level newLevel) {
  logThreshold = newLevel;
}

logging::log_level ErrLogger::threshold() {
  return logThreshold;
}

logging::log_level ErrLogger::level() {
  return logLevel;
}

