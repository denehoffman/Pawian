#pragma once

#include "Logger.hh"
#include <cstring>
#include <iostream>

class ErrLogger {
public:
  static ErrLogger& instance();
  ~ErrLogger();

  void setLevel(logging::log_level newLevel);
  void setThreshold(logging::log_level newLevel);
  logging::log_level threshold();
  logging::log_level level();

  static std::ostringstream buffer;

  template<class T>
  inline ErrLogger& operator<< (T val) {
    buffer << val;
    return *this;
  }

  // function that takes a custom stream, and returns it
  typedef ErrLogger& (*ErrLoggerManipulator)(ErrLogger&);

  // take in a function with the custom signature
  ErrLogger& operator<<(ErrLoggerManipulator manip) {
    // call the function, and return it's value
    return manip(*this);
  }

  // define the custom endl for this stream.
  // note how it matches the `ErrLoggerManipulator`
  // function signature
  //  static ErrLogger& endl(ErrLogger& stream) {
  static ErrLogger& endl(ErrLogger& stream) {
    // print a log line
    if (ErrLogger::instance().level() >= ErrLogger::instance().threshold())
      logging::get_logger().log(buffer.str(), ErrLogger::instance().level());
    buffer.str("");
    buffer.clear();
    return stream;
  }

  // this is the type of std::cout
  typedef std::basic_ostream<char, std::char_traits<char> > CoutType;

  // this is the function signature of std::endl
  typedef CoutType& (*StandardEndLine)(CoutType&);

  // define an operator<< to take in std::endl
  ErrLogger& operator<<(StandardEndLine manip) {
    // call the function, but we cannot return it's value
    // manip(std::cout);
    return *this;
  }

private:
  static ErrLogger* theErrLogger;
  logging::log_level logLevel;
  logging::log_level logThreshold;

  ErrLogger();

};



#define endmsg ErrLogger::endl

#define InfoMsg     ErrLogger::instance().setLevel(logging::log_level::INFO);ErrLogger::instance()
#define DebugMsg    ErrLogger::instance().setLevel(logging::log_level::DEBUG);ErrLogger::instance()<<__FILE__<<" line "<<__LINE__<<": "
#define WarningMsg  ErrLogger::instance().setLevel(logging::log_level::WARN);ErrLogger::instance()<<__FILE__<<" line "<<__LINE__<<": "
#define ErrMsg      ErrLogger::instance().setLevel(logging::log_level::ERROR);ErrLogger::instance()<<__FILE__<<" line "<<__LINE__<<": "
// #define Alert       ErrLogger::instance()<<strrchr(__FILE__ ,(int)'/')<< ":"<<__LINE__<<": "
#define Alert       ErrLogger::instance().setLevel(logging::log_level::ERROR);ErrLogger::instance()<<__FILE__<<" line "<<__LINE__<<": "

