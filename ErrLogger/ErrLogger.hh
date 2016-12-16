#pragma once

#include "log4cpp/Category.hh"
#include "log4cpp/OstreamAppender.hh"
#include "log4cpp/PatternLayout.hh"
#include "memory"

#include <string.h>

class ErrLogger
{
public:
  static ErrLogger* instance();
  ~ErrLogger();

  log4cpp::Category& logger();
  void setLevel(log4cpp::Priority::PriorityLevel level);

private:
  static ErrLogger* theErrLogger;

  ErrLogger();

  std::shared_ptr<log4cpp::Category> category;
  log4cpp::Appender* app;
  log4cpp::PatternLayout* layout;
};

#define endmsg log4cpp::eol
#define DebugMsg ErrLogger::instance()->logger()<<log4cpp::Priority::DEBUG<<strrchr(__FILE__ ,(int)'/')<< ":"<<__LINE__<<": "
#define InfoMsg     ErrLogger::instance()->logger()<<log4cpp::Priority::INFO
#define WarningMsg  ErrLogger::instance()->logger()<<log4cpp::Priority::WARN<<strrchr(__FILE__ ,(int)'/')<< ":"<<__LINE__<<": "
#define ErrMsg   ErrLogger::instance()->logger()<<log4cpp::Priority::ERROR<<strrchr(__FILE__ ,(int)'/')<< ":"<<__LINE__<<": "
#define Alert    ErrLogger::instance()->logger()<<log4cpp::Priority::ALERT<<strrchr(__FILE__ ,(int)'/')<< ":"<<__LINE__<<": "
