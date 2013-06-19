#include "ErrLogger.hh"

ErrLogger* ErrLogger::theErrLogger = NULL;

ErrLogger* ErrLogger::instance()
{
  if (0==theErrLogger)
    theErrLogger = new ErrLogger();
  return theErrLogger;
}

ErrLogger::ErrLogger()
{
  app = new log4cpp::OstreamAppender("FileAppender", &std::cout);
  layout = new log4cpp::PatternLayout();
  layout->setConversionPattern(std::string("%p: %m%n"));
  app->setLayout(layout);
  category = std::shared_ptr<log4cpp::Category>(&(log4cpp::Category::getInstance("ErrLogger")));
  category->setAdditivity(false);
  category->setAppender(app);
  category->setPriority(log4cpp::Priority::INFO);
}

ErrLogger::~ErrLogger()
{
//   delete category;
//   delete app;
//   delete layout;
}

log4cpp::Category& ErrLogger::logger()
{
  return *category.get();
}

void ErrLogger::setLevel(log4cpp::Priority::PriorityLevel level)
{
  category->setPriority(level);
}
