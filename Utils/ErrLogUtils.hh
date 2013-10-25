#pragma once

#include "ConfigParser/ParserBase.hh"
#include "ErrLogger/ErrLogger.hh"

void setErrLogMode( const ParserBase::enErrLogMode& erlMode ) {
	switch(erlMode) {
	case ParserBase::debug :
		ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);
		break;
	case ParserBase::trace :
		ErrLogger::instance()->setLevel(log4cpp::Priority::INFO);
		break;
	case ParserBase::routine :
		ErrLogger::instance()->setLevel(log4cpp::Priority::INFO);
		break;
	case ParserBase::warning :
		ErrLogger::instance()->setLevel(log4cpp::Priority::WARN);
		break;
	case ParserBase::error :
		ErrLogger::instance()->setLevel(log4cpp::Priority::ERROR);
		break;
	case ParserBase::alert :
		ErrLogger::instance()->setLevel(log4cpp::Priority::ALERT);
		break;
	default:
		ErrLogger::instance()->setLevel(log4cpp::Priority::DEBUG);
	}
}

