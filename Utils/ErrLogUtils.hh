#pragma once

#include "ConfigParser/ParserBase.hh"
#include "ErrLogger/ErrLogger.hh"

void setErrLogMode( const ParserBase::enErrLogMode& erlMode ) {
	switch(erlMode) {
	case ParserBase::debug :
		ErrLogger::instance().setThreshold(logging::log_level::DEBUG);
		break;
	case ParserBase::trace :
		ErrLogger::instance().setThreshold(logging::log_level::INFO);
		break;
	case ParserBase::notice :
		ErrLogger::instance().setThreshold(logging::log_level::NOTICE);
		break;
	case ParserBase::warning :
		ErrLogger::instance().setThreshold(logging::log_level::WARN);
		break;
	case ParserBase::error :
		ErrLogger::instance().setThreshold(logging::log_level::ERROR);
		break;
	case ParserBase::alert :
		ErrLogger::instance().setThreshold(logging::log_level::ERROR);
		break;
	default:
		ErrLogger::instance().setThreshold(logging::log_level::DEBUG);
	}
}

