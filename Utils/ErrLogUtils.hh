#pragma once

#include "ConfigParser/ParserBase.hh"
#include "ErrLogger/ErrLogger.hh"

void setErrLogMode( const ParserBase::enErrLogMode& erlMode ) {
	switch(erlMode) {
	case ParserBase::debug :
		ErrLogger::instance().setLevel(logging::log_level::DEBUG);
		break;
	case ParserBase::trace :
		ErrLogger::instance().setLevel(logging::log_level::INFO);
		break;
	case ParserBase::routine :
		ErrLogger::instance().setLevel(logging::log_level::INFO);
		break;
	case ParserBase::warning :
		ErrLogger::instance().setLevel(logging::log_level::WARN);
		break;
	case ParserBase::error :
		ErrLogger::instance().setLevel(logging::log_level::ERROR);
		break;
	case ParserBase::alert :
		ErrLogger::instance().setLevel(logging::log_level::ERROR);
		break;
	default:
		ErrLogger::instance().setLevel(logging::log_level::DEBUG);
	}
}

