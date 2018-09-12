
#include <iostream>
#include <string>
#include <functional>
#include <cassert>

#include "logger.hpp"






Logger::Logger( std::ostream& os, const std::string& prefix, const std::string& postfix )
: Logger( os, Logger::affix_write( prefix ), Logger::affix_write( postfix ) )
{
    
}


Logger::Logger( Logger& parent, const std::string& prefix, const std::string& postfix )
: Logger( parent, Logger::affix_write( prefix ), Logger::affix_write( postfix ) )
{
    
}



Logger::Logger( std::ostream& os, const std::function<void(Logger&)>& prefix, const std::function<void(Logger&)>& postfix )
: internalstream( os ), prefix( prefix ), postfix( postfix )
{
    
}


Logger::Logger( Logger& parent, const std::function<void(Logger&)>& prefix, const std::function<void(Logger&)>& postfix )
: internalstream( parent.getstream() ), prefix( prefix ), postfix( postfix )
{
    prefix( *this );
}



Logger::~Logger()
{
    postfix( *this );
    internalstream.flush();
}



std::ostream& Logger::getstream()
{
    return internalstream;
}

Logger& Logger::write( const std::string& str )
{
    internalstream << str;
}
    
    





