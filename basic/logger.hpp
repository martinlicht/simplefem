#ifndef INCLUDEGUARD_LOGGER_HPP
#define INCLUDEGUARD_LOGGER_HPP

#include <iostream>
#include <string>
#include <functional>
#include <cassert>





class Logger
{
public:
    
    explicit Logger( std::ostream&, const std::string& prefix = "", const std::string& postfix = "" );
    explicit Logger( Logger& parent, const std::string& prefix = "", const std::string& postfix = "" );
    
    explicit Logger( std::ostream&, const std::function<void(Logger&)>& prefix = Logger::affix_do_nothing(), const std::function<void(Logger&)>& postfix = Logger::affix_do_nothing() );
    explicit Logger( Logger& parent, const std::function<void(Logger&)>& prefix = Logger::affix_do_nothing(), const std::function<void(Logger&)>& postfix = Logger::affix_do_nothing() );
    
    ~Logger();
    
    std::ostream& getstream();
    
    Logger& write( const std::string& str );
    
    template<class T>
    Logger& operator<<( const T& t )
    {
        getstream() << t;
        return *this;
    }

    
    Logger& operator<<(std::ostream& (*f)(std::ostream&))
    {
        f( getstream() );
        return *this;
    }
    
    static const std::function<void(Logger&)> affix_do_nothing(){ return [](Logger&){ return; }; };

    static const std::function<void(Logger&)> affix_write( const std::string& str )
    {
        return [str](Logger& logger ){ logger.write( str ); };
    }
    

    
private:
    
    std::ostream& internalstream;
    std::function<void(Logger&)> prefix;
    std::function<void(Logger&)> postfix;

    
};








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
    return *this;
}
    
    









#endif
