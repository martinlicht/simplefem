#ifndef INCLUDEGUARD_LOGGER_HPP
#define INCLUDEGUARD_LOGGER_HPP

#include <iostream>
#include <string>
#include <functional>
#include <cassert>





class Logger
{
public:
    
    Logger( std::ostream&, const std::string& prefix = "", const std::string& postfix = "" );
    Logger( Logger& parent, const std::string& prefix = "", const std::string& postfix = "" );
    
    Logger( std::ostream&, const std::function<void(Logger&)>& prefix = Logger::affix_do_nothing(), const std::function<void(Logger&)>& postfix = Logger::affix_do_nothing() );
    Logger( Logger& parent, const std::function<void(Logger&)>& prefix = Logger::affix_do_nothing(), const std::function<void(Logger&)>& postfix = Logger::affix_do_nothing() );
    
    ~Logger();
    
    std::ostream& getstream();
    
    Logger& write( const std::string& str );
    
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



template<typename T>
Logger& operator<<( Logger& logger, const T& t )
{
    logger.getstream() << t;
    return logger;
}




// returns a temporary logger to write stuff to, and line breaks on destruction 
// Example usage:
//     LOG << "This is a short message with a number: " << 5;      
//     ERR << "This is an error message.";      

#define LOG     Logger( clog, "", "\n" )
#define ERR     Logger( cerr, "", "\n" )




// treat the following macros as PRINT 'str' commands
// Example usage:
//     NOTICE "This is a short information"

#define NOTE    Logger( clog, "", "\n" ) <<

#define WARN    Logger( cerr, "", "\n" ) <<
#define ALERT   Logger( cerr, "", "\n" ) <<
#define ERROR   Logger( cerr, "", "\n" ) <<





#endif
