#ifndef INCLUDEGUARD_LOGGER_HPP
#define INCLUDEGUARD_LOGGER_HPP

#include <cassert>
#include <functional>
#include <ostream>
#include <string>



class NilLogger
{
public:
    
    template<class T>
    NilLogger& operator<<( const T& )
    {
        return *this;
    }

    
    NilLogger& operator<<(std::ostream& (*)(std::ostream&))
    {
        return *this;
    }
    
};




class Logger
{
    public:
        
        explicit inline Logger( Logger& parent, const std::string& prefix_str = "", const std::string& postfix_str = "" )
        : Logger( parent, Logger::affix_write( prefix_str ), Logger::affix_write( postfix_str ) )
        {
            
        }

        explicit inline Logger( std::ostream& os, const std::string& prefix_str = "", const std::string& postfix_str = "" )
        : Logger( os, Logger::affix_write( prefix_str ), Logger::affix_write( postfix_str ) )
        {
            
        }

        explicit inline Logger( Logger& parent, const std::function<void(Logger&)>& prefix = Logger::affix_do_nothing(), const std::function<void(Logger&)>& postfix = Logger::affix_do_nothing() )
        : internalstream( parent.getstream() ), prefix( prefix ), postfix( postfix )
        {
            
        }

        explicit inline Logger( std::ostream& os, const std::function<void(Logger&)>& prefix = Logger::affix_do_nothing(), const std::function<void(Logger&)>& postfix = Logger::affix_do_nothing() )
        : internalstream( os ), prefix( prefix ), postfix( postfix )
        {
            prefix( *this );
        }




        inline ~Logger()
        {
            postfix( *this );
            internalstream.flush();
        }

        inline std::ostream& getstream()
        {
            return internalstream;
        }
        
        
        
        inline Logger& write( const std::string& str )
        {
            internalstream << str;
            return *this;
        }
        
        
        template<class T>
        Logger& operator<<( const T& t )
        {
            getstream() << t;
            return *this;
        }

        
        Logger& operator<<( std::ostream& (*const f)(std::ostream&) )
        {
            f( getstream() );
            return *this;
        }
        
        
        
        static const std::function<void(Logger&)> affix_do_nothing()
        { 
            return [](Logger&){ return; };
        }

        static const std::function<void(Logger&)> affix_write( const std::string& str )
        {
            return [=](Logger& logger ){ logger.write( str ); };
        }
        

        
    private:
        
        std::ostream& internalstream;
        std::function<void(Logger&)> prefix;
        std::function<void(Logger&)> postfix;

        
};








    
    

#endif
