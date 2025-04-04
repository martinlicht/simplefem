#ifndef INCLUDEGUARD_LOGGER_HPP
#define INCLUDEGUARD_LOGGER_HPP

#include <cassert>
#include <functional>
#include <ostream>
#include <string>
#include <sstream>



#include <cstdio>

template< typename L, typename... Params >
void printf_into_logger( L logger, const char* formatstring, Params... args )
{
    std::size_t length = std::snprintf(nullptr, 0, formatstring, args... ) + 1;
    char* str = new char[length];
    std::snprintf( str, length, formatstring, args... );
    
    logger << str;

    delete[] str;
}



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




class DirectLogger
{
    public:
        
        explicit inline DirectLogger( DirectLogger& parent, const std::string& prefix_str = "", const std::string& postfix_str = "" )
        : DirectLogger( parent, DirectLogger::affix_write( prefix_str ), DirectLogger::affix_write( postfix_str ) )
        {
            
        }

        explicit inline DirectLogger( std::ostream& os, const std::string& prefix_str = "", const std::string& postfix_str = "" )
        : DirectLogger( os, DirectLogger::affix_write( prefix_str ), DirectLogger::affix_write( postfix_str ) )
        {
            
        }

        explicit inline DirectLogger( DirectLogger& parent, const std::function<void(DirectLogger&)>& prefix = DirectLogger::affix_do_nothing(), const std::function<void(DirectLogger&)>& postfix = DirectLogger::affix_do_nothing() )
        : internalstream( parent.getstream() ), prefix( prefix ), postfix( postfix )
        {
            
        }

        explicit inline DirectLogger( std::ostream& os, const std::function<void(DirectLogger&)>& prefix = DirectLogger::affix_do_nothing(), const std::function<void(DirectLogger&)>& postfix = DirectLogger::affix_do_nothing() )
        : internalstream( os ), prefix( prefix ), postfix( postfix )
        {
            prefix( *this );
        }




        inline ~DirectLogger()
        {
            postfix( *this );
            internalstream.flush();
        }

        inline std::ostream& getstream()
        {
            return internalstream;
        }
        
        
        
        inline DirectLogger& write( const std::string& str )
        {
            internalstream << str;
            return *this;
        }
        
        
        template<class T>
        DirectLogger& operator<<( const T& t )
        {
            getstream() << t;
            return *this;
        }

        
        DirectLogger& operator<<( std::ostream& (*const f)(std::ostream&) )
        {
            f( getstream() );
            return *this;
        }
        
        
        
        static const std::function<void(DirectLogger&)> affix_do_nothing()
        { 
            return []( DirectLogger& ){ return; };
        }

        static const std::function<void(DirectLogger&)> affix_write( const std::string& str )
        {
            return [=]( DirectLogger& logger ){ logger.write( str ); };
        }
        
    private:
        
        std::ostream& internalstream;
        std::function<void(DirectLogger&)> prefix;
        std::function<void(DirectLogger&)> postfix;
        
};








class Logger
{
    public:
        
        explicit inline Logger( 
            Logger& parent, const std::string& prefix_str = "", const std::string& postfix_str = ""
            , const char* filename = "UNKNOWN", const int linenumber = -1 )
        : Logger( parent, Logger::affix_write( prefix_str ), Logger::affix_write( postfix_str ), filename, linenumber )
        {}
        
        explicit inline Logger( 
            std::ostream& os, const std::string& prefix_str = "", const std::string& postfix_str = ""
            , const char* filename = "UNKNOWN", const int linenumber = -1 )
        : Logger( os, Logger::affix_write( prefix_str ), Logger::affix_write( postfix_str ), filename, linenumber )
        {}
        
        explicit inline Logger( 
            Logger& parent, const std::function<void(std::ostream&)>& prefix = Logger::affix_do_nothing(), const std::function<void(std::ostream&)>& postfix = Logger::affix_do_nothing()
            , const char* filename = "UNKNOWN", const int linenumber = -1 )
        : Logger( parent.getstream(), prefix, postfix, filename, linenumber )
        {}
        
        explicit inline Logger( 
            std::ostream& os, const std::function<void(std::ostream&)>& prefix = Logger::affix_do_nothing(), const std::function<void(std::ostream&)>& postfix = Logger::affix_do_nothing()
            , const char* filename = "UNKNOWN", const int linenumber = -1 )
        : internalstream( os ), prefix( prefix ), postfix( postfix ), filename(filename), linenumber(linenumber)
        {
//             internalbuffer.reserve(256);
        }
        



        inline ~Logger()
        {
            
            const auto str = internalbuffer.str();
            
            bool use_prefix_next  = true;
            
            if(true)
            {
                
                for( int c = 0; c < str.size(); c++ )
                {

                    if( use_prefix_next ) { 
                        use_prefix_next = false;
                        prefix( internalstream );
                    }
                    
                    auto character = str.at(c);
                    
                    internalstream << character;

                    if( character == '\n' ) 
                    { 
                        internalstream.flush();
                        use_prefix_next = true;
                    }
                    
                }
                
                postfix( internalstream );

            } else {
                
                prefix( internalstream );
                internalstream << str;
                postfix( internalstream );
                
            }
                    
//             str::size_t curr = 0;
//             prefix( internalstream );
//             str::size_t next = str::find_first_of( '\n', 0 );
//             internalstream << str.substr( 0, next );
//             curr = next;
//             while( curr < str.length() )
//             {
//                 // curr points to a newline character
//                 internalstream << '\n';
//                 prefix( internalstream );
//                 curr++;
//                 if( not curr < str.length() ) break;
//                 str::size_t next = str::find_first_of( '\n', curr );
//                 internalstream << str.substr( curr, next );
//                 next = curr;
//             }
            
            
            if( print_file_and_line ) {
                prefix( internalstream );
                // internalstream << "\e[91m" << filename << ':' << linenumber << "\e[39m" << '\n';
                internalstream << "" << filename << ':' << linenumber << '\n';
                postfix( internalstream );
            }

            #ifndef NDEBUG
            internalstream.flush();
            #endif
            
        }

        inline std::ostream& getstream()
        {
            return internalstream;
        }
        
        
        
        inline Logger& write( const std::string& str )
        {
            internalbuffer << str;
            return *this;
        }
        
        template<class T>
        Logger& operator,( const T& t )
        {
            return ( *this << t );
        }
        
        Logger& operator,( std::ostream& (*const f)(std::ostream&) )
        {
            return ( *this << f );
        }
        
        template<class T>
        Logger& operator<<( const T& t )
        {
            internalbuffer << t;
            return *this;
        }
        
        Logger& operator<<( std::ostream& (*const f)(std::ostream&) )
        {
            f( internalbuffer );
            return *this;
        }
        
        
        
        static const std::function<void(std::ostream&)> affix_do_nothing()
        { 
            return []( std::ostream& ){ return; };
        }

        static const std::function<void(std::ostream&)> affix_write( const std::string& str )
        {
            return [=]( std::ostream& os ){ os << str; };
        }
        
    private:
        
        std::ostringstream internalbuffer;
        std::ostream& internalstream;
        std::function<void(std::ostream&)> prefix;
        std::function<void(std::ostream&)> postfix;
        
        bool print_file_and_line = false;
        const std::string filename;
        const int linenumber;
        
};







    
    

#endif
