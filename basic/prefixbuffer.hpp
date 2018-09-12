
#ifndef INCLUDEGUARD_PREFIXBUFFER_HPP
#define INCLUDEGUARD_PREFIXBUFFER_HPP

#include <iostream>
#include <streambuf>
#include <string>
#include <locale>
#include <cstdio>
#include <cassert>



class prefixbuffer : public std::streambuf
{

    private: 
        
        bool newline;
        std::string lineprefix;
        std::streambuf * pbuf;
        
    public:
        
        prefixbuffer( std::streambuf* pbuf, const std::string& prefix = "=>" )
        : pbuf( pbuf ), newline( true ), lineprefix( prefix )
        {
            assert( pbuf != nullptr );
            std::string welcome = "Welcome!\n";
            pbuf->sputn( welcome.c_str(), welcome.length() );
        }

        std::streambuf* getbuffer()
        {
            return pbuf;
        }
        
        const std::string& getprefix()
        {
            return lineprefix;
        }
        
        const std::string& setprefix( const std::string& str )
        {
            return lineprefix = str;
        }
    
        virtual int sync()
        {
            return pbuf->pubsync();
        }
        
    protected:
        
        virtual int_type overflow (int_type c) {
        
            if ( c != std::char_traits<char>::eof() ) {
                
                if( newline ) {
                
                    newline = false;
                    
                    if ( pbuf->sputn( lineprefix.c_str(), lineprefix.length() ) == std::char_traits<char>::eof() )
                        return std::char_traits<char>::eof();
                    
                }
                
                if( c == '\n' )
                    newline = true;
            
                c = std::toupper(static_cast<char>(c),getloc());

                if ( pbuf->sputc(c) == std::char_traits<char>::eof() )
                    return std::char_traits<char>::eof();
                
                
            }
            
            return c;
        }
        
    
    
};



inline void contextOpen( const std::string& prefix = "=>" ) 
{
    prefixbuffer* pb_log = new prefixbuffer( std::clog.rdbuf(), prefix );
    prefixbuffer* pb_err = new prefixbuffer( std::cerr.rdbuf(), prefix );
    std::clog.rdbuf( pb_log );
    std::cerr.rdbuf( pb_err );
}

inline void contextClose()
{
    prefixbuffer* pb_log = dynamic_cast<prefixbuffer*>(std::clog.rdbuf());
    prefixbuffer* pb_err = dynamic_cast<prefixbuffer*>(std::cerr.rdbuf());
    std::streambuf* old_log = pb_log->getbuffer();
    std::streambuf* old_err = pb_err->getbuffer();
    std::clog.rdbuf( old_log );
    std::cerr.rdbuf( old_err );
    delete pb_log;
    delete pb_err;
}



// 
// Example usage:
//     prefixbuffer* pb_log = new prefixbuffer( std::clog.rdbuf() );
//     std::clog.rdbuf( pb_log );
//     ...
//     clog << "log entry" << std::endl;
//     ...
//     prefixbuffer* pb_log = clog.rdbuf();
//     std::clog.rdbuf( pb_log.getbuffer() );
//     delete pb_log;
// 
// The preamble and postamble can be paired up as function calls,
// and they can be stacked up to impose several levels.
// However, they must evolve as pairs.






// class prefixbuf : public std::streambuf
// {
//     std::string     prefix;
//     std::streambuf* sbuf;
//     bool            need_prefix;
// 
//     int sync()
//     {
//         return this->sbuf->pubsync();
//     }
//     
//     int overflow(int c)
//     {
//         
//         if (c != std::char_traits<char>::eof())
//         {
//             
//             if ( this->need_prefix
//                  && !this->prefix.empty()
//                  && this->prefix.size() != this->sbuf->sputn( &this->prefix[0], this->prefix.size() )
//                ) return std::char_traits<char>::eof();
//             
//             this->need_prefix = c == '\n';
//         }
//         
//         return this->sbuf->sputc(c);
//         
//     }
//     
//     public:
//         
//         prefixbuf( std::streambuf* sbuf )
//             : prefix(""), sbuf(sbuf), need_prefix(true) {
//         }
// };



#endif
