

/**/

#include <ostream>
#include "../../basic.hpp"
#include "../../combinatorics/generatemultiindices.hpp"


using namespace std;

extern const char* TestName;
#define TESTNAME( cstr ) const char* TestName = cstr

TESTNAME( "" );

int main()
{
    LOG << "Unit Test: " << TestName << endl;
    
    if(true)
    {
        
        LOG << "First bulk" << endl;
    
        IndexRange ir( 1, 4 );
        std::vector<MultiIndex> vmi = generateMultiIndices( ir, 2);
        
        for( const auto& mi : vmi )
            LOG << mi;
        LOG << endl;
        
    }
    
    if(true)
    {
        
        LOG << "Second bulk" << endl;
    
        IndexRange ir( 0, 2 );
        std::vector<MultiIndex> vmi = generateMultiIndices( ir, 2);
        
        for( const auto& mi : vmi )
            LOG << mi;
        LOG << endl;
    
    }
    
    if(true)
    {
        
        LOG << "Third bulk" << endl;
    
        IndexRange ir( 0, 2 );
        std::vector<MultiIndex> vmi = generateMultiIndices( ir, 3);

        for( const auto& mi : vmi )
            LOG << mi;
        LOG << endl;
        
    }
    
    if(true)
    {
      
        LOG << "Zero degree" << endl;
    
        IndexRange ir( 0, 2 );
        std::vector<MultiIndex> vmi = generateMultiIndices( ir, 0);

        for( const auto& mi : vmi )
            LOG << mi;
        LOG << endl;
      
    }
    
    if(true)
    {
      
        LOG << "First degree" << endl;
    
        IndexRange ir( 0, 2 );
        std::vector<MultiIndex> vmi = generateMultiIndices( ir, 1);

        for( const auto& mi : vmi )
            LOG << mi;
        LOG << endl;
      
    }
    
    if(true)
    {
      
        LOG << "Degree 7 in one variable" << endl;
    
        IndexRange ir( 2, 2 );
        std::vector<MultiIndex> vmi = generateMultiIndices( ir, 7);

        for( const auto& mi : vmi )
            LOG << mi;
        LOG << endl;
      
    }
    
    if(true)
    {

        LOG << "MultiIndex over empty range " << endl;
    
        IndexRange ir( 2, -3 );
        std::vector<MultiIndex> vmi = generateMultiIndices( ir, 3);

        for( const auto& mi : vmi )
            LOG << mi;
        LOG << endl;
    
    }
    
    
    
    if(true)
    {
            
        LOG << "Test generator for general multiindices" << endl;
        
        const std::vector<int> Ns = { /*-1, 0,*/ 1, 2, 3, 4 };
        const std::vector<int> Rs = { 0, 1, 2, 3, 4, 5, 6 };
        
        
        for( const int& N : Ns ) 
        for( const int& R : Rs ) 
        {
        
            LOG << N << space << R << '$' << nl;
            
            const IndexRange bereich( 1, N );
        
            const std::vector<MultiIndex> all = generateMultiIndices( bereich, R );
            
            assert( all.size() == binomial_integer( bereich.cardinality() + R - 1, R ) );
            
            for( int i = 0; i < all.size(); i++ )
            for( int j = 0; j < all.size(); j++ )
                if( i != j )
                    assert( all[i] != all[j] );
            
            for( const MultiIndex& mi : all ) {
                assert( mi.absolute() == R );
            }
                
            for( const MultiIndex& mi : all )
                LOG << mi << endl;
            
            
        }
        
        LOG << "Tested" << endl;
            
    }
    
    LOG << "Finished Unit Test: " << TestName << endl;
    
    return 0;
}
