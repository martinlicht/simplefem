

/**/

#include <ostream>
#include "../../basic.hpp"
#include "../../combinatorics/generatemultiindices.hpp"


using namespace std;

int main()
{
    LOG << "Unit Test for Multiindex Generators" << nl;
    
    if(true)
    {
        
        LOG << "First bulk" << nl;
    
        IndexRange ir( 1, 4 );
        std::vector<MultiIndex> vmi = generateMultiIndices( ir, 2);
        
        for( const auto& mi : vmi )
            LOG << mi;
        LOG << nl;
        
    }
    
    if(true)
    {
        
        LOG << "Second bulk" << nl;
    
        IndexRange ir( 0, 2 );
        std::vector<MultiIndex> vmi = generateMultiIndices( ir, 2);
        
        for( const auto& mi : vmi )
            LOG << mi;
        LOG << nl;
    
    }
    
    if(true)
    {
        
        LOG << "Third bulk" << nl;
    
        IndexRange ir( 0, 2 );
        std::vector<MultiIndex> vmi = generateMultiIndices( ir, 3);

        for( const auto& mi : vmi )
            LOG << mi;
        LOG << nl;
        
    }
    
    if(true)
    {
      
        LOG << "Zero degree" << nl;
    
        IndexRange ir( 0, 2 );
        std::vector<MultiIndex> vmi = generateMultiIndices( ir, 0);

        for( const auto& mi : vmi )
            LOG << mi;
        LOG << nl;
      
    }
    
    if(true)
    {
      
        LOG << "First degree" << nl;
    
        IndexRange ir( 0, 2 );
        std::vector<MultiIndex> vmi = generateMultiIndices( ir, 1);

        for( const auto& mi : vmi )
            LOG << mi;
        LOG << nl;
      
    }
    
    if(true)
    {
      
        LOG << "Degree 7 in one variable" << nl;
    
        IndexRange ir( 2, 2 );
        std::vector<MultiIndex> vmi = generateMultiIndices( ir, 7);

        for( const auto& mi : vmi )
            LOG << mi;
        LOG << nl;
      
    }
    
    if(true)
    {

        LOG << "MultiIndex over empty range " << nl;
    
        IndexRange ir( 2, -3 );
        std::vector<MultiIndex> vmi = generateMultiIndices( ir, 3);

        for( const auto& mi : vmi )
            LOG << mi;
        LOG << nl;
    
    }
    
    
    
    if(true)
    {
            
        LOG << "Test generator for general multiindices" << nl;
        
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
                LOG << mi << nl;
            
            
        }
        
        LOG << "Tested" << nl;
            
    }
    
    LOG << "Finished Unit Test" << nl;
    
    return 0;
}
