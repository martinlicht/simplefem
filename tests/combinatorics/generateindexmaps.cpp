

/**/

#include <ostream>
#include "../../basic.hpp"
#include "../../combinatorics/indexrange.hpp"
#include "../../combinatorics/indexmap.hpp"
#include "../../combinatorics/generateindexmaps.hpp"


using namespace std;

int main()
{
        LOG << "Unit Test for Index Map Generators" << endl;
        
        if(true)
        {
                
                LOG << "Teste generator for empty index map" << endl;
                
                IndexRange from( 1,-3 );
                IndexRange targ( 2,5 );
                
                std::vector<IndexMap>  all = generateEmptyMap( from, targ );
                
                for( const IndexMap& im : all )
                        LOG << im << endl;
                
                LOG << "Tested" << endl;
                
        }
        
        if(true)
        {
                
                LOG << "Teste generator for all: [1,3] -> [1,3]" << endl;
                
                IndexRange bereich( 1,3 );
                std::vector<IndexMap> all = generateIndexMaps( bereich, bereich );

                for( const IndexMap& im : all )
                        LOG << im << endl;
                
                LOG << "Tested" << endl;
                
        }
        
        if(true)
        {
                
                LOG << "Teste generator for permutations of [1,3]" << endl;
                
                IndexRange bereich( 1,3 );
                std::vector<IndexMap>  all = generatePermutations( bereich );

                for( const IndexMap& im : all )
                        LOG << im << endl << "signum: " << signPermutation( im ) << endl;
                
                LOG << "Tested" << endl;
                
        }
        
        if(true)
        {
                
                LOG << "Teste generator for sigmas [1,3] -> [2,5]" << endl;
                
                IndexRange from( 1,3 );
                IndexRange targ( 2,5 );
                std::vector<IndexMap>  all = generateSigmas( from, targ );

                for( const IndexMap& im : all )
                        LOG << im << endl;
                
                LOG << "Tested" << endl;
                
        }
        
        
        
        
        if(true)
        {
                
                LOG << "Test generator for empty index map" << endl;
                
                const IndexRange from( 1,-3 );
                const IndexRange targ( 2,5 );
                
                const std::vector<IndexMap>  all = generateEmptyMap( from, targ );
                
                assert( all.size() == 1 );
                
                for( const IndexMap& im : all )
                    LOG << im << endl;
                
                LOG << "Tested" << endl;
                
        }
        
        if(true)
        {
                
                LOG << "Test generator for general index mappings" << endl;
                
                const std::vector<int> Ns = { -1, 0, 1, 2, 3, 4 };
                const std::vector<int> Ks = { -1, 0, 1, 2, 3, 4 };
                
                
                for( const int& N : Ns ) 
                for( const int& K : Ks ) 
                {
                
                    const IndexRange source( 1, N );
                    const IndexRange target( 1, K );
                    
                    const std::vector<IndexMap> all = generateIndexMaps( source, target );
                    
                    if( source == target )
                        assert( all == generateIndexMaps( source) );

                    if( source.cardinality() == 0 )
                        assert( all.size() == 1 );
                    else if( target.cardinality() == 0 )
                        assert( all.size() == 0 );
                    else
                        assert( all.size() == power_integer( target.cardinality(), source.cardinality() ) );
                    
                    for( int i = 0; i < all.size(); i++ )
                    for( int j = 0; j < all.size(); j++ )
                        if( i != j )
                            assert( all[i] != all[j] );
                    
                    for( const IndexMap& im : all )
                        LOG << im << endl;
                    
                }
                
                LOG << "Tested" << endl;
                
        }
        
        if(true)
        {
                
                LOG << "Test generator for permutations of [1,3]" << endl;
                
                const std::vector<int> Ns = { -1, 0, 1, 2, 3, 5 };
                
                
                for( const int& N : Ns ) {
                    
                    const IndexRange bereich( 1, N );
                    const std::vector<IndexMap>  all = generatePermutations( bereich );

                    if( bereich.cardinality() > 0 )
                        assert( all.size() > 0 );
                    
                    assert( all.size() == factorial_integer( bereich.cardinality() ) );
                    
                    for( int i = 0; i < all.size(); i++ )
                    for( int j = 0; j < all.size(); j++ )
                        if( i != j )
                            assert( all[i] != all[j] );
                    
                    for( const IndexMap& im : all )
                            LOG << im << endl << "signum: " << signPermutation( im ) << endl;
                    
                    int number_of_even_permutations = 0;
                    
                    for( const IndexMap& im : all ) {
                        assert( im.getSourceRange() == im.getTargetRange() );
                        for( const int& i : bereich ) {
                            assert( bereich.min() <= im[i] and im[i] <= bereich.max() );
                            for( const int& j : bereich  )
                                if( i != j )
                                    assert( im[i] != im[j] );
                        }
                        assert( im.isbijective() );
                        
                        if( signPermutation( im ) == 1 ) 
                            number_of_even_permutations++;
                    }
                    
                    if( bereich.cardinality() > 1 )
                        assert( number_of_even_permutations * 2 == all.size() );
                    else
                        assert( number_of_even_permutations     == all.size() );
                    
                    
                    
                }
                
                LOG << "Tested" << endl;
                
        }
        
        if(true)
        {
                
                LOG << "Test generator for sigmas [1,3] -> [2,5]" << endl;
                
                const std::vector<int> Ns = { -1, 0, 1, 2, 3, 5 };
                const std::vector<int> Ls = {  1, 2, 3, 4, 5, 7, 9 };
                
                
                for( const int& N : Ns ) 
                for( const int& L : Ls ) 
                {
                    
                    IndexRange source( 1, N );
                    IndexRange target( 2, L );
                    
                    if( source.cardinality() > target.cardinality() )
                        continue;
                    
                    std::vector<IndexMap>  all = generateSigmas( source, target );
                    
                    assert( all.size() == binomial_integer( target.cardinality(), source.cardinality() ) );
                    
                    for( const IndexMap& sigma : all ) {
                        assert( sigma.getSourceRange() == source );
                        assert( sigma.getTargetRange() == target );
                        for( const int& i : source ) {
                            assert( target.min() <= sigma[i] and sigma[i] <= target.max() );
                            if( i != source.max() )
                                assert( sigma[i] < sigma[i+1] );
                        }
                        assert( sigma.isstrictlyascending() );
                        
                    }
                    
                    for( const IndexMap& sigma : all )
                        LOG << sigma << endl;
                    
                }
                
                LOG << "Tested" << endl;
                
        }
        
        LOG << "Finished Unit Test" << endl;
        
        return 0;
}
