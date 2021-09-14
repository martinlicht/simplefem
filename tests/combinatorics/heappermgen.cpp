

/**/

#include <iostream>
#include "../../basic.hpp"
#include "../../combinatorics/indexrange.hpp"
#include "../../combinatorics/indexmap.hpp"
#include "../../combinatorics/heappermgen.hpp"


using namespace std;

extern const char* TestName;
#define TESTNAME( cstr ) const char* TestName = cstr

TESTNAME( "producing permutations via Heap's algorithm" );

int main()
{
    LOG << "Unit Test: " << TestName << endl;
    
    
    if(true)
    {
        
        const int N = 5;
        
        int seed = 77;
        std::vector<int> memo(N);
        std::vector<int> perm(N);
        for( int i = 0; i < N; i++ ) perm[i] = i;
        
        std::vector< std::vector<int> > perms;
        
        HeapsAlgorithmInit( seed, memo, perm );
        
        do {
        
            perms.push_back( perm );
            
        } while ( HeapsAlgorithmStep( seed, memo, perm ) );
        
        
        
        assert( perms.size() == factorial_integer(N) );
        
        for( const auto& some_perm : perms )
            assert( IndexMap( IndexRange(0,N-1), some_perm ).isbijective() );
        
        for( int i = 0; i < perms.size(); i++ )
        for( int j = 0; j < perms.size(); j++ )
            if( i != j )
                assert( perms[i] != perms[j] );
        
        for( const auto& some_perm : perms ) {
            for ( int entry : some_perm ) LOG << entry << space;
            LOG << endl;
        }
        
        LOG << "--------------------" << endl;
        
    }
    
    
    if(true)
    {
    
        std::vector<int> Ns = {0,1,2,3,4,5};
        
        for( const int& N : Ns )
        {
            
            LOG << N << nl;
            
            int seed = 77;
            std::vector<int> memo(N);
            std::vector<int> perm(N);
            for( int i = 0; i < N; i++ ) perm[i] = i;
            
            std::vector< std::vector<int> > perms;
            
            HeapsAlgorithmInit( seed, memo, perm );
            
            do {
            
                perms.push_back( perm );
                
            } while ( HeapsAlgorithmStep( seed, memo, perm ) );
            
            
            
            assert( perms.size() == factorial_integer(N) );
            
            for( const auto& some_perm : perms )
                assert( IndexMap( IndexRange(0,N-1), some_perm ).isbijective() );
            
            for( int i = 0; i < perms.size(); i++ )
            for( int j = 0; j < perms.size(); j++ )
                if( i != j )
                    assert( perms[i] != perms[j] );
            
            
            
            for( const auto& some_perm : perms ) {
                for ( int entry : some_perm ) LOG << entry << space;
                LOG << endl;
            }
            
            LOG << "--------------------" << endl;
            
            
        }
    
    }

    LOG << "Finished Unit Test: " << TestName << endl;
    
    return 0;
}
