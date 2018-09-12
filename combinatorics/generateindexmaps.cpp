
#include "generateindexmaps.hpp"

#include <algorithm>
#include <vector>
#include <iterator>

#include "../basic.hpp"
#include "indexrange.hpp"
#include "indexmap.hpp"


std::vector<IndexMap> 
generateEmptyMap( const IndexRange& from, const IndexRange& to )
{
    
    from.check();
    to.check();
    assert( from.isempty() );
    
    IndexMap myempty = IndexMap( from, to );
    myempty.check();
    
    std::vector<IndexMap> ret;
    ret.push_back( myempty );
    return ret;
    
}


std::vector<IndexMap> 
generateIndexMaps( const IndexRange& from, const IndexRange& to )
{
        
        from.check();
        to.check();
        
        if( from.isempty() )
            return generateEmptyMap( from, to );
        
        if( to.isempty() )
            return std::vector<IndexMap>();
        
        assert( from.cardinality() > 0 );
        assert( to.cardinality() > 0 );
        
        const int num = integerpower( to.cardinality(), from.cardinality() );
        std::vector<IndexMap> ret;
        ret.reserve( num );
        
        for( int it = 0; it < num; it++ ) {
          
          std::vector<int> values( from.cardinality() );
          int it_clone = it;
          for( int digit = 0; digit < from.cardinality(); digit++ ) {
              values[digit] = to.min() + it_clone % to.cardinality();
              it_clone /= to.cardinality();
          }

          ret.push_back( IndexMap( from, to, values ) );
          
        }
        
        for( const auto& foo : ret )
            foo.check();
        
        return ret;
        
}

std::vector<IndexMap> 
generatePermutations( const IndexRange& ir )
{
        
        ir.check();
        std::vector<IndexMap> allmappings = generateIndexMaps( ir, ir );
        
        std::vector<IndexMap> ret;
        ret.reserve( factorial( ir.cardinality() ) );
        
        for( const auto& mapping : allmappings )
          if( mapping.isbijective() )
            ret.push_back( mapping );
        
//         std::vector<IndexMap> ret( factorial( ir.cardinality() ), IndexMap(ir) );
//         
//         std::vector<IndexMap>
//         std::copy_if( allmappings.begin(), allmappings.end(), ret.begin(),
//                 []( const IndexMap& im ) -> bool { return im.isbijective(); }
//                 );
        
        for( const auto& perm : ret )
            assert( (perm.check(),true) && perm.isbijective() );
        
        return ret;
        
}

int signPermutation( const IndexMap& im )
{
        
    im.check();
    assert( im.isbijective() );
    assert( im.getSourceRange() == im.getDestRange() );
    
    const IndexRange& ir = im.getSourceRange();
    int zaehler = 1;
    int nenner = 1;
    
    /* TODO: checks for integer limits */
    
    for( int s = ir.min(); s <= ir.max(); s++ )
    for( int t = s+1; t <= ir.max(); t++ )
    {
        nenner *= ( t - s );
        zaehler *= ( im[ t ] - im[ s ] );
    }
    
    assert( zaehler % nenner == 0 );
    int ret = zaehler / nenner;
    assert( ret == 1 || ret == -1 );
    return ret;
    
}

std::vector<IndexMap> 
generateSigmas( const IndexRange& from, const IndexRange& to )
{
        from.check();
        to.check();
        
        std::vector<IndexMap> allmappings = generateIndexMaps( from, to );
        
        std::vector<IndexMap> ret;
        assert( to.cardinality() >= from.cardinality() ); // TODO: Handle this warning
        ret.reserve( binomial<int>( to.cardinality(), from.cardinality() ) );
        
        for( const auto& mapping : allmappings )
          if( mapping.isstrictlyascending() )
            ret.push_back( mapping );
        
//         std::vector<IndexMap> ret( , IndexMap(from,to) );
//         
//         std::copy_if( allmappings.begin(), allmappings.end(), ret.begin(),
//                 []( const IndexMap& im ) -> bool { return im.isstrictlyascending(); }
//                 );
        
        for( const auto& sigma : ret )
            assert( (sigma.check(),true) && sigma.isstrictlyascending() );
        
        return ret;
}



