

/**/

#include <iostream>
#include "../../basic.hpp"
#include "../../combinatorics/generatemultiindices.hpp"


using namespace std;

int main()
{
    cout << "Unit Test for Multiindex Generators" << endl;
    
    cout << "First bulk" << endl;
    
    {
      
      IndexRange ir( 1, 4 );
      std::vector<MultiIndex> vmi = generateMultiIndices( ir, 2);
      
      for( const auto& mi : vmi )
	cout << mi;
      
    }
    
    cout << "Second bulk" << endl;
    
    {
      
      IndexRange ir( 0, 2 );
      std::vector<MultiIndex> vmi = generateMultiIndices( ir, 2);
      
      for( const auto& mi : vmi )
	cout << mi;
      
    }
    
    cout << "Third bulk" << endl;
    
    {
      
      IndexRange ir( 0, 2 );
      std::vector<MultiIndex> vmi = generateMultiIndices( ir, 3);
      
      for( const auto& mi : vmi )
	cout << mi;
      
    }
    
    cout << "Zero degree" << endl;
    
    {
      
      IndexRange ir( 0, 2 );
      std::vector<MultiIndex> vmi = generateMultiIndices( ir, 0);
      
      for( const auto& mi : vmi )
	cout << mi;
      
    }
    
    cout << "First degree" << endl;
    
    {
      
      IndexRange ir( 0, 2 );
      std::vector<MultiIndex> vmi = generateMultiIndices( ir, 1);
      
      for( const auto& mi : vmi )
	cout << mi;
      
    }
    
    cout << "Degree 7 in one variable" << endl;
    
    {
      
      IndexRange ir( 2, 2 );
      std::vector<MultiIndex> vmi = generateMultiIndices( ir, 7);
      
      for( const auto& mi : vmi )
	cout << mi;
      
    }
    
    cout << "MultiIndex over empty range " << endl;
    
    {
      
      IndexRange ir( 2, -3 );
      std::vector<MultiIndex> vmi = generateMultiIndices( ir, 7);
      
      for( const auto& mi : vmi )
	cout << mi;
      
    }
    
    cout << "Finished Unit Test" << endl;
    
    return 0;
}
