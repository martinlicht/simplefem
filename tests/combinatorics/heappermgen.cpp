

/**/

#include <iostream>
#include "../../basic.hpp"
#include "../../combinatorics/heappermgen.hpp"


using namespace std;

int main()
{
        cout << "Unit Test for producting permutations via Heap's algorithm" << endl;
        
        const int N = 5;
        
        int i = 77;
        std::vector<int> c(N);
        std::vector<int> a(N);
        for( int i = 0; i < N; i++ ) a[i] = i;
        
        for ( int entry : a )
          cout << entry << space;
        cout << endl;
        
        cout << "--------------------" << endl;
        
        HeapsAlgorithmInit( i, c, a );
        
        do {
          
          for ( int entry : a )
            cout << entry << space;
          cout << endl;
          
        } while ( HeapsAlgorithmStep( i, c, a ) );
        
        cout << "Finished Unit Test" << endl;
        
        return 0;
}
