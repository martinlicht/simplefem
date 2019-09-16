
// #include "heappermgen.hpp"

#include <algorithm>
#include <vector>
#include <iostream>
#include <iterator>

#include "../basic.hpp"



void SigmaAlgorithmInit( int n, int k, std::vector<int>& s )
{
  assert( 0 <= k && k <= n && k == s.size() );
  for( int j = 0; j < s.size(); j++ ) s[j] = j;
}


bool SigmaAlgorithmStep( int n, int k, std::vector<int>& s )
{
  assert( 0 <= k && k <= n && k == s.size() );
  
  int i = 0;
  
  while( i < k-1 )
      if( s[i] < s[i+1]-1 )
          break;
      else
          i++;
  
  if( i == k-1 && s[i] == n ) return false;
  
  assert( ( i==k-1 and s[i] < n ) or s[i] < s[i+1] );
  
  s[i]++;
  
  for( int j = 0; j < i; j++ ) s[j] = j;
  
  return true;
}





void RevSigmaAlgorithmInit( int n, int k, std::vector<int>& s )
{
  assert( 0 <= k && k <= n && k == s.size() );
  for( int j = 0; j < s.size(); j++ ) s[j] = j;
}


bool RevSigmaAlgorithmStep( int n, int k, std::vector<int>& s )
{
  assert( 0 <= k && k <= n && k == s.size() );
  
  int i = k-1;
  
  while( i > 0 )
      if( i==k-1 and s[i]<n )
          break;
      else if( s[i] < s[i+1]-1 )
          break;
      else
          i--;
  
  if( i == 0 && s[i] == n-k+1 ) return false;
  
  assert( ( i==k-1 and s[i] < n ) or s[i] < s[i+1] );
  
  s[i]++;
  
  for( int j = i+1; j < k; j++ ) s[j] = s[j-1]+1;
  
  return true;
}



int main()
{
    int n = 5; 
    int k = 4;
    std::vector<int> s(k);
    
    int iter = 0;
    RevSigmaAlgorithmInit(n,k,s);
    do{
        std::cout << iter++ << ": ";
        for( auto i : s ) std::cout << i << ' ';
        std::cout << '\n';
    }while(RevSigmaAlgorithmStep(n,k,s));
    return 0;
}

