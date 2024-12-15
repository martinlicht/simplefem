#include <iostream>
#include <vector>

using namespace std;





// Internal function to merge lists using backtracking
void mergeLists(
    const std::vector<int>& list1,
    const std::vector<int>& list2,
    std::vector<int>& current,
    int index1,
    int index2,
    std::vector<std::vector<int>>& results) {

    // If we've used all elements from both lists, store the result
    if (index1 == list1.size() && index2 == list2.size()) {
        results.push_back(current);
        return;
    }

    // If there are remaining elements in list1, add the next one
    if (index1 < list1.size()) {
        current.push_back(list1[index1]);
        mergeLists(list1, list2, current, index1 + 1, index2, results);
        current.pop_back(); // Backtrack
    }

    // If there are remaining elements in list2, add the next one
    if (index2 < list2.size()) {
        current.push_back(list2[index2]);
        mergeLists(list1, list2, current, index1, index2 + 1, results);
        current.pop_back(); // Backtrack
    }
}





std::vector<std::vector<int>> mergeLists( 
    const std::vector<int>& list1,
    const std::vector<int>& list2
){
    std::vector<std::vector<int>> results;
    std::vector<int> current;
    mergeLists( list1, list2, current, 0, 0, results );
    return results;
}



std::vector<std::vector<int>> mergeLists( 
    const std::vector<std::vector<int>>& listlists1,
    const std::vector<std::vector<int>>& listlists2
){
    std::vector<std::vector<int>> results;

    for( const auto& list1 : listlists1 )
    for( const auto& list2 : listlists2 )
    {
        const auto temp = mergeLists( list1, list2 );
        results.insert( results.end(), temp.begin(), temp.end() );
    }

    return results;
}


std::vector<std::vector<int>> mergeLists( 
    const std::vector<std::vector<std::vector<int>>>& listlistlists
){
    if( listlistlists.size() == 0 ) return {};

    std::vector<std::vector<int>> result = listlistlists[0];
    
    for( int t = 1; t < listlistlists.size(); t++ )
    {
        result = mergeLists( result, listlistlists[t] );
    }

    return result;
}




int main() {
    // Example input
    std::vector<int> list1 = {-1, -2 };
    std::vector<int> list2 = { 1,  2 };
    std::vector<int> list3 = {10, 20 };

    std::vector<std::vector<std::vector<int>>> initial = { { list1 }, { list2 }, { list3 } };

    std::vector<std::vector<int>> results = mergeLists( initial );
    
    // Print results
    for( int c = 0; c < results.size(); c++ )
    {
        const auto& combination = results[c];
        cout << c << '\t';
        for( int num : combination ){
            cout << num << " ";
        }
        cout << endl;
    }

    return 0;
}
