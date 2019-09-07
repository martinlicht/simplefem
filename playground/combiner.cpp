
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

template<class ForwardIt, class BinaryPredicate, class BinaryOp>
ForwardIt combine_equivalents(ForwardIt first, ForwardIt last, BinaryPredicate p, BinaryOp op)
{
    if (first == last)
        return last;

    ForwardIt dest = first;
    ForwardIt src  = first;
    while (++src != last) {
        if ( p( *dest, *src ) ) {
            *dest = op( std::move( *dest ), *src );
        } else {
            ++dest;
            *dest = std::move( *src );
	  }
    }
    return ++dest;
}



struct Payment{
    std::string name;
    int amount;
};

int main()
{
    
    std::vector<Payment> ps = {
        { "Alice", 100 },
        { "Bob", 50 },
        { "Charles", 110 },
        { "Dopinder", 70 },
        { "Eve", 40 },
        { "Charles", 60 },
        { "Alice", 30 },
        { "Eve", 90 },
        { "Alice", 110 },
        { "Fangyin", 30 }
    };
    
    for( auto& p : ps ) std::cout << p.name << ":" << p.amount << std::endl;
    std::cout << std::endl;
    
    std::sort( ps.begin(), ps.end(), 
        []( Payment a, Payment b )  -> bool{ return a.name < b.name; }
    );
    
    for( auto& p : ps ) std::cout << p.name << ":" << p.amount << std::endl;
    std::cout << std::endl;
    
    auto last = combine_equivalents( ps.begin(), ps.end(), 
        []( const Payment& a, const Payment& b ) -> bool{ 
            return a.name == b.name;
        },
        []( const Payment& a, const Payment& b ) -> Payment{ 
            std::cout << "Combine " 
                      << a.name << ":" << a.amount << " | "
                      << b.name << ":" << b.amount << std::endl;
            return { a.name, a.amount + b.amount };
        } 
    );
    
    std::cout << std::endl;
    
    for( auto& p : ps ) std::cout << p.name << ":" << p.amount << std::endl;
    std::cout << std::endl;
    
    ps.erase( last, ps.end() );
    
    for( auto& p : ps ) std::cout << p.name << ":" << p.amount << std::endl;
    std::cout << std::endl;
    
    return 0;
}
