#ifndef INCLUDEGUARD_MATCSR
#define INCLUDEGUARD_MATCSR

#include <vector>
#include <algorithm>
#include <cmath>

#include "../basic.hpp"
#include "sparsematrix.hpp"






class MatrixCSR:
public LinearOperator /* every matrix is a linear operator */
{

    public:

        explicit MatrixCSR( int rows, int columns, 
                            std::vector<int> A, std::vector<int> C, std::vector<Float> V );
    
        virtual ~MatrixCSR( );

//         explicit MatrixCSR( const ScalingOperator& matrix );
//         explicit MatrixCSR( const DiagonalOperator& matrix );
//         explicit MatrixCSR( const DenseMatrix& );
        
        virtual void check() const override;
        virtual void print( std::ostream& ) const override;
        virtual void printplain( std::ostream& ) const;

        int getnumberofentries() const;
        void sortentries() const;
        void sortandcompressentries() const;

        virtual FloatVector apply( const FloatVector& add, Float scaling ) const override;

    private:

        mutable std::vector<int> A;
        mutable std::vector<int> C;   // column index of each term 
        mutable std::vector<Float> V; // numerical value of each term
    
};



// MatrixCSR operator&( const MatrixCSR& left, const MatrixCSR& right );




#endif
