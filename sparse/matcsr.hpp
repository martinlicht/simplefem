#ifndef INCLUDEGUARD_MATCSR
#define INCLUDEGUARD_MATCSR

#include <cmath>
#include <algorithm>
#include <memory>
#include <utility>
#include <vector>

#include "../basic.hpp"
#include "sparsematrix.hpp"






class MatrixCSR:
public LinearOperator /* every matrix is a linear operator */
{

    public:

        explicit MatrixCSR( int rows, int columns, 
                            const std::vector<int>& A, 
                            const std::vector<int>& C, 
                            const std::vector<Float>& V );

        explicit MatrixCSR( const SparseMatrix& mat );

        virtual ~MatrixCSR( );

        MatrixCSR( const MatrixCSR& );
        MatrixCSR& operator=( const MatrixCSR& );
        MatrixCSR( MatrixCSR&& );
        MatrixCSR& operator=( MatrixCSR&& );

        virtual std::shared_ptr<LinearOperator> get_shared_pointer_to_clone() const& override {
            std::shared_ptr<MatrixCSR> cloned = std::make_shared<MatrixCSR>( *this );
            return cloned;
        }
        
        virtual std::unique_ptr<LinearOperator> get_unique_pointer_to_heir() && override {
            std::unique_ptr<MatrixCSR> heir = std::make_unique<MatrixCSR>( std::move(*this) );
            return heir;
        }
        
        

//         explicit MatrixCSR( const ScalingOperator& matrix );
//         explicit MatrixCSR( const DiagonalOperator& matrix );
//         explicit MatrixCSR( const DenseMatrix& );
        
        virtual void check() const override;
        virtual void print( std::ostream& ) const override;
        virtual void printplain( std::ostream& ) const;

        int getnumberofentries() const;
        void sortentries() const;
        
        using LinearOperator::apply;
        virtual void apply( FloatVector& dest, const FloatVector& add, Float scaling ) const override;

    private:

        mutable std::vector<int> A;
        mutable std::vector<int> C;   // column index of each term 
        mutable std::vector<Float> V; // numerical value of each term
    
};




#endif
