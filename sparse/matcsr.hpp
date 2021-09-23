#ifndef INCLUDEGUARD_SPARSE_MATCSR_HPP
#define INCLUDEGUARD_SPARSE_MATCSR_HPP

#include <memory>
#include <utility>
#include <vector>

#include "../basic.hpp"
#include "../operators/floatvector.hpp"
#include "../operators/linearoperator.hpp"
#include "sparsematrix.hpp"




/************************
****
****  Class for Sparse Matrices in CSR format
****  - instantiates LinearOperator
****  
************************/




class MatrixCSR:
public LinearOperator /* every matrix is a linear operator */
{

    public:

        explicit MatrixCSR( int rows, int columns, 
                            const std::vector<int>& A, 
                            const std::vector<int>& C, 
                            const std::vector<Float>& V );

        explicit MatrixCSR( const SparseMatrix& mat );

        explicit MatrixCSR( int rows, int columns );

        virtual ~MatrixCSR( );

        MatrixCSR( const MatrixCSR& );
        MatrixCSR& operator=( const MatrixCSR& );
        MatrixCSR( MatrixCSR&& );
        MatrixCSR& operator=( MatrixCSR&& );

        virtual std::shared_ptr<LinearOperator> get_shared_pointer_to_clone() const& override
        {
            std::shared_ptr<MatrixCSR> cloned = std::make_shared<MatrixCSR>( *this );
            return cloned;
        }
        
        virtual std::unique_ptr<LinearOperator> get_unique_pointer_to_heir() && override
        {
            std::unique_ptr<MatrixCSR> heir = std::make_unique<MatrixCSR>( std::move(*this) );
            return heir;
        }
        
        virtual void check() const override;
        virtual std::string text() const override;
        virtual void print( std::ostream& ) const override;
        virtual void printplain( std::ostream& ) const;

        using LinearOperator::apply;
        virtual void apply( FloatVector& dest, const FloatVector& add, Float scaling ) const override;
        
        
        /* manipulation and information */
        
        void scale ( Float s ); 

        bool isfinite() const;
        
        FloatVector diagonal() const;

        
        /* access and information to internal data */
        
        const int*   getA() const;
        
        const int*   getC() const;
        
        const Float* getV() const;

        int getnumberofentries() const;

        int getnumberofzeroentries() const;

        Float eigenvalueupperbound() const;

//         void sortentries() const;
        
        
        
    private:

        std::vector<int>   A;
        std::vector<int>   C; // column index of each term 
        std::vector<Float> V; // numerical value of each term
    
};





DiagonalOperator InverseDiagonalPreconditioner( const MatrixCSR& mat );





inline MatrixCSR operator*( const MatrixCSR& mat, Float s )
{
    auto foo = mat;
    foo.scale(s);
    return foo;
}

inline MatrixCSR operator*( Float s, const MatrixCSR& mat )
{
    auto foo = mat;
    foo.scale(s);
    return foo;
}



#endif
