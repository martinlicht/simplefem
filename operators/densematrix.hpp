#ifndef INCLUDEGUARD_DENSEMATRIX
#define INCLUDEGUARD_DENSEMATRIX

#include <vector>
#include "../basic.hpp"
#include "../combinatorics/indexmap.hpp"
#include "linearoperator.hpp"


class DenseMatrix:
public LinearOperator /* every matrix is a linear operator */
{

    public:
	
		DenseMatrix(int,int);
		virtual ~DenseMatrix();
		
		virtual void check() const override;
		virtual void print( std::ostream& ) const override;
		
		virtual void applyadd( FloatVector& dest, const FloatVector& add, Float s, Float t ) const override;
        
        Float get(int,int) const;
        void set(int,int,Float);
        Float& operator()( int, int );
        const Float& operator()( int, int ) const;
        
        void zeromatrix();
        void randommatrix();
        void unitmatrix();
        void scale( Float );
        void set( Float );
        void indexmapping( const IndexMap& );
        
        DenseMatrix submatrix( const IndexMap& rows, const IndexMap& columns ) const;
        Float determinant() const;
        
        DenseMatrix adjunctMatrix() const;
        DenseMatrix transpose() const;
        // DenseMatrix inverse() const;
        // void gaussfactorization( DenseMatrix&, DenseMatrix&, DenseMatrix&, DenseMatrix&, Float& ) const;

		FloatVector getrow( int r ) const;
		void setrow( int r, const FloatVector& row );
		FloatVector getcolumn( int c ) const;
		void setcolumn( int c, const FloatVector& column );

        void add( const DenseMatrix& );
        void add( Float, const DenseMatrix& );
        void add( Float, Float, const DenseMatrix& );
        static DenseMatrix MatrixMult( const DenseMatrix&, const DenseMatrix& );
        
        bool issquare() const;
        
    private:
	
        std::vector<Float> entries;
	
};
  
  

inline DenseMatrix& operator+=( DenseMatrix& left, const DenseMatrix& right )
{
    left.add( right );
    return left;
}

inline DenseMatrix& operator-=( DenseMatrix& left, const DenseMatrix& right )
{
    left.add( -1., right );
    return left;
}

inline DenseMatrix& operator*=( DenseMatrix& left, const DenseMatrix& right )
{
    DenseMatrix temp( left );
    temp = DenseMatrix::MatrixMult( left, right );
    left = temp;
    return left;
}

inline DenseMatrix& operator*=( DenseMatrix& left, Float right )
{
    left.scale( right );
    return left;
}

inline DenseMatrix& operator/=( DenseMatrix& left, Float right )
{
    left.scale( 1. / right );
    return left;
}




inline DenseMatrix operator+( const DenseMatrix& left, const DenseMatrix& right )
{
    DenseMatrix ret( left );
    ret += right;
    return ret;
}

inline DenseMatrix operator-( const DenseMatrix& left, const DenseMatrix& right )
{
    DenseMatrix ret( left );
    ret -= right;
    return ret;
}

inline DenseMatrix operator*( const DenseMatrix& left, const DenseMatrix& right )
{
    DenseMatrix ret( left );
    ret *= right;
    return ret;
}

inline DenseMatrix operator*( Float left, const DenseMatrix& right )
{
    DenseMatrix ret( right );
    ret *= left;
    return ret;
}

inline DenseMatrix operator*( const DenseMatrix& left, Float right )
{
    DenseMatrix ret( left );
    ret *= right;
    return ret;
}

inline DenseMatrix operator/( const DenseMatrix& left, Float right )
{
    DenseMatrix ret( left );
    ret /= right;
    return ret;
}



#endif