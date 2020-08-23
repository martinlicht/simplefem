
#include "herzogsoodhalter.hpp"

#include "../operators/floatvector.hpp"



HerzogSoodhalterMethod::HerzogSoodhalterMethod( const LinearOperator& op )
: IterativeSolver(), A( op )
{
    this->max_iteration_count = op.getdimin();
    HerzogSoodhalterMethod::check();
}

HerzogSoodhalterMethod::~HerzogSoodhalterMethod()
{
    HerzogSoodhalterMethod::check();
}

void HerzogSoodhalterMethod::check() const
{
    IterativeSolver::check();
    assert( A.getdimin() == A.getdimout() );
}

void HerzogSoodhalterMethod::print( std::ostream& os ) const
{
    os << "Print Herzog-Soodhalter Method." << std::endl;
}



  
void HerzogSoodhalterMethod::solve( FloatVector& x, const FloatVector& b ) const
{
    check();
    x.check();
    b.check();
    
    assert( x.getdimension() == b.getdimension() );
    assert( A.getdimin()  == x.getdimension() );
    assert( A.getdimout() == b.getdimension() );
    
    const int dimension = A.getdimin();

    FloatVector v0( dimension, 0. );
    FloatVector v1( dimension, 0. );
    FloatVector w0( dimension, 0. );
    FloatVector w1( dimension, 0. );
    
    FloatVector vn( dimension, 0. );
    FloatVector wn( dimension, 0. );
    FloatVector  p( dimension, 0. );
    
    Float gamma = 0.;
    Float eta;
    Float s0, s1, c0, c1;
    
    recent_iteration_count = 0;

    while( recent_iteration_count < max_iteration_count ){
        
        
        bool restart_condition = (recent_iteration_count == 0);
        
        bool residual_seems_small = (eta < tolerance);
        
        if( restart_condition or residual_seems_small ) {
            
            v0.zero();
            w0.zero();
            v1 = b - A * x;
            w1.zero();
            
            gamma = v1.norm();
            v1 /= gamma;
            
            assert( gamma > 0. );
            
            s0 = s1 = 0;
            c0 = c1 = 1;
            
            eta = gamma;
            
        }
        
        bool residual_is_small = (eta < tolerance);
        
        if( residual_is_small )
            break;

            Float temp;
            temp = (b-A*x).norm();
//             temp = gamma;
            LOG << recent_iteration_count << space << temp << space << eta << space << temp/eta;
            
        {
            

            FloatVector p = A * v1;
 
            Float delta = p * v1;
 
            FloatVector vn = p - delta * v1 - gamma * v0;
 
            Float gamma_n = vn.norm();
            vn /= gamma_n;
            assert( gamma_n > 0. );
 
            Float alpha_0 = c1 * delta - c0 * s1 * gamma;
            Float alpha_1 = std::sqrt( alpha_0 * alpha_0 + gamma_n * gamma_n );
            Float alpha_2 = s1 * delta + c0 * c1 * gamma;
            Float alpha_3 = s0 * gamma;
 
            assert( alpha_1 > 0. );

            Float cn = alpha_0 / alpha_1;
            Float sn = gamma_n / alpha_1;
            
            FloatVector wn = ( v1 - alpha_2 * w1 - alpha_3 * w0 ) / alpha_1;
            x = x + cn * eta * wn;
 
            eta = - sn * eta;
            
            LOG << "\t" << alpha_0 << space << alpha_1 << space << alpha_2 << space << alpha_3;
            LOG << "\t" << cn << space << sn << space << eta;
            LOG << "\t" << p.norm() << space << wn.norm();
            
            
            v0 = v1;
            w0 = w1;
            v1 = vn;
            w1 = wn;
            
            gamma = gamma_n;
            
            c0 = c1;
            s0 = s1;
            c1 = cn;
            s1 = sn;

        }

        
        
        
        
        
        
        recent_iteration_count++;
        
    };
    
    /* HOW DID WE FINISH ? */
    
    Float resnorm = ( b - A * x ).norm();
    recent_deviation = resnorm;
    
    if( recent_deviation > tolerance ) {
        LOG << "MINRES process has failed. (" << recent_iteration_count << "/" << max_iteration_count << ") : " << recent_deviation << "/" << tolerance;
    } else { 
        LOG << "MINRES process has succeeded. (" << recent_iteration_count << "/" << max_iteration_count << ") : " << recent_deviation << "/" << tolerance;
    }

    
}
  

  


/*
 * 
 * v0 = 0
 * w0 = 0
 * v1 = b - A x
 * w1 = 0
 * 
 * gamma = norm(v1)
 * v1 = v1/gamma
 * 
 * s0 = 0
 * s1 = 0
 * c0 = 1
 * c1 = 1
 * 
 * eta = gamma
 * 
 * FOR j = 1 ..... 
 * 
 *  Av = A * v
 *  
 *  delta = Av * v
 *  
 *  v_new = Av - delta * v1 - gamma_1 v0
 *  
 *  gamma_new = norm( v_new )
 *  v_new /= gamma
 *  
 *  alpha_0 = c1 delta - c0 s1 gamma
 *  alpha_1 = sqrt( alpha_0 * alpha_0 + gamma_new * gamma_new )
 *  alpha_2 = s1 delta - c0 c1 gamma 
 *  alpha_3 = s0 gamma 
 * 
 *  c_new = alpha0 / alpha1
 *  s_new = gamma_new / alpha1
 *  
 *  w_new = ( v1 - alpha2 w1 - alpha3 w0 ) / alpha1
 *  x = x - c_new eta w_new
 *  
 *  eta = - s_new eta 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
 * 
*/
