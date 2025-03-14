
#include <cmath>

#include "../operators/floatvector.hpp"

#include "systemsolver.hpp"


static const bool cppsys_restart_on_full_dimension = false;
static const bool cppsys_restart_before_finish     = false;



  
int BlockHerzogSoodhalterMethod( 
    FloatVector& x_A, 
    FloatVector& x_C, 
    const FloatVector& b_A, 
    const FloatVector& b_C, 
    const LinearOperator& A, const LinearOperator& Bt, const LinearOperator& B, const LinearOperator& C, 
    Float tolerance,
    int print_modulo,
    const LinearOperator& PAinv, const LinearOperator& PCinv
) {

    x_A.check();
    x_C.check();
    
    b_A.check();
    b_C.check();
    
    A.check(); C.check(); B.check(); Bt.check(); 

    assert( A.getdimin() == A.getdimout() );
    assert( C.getdimin() == C.getdimout() );
    
    const int dimension_A = A.getdimin();
    const int dimension_C = C.getdimin();

    assert( B.getdimin()  == dimension_A );
    assert( B.getdimout() == dimension_C );
    
    assert( Bt.getdimin()  == B.getdimout() );
    assert( Bt.getdimout() == B.getdimin()  );
    
    assert( b_A.getdimension() == dimension_A );
    assert( b_C.getdimension() == dimension_C );
    assert( x_A.getdimension() == dimension_A );
    assert( x_C.getdimension() == dimension_C );
    
    assert( PAinv.getdimin() == PAinv.getdimout() );
    assert( PCinv.getdimin() == PCinv.getdimout() );
    assert( dimension_A == PAinv.getdimout() );
    assert( dimension_C == PCinv.getdimout() );

    /* Determine the print flags */

    const bool do_print_begin     = print_modulo >= -1;
    const bool do_print_interim   = print_modulo >=  1;
    const bool do_print_restart   = print_modulo >=  0;
    // const bool do_print_breakdown = print_modulo >=  0;
    // const bool do_print_warning   = print_modulo >=  0;
    const bool do_print_finish    = print_modulo >= -1;
    
    /* Build up data */
    
    FloatVector v0_A( dimension_A, 0. );
    FloatVector v1_A( dimension_A, 0. );
    FloatVector w0_A( dimension_A, 0. );
    FloatVector w1_A( dimension_A, 0. );
    FloatVector  z_A( dimension_A, 0. );
    
    FloatVector v0_C( dimension_C, 0. );
    FloatVector v1_C( dimension_C, 0. );
    FloatVector w0_C( dimension_C, 0. );
    FloatVector w1_C( dimension_C, 0. );
    FloatVector  z_C( dimension_C, 0. );
    
    FloatVector vn_A( dimension_A, 0. );
    FloatVector wn_A( dimension_A, 0. );
    FloatVector zn_A( dimension_A, 0. );
    
    FloatVector vn_C( dimension_C, 0. );
    FloatVector wn_C( dimension_C, 0. );
    FloatVector zn_C( dimension_C, 0. );
    
    FloatVector  m_A( dimension_A, 0. );
    FloatVector  m_C( dimension_C, 0. );
    
    FloatVector p_A( dimension_A, 0. );
    FloatVector p_C( dimension_C, 0. );
    
    
    Float mu_A = notanumber;
    Float mu_C = notanumber;
    
    Float gamma = notanumber;
    
    Float eta   = notanumber;
    Float eta_A = notanumber;
    Float eta_C = notanumber;
    
    Float s0 = notanumber;
    Float s1 = notanumber;
    Float c0 = notanumber;
    Float c1 = notanumber;
    
    int max_iteration_count = dimension_A + dimension_C;
    int recent_iteration_count = 0;

    if( do_print_begin ) LOGPRINTF( "(%d/%d)      BEGIN: Block Herzog-Soodhalter CSR\n", recent_iteration_count, max_iteration_count );

    while( recent_iteration_count < max_iteration_count ){
        
        bool restart_condition = ( recent_iteration_count == 0 ) or ( cppsys_restart_on_full_dimension and recent_iteration_count != 0 );
        
        bool residual_seems_small = ( recent_iteration_count != 0 ) and ( absolute(eta) < tolerance );
        
        if( restart_condition or ( residual_seems_small and cppsys_restart_before_finish ) ) UNLIKELY {
            
            // 1
            v0_A.zero(); w0_A.zero(); w1_A.zero();
            v0_C.zero(); w0_C.zero(); w1_C.zero();

            // 2 
            v1_A = b_A - A * x_A - Bt * x_C;
            z_A = PAinv * v1_A;
            v1_C = b_C - B * x_A - C  * x_C;
            z_C = PCinv * v1_C;
            
            // 3
            gamma = std::sqrt( v1_A * z_A + v1_C * z_C );
            
            assert( gamma >= 0. );
            
            v1_A /= gamma;
            v1_C /= gamma;
            z_A /= gamma;
            z_C /= gamma;
            
            // 4 -- 5 
            Float psi_A = z_A * v1_A; 
            Float psi_C = z_C * v1_C;
            mu_A = psi_A; 
            mu_C = psi_C; 
            
            // 6 
            m_A = v1_A;
            m_C = v1_C;

            // 7 
            eta = gamma; 
            eta_A = gamma * std::sqrt( psi_A );
            eta_C = gamma * std::sqrt( psi_C );
            
            s0 = s1 = 0;
            c0 = c1 = 1;
            
            if( do_print_restart ) {
                LOGPRINTF( "(%d/%d)   RESTART: Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble) absolute(eta), (double)(safedouble)tolerance );
                LOGPRINTF( "(%d/%d)            Gamma: %.9le Eta_A %.9le Eta_C %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble)gamma, (double)(safedouble)eta_A, (double)(safedouble)eta_C );
            }

        }
        
        bool residual_is_small = ( absolute(eta) < tolerance );
        
        if( residual_is_small ) UNLIKELY 
            break;

            
        {
            
            // 8 -- 9
            p_A = A * z_A + Bt * z_C;
            p_C = B * z_A + C  * z_C;
 
            Float delta = z_A * p_A + z_C * p_C;

            // 10 
            vn_A = p_A - delta * v1_A - gamma * v0_A;
            vn_C = p_C - delta * v1_C - gamma * v0_C;
            
            // 11 
            zn_A = PAinv * vn_A;
            zn_C = PCinv * vn_C;

            // 12 
            Float gamma_n = std::sqrt( zn_A * vn_A + zn_C * vn_C );
            // Assert( gamma_n > 0., gamma_n );

            // 13 -- 14
            vn_A /= gamma_n; zn_A /= gamma_n;
            vn_C /= gamma_n; zn_C /= gamma_n;

            // 15 -- 18
            Float alpha_0 = c1 * delta - c0 * s1 * gamma;
            assert( alpha_0 * alpha_0 + gamma_n * gamma_n > 0. );
            Float alpha_1 = std::sqrt( alpha_0 * alpha_0 + gamma_n * gamma_n );
            Float alpha_2 = s1 * delta + c0 * c1 * gamma;
            Float alpha_3 = s0 * gamma;
 
            assert( alpha_1 > 0. );

            // 19 
            Float cn = alpha_0 / alpha_1;
            Float sn = gamma_n / alpha_1;
            
            // 20 -- 21 
            Float theta_A = m_A * zn_A;
            Float theta_C = m_C * zn_C; 
            
            Float psi_A = zn_A * vn_A;
            Float psi_C = zn_C * vn_C; 

            // 22 
            m_A = - sn * m_A + cn * vn_A;
            m_C = - sn * m_C + cn * vn_C;

            // 23 
            wn_A = ( z_A - alpha_3 * w0_A - alpha_2 * w1_A ) / alpha_1;
            wn_C = ( z_C - alpha_3 * w0_C - alpha_2 * w1_C ) / alpha_1;

            // 24 
            x_A = x_A + cn * eta * wn_A;
            x_C = x_C + cn * eta * wn_C;
 
            // 25 -- 26
            mu_A = sn * sn * mu_A - 2 * sn * cn * theta_A + cn * cn * psi_A;
            mu_C = sn * sn * mu_C - 2 * sn * cn * theta_C + cn * cn * psi_C;

            // 27 -- 28 
            eta   = -sn * eta;
            eta_A = eta * std::sqrt( psi_A );
            eta_C = eta * std::sqrt( psi_C );


            // update 
            v0_A = v1_A; v1_A = vn_A;
            v0_C = v1_C; v1_C = vn_C;
            w0_A = w1_A; w1_A = wn_A;
            w0_C = w1_C; w1_C = wn_C;
            
            z_A = zn_A;
            z_C = zn_C;
            
            
            gamma = gamma_n;
            
            c0 = c1; c1 = cn;
            s0 = s1; s1 = sn;

            FloatVector r_A = b_A - A * x_A - Bt * x_C;
            FloatVector r_C = b_C - B * x_A - C  * x_C;
            Float r = std::sqrt( r_A * r_A + r_C * r_C );

            bool print_condition = ( print_modulo > 0 and recent_iteration_count % print_modulo == 0 );
            
            if( print_condition and do_print_interim ) {
                LOGPRINTF( "(%d/%d)   INTERIM: Full Residual norm is %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble)r );
                LOGPRINTF( "(%d/%d)            eta_A=%.9le eta_C=%.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble)eta_A, (double)(safedouble)eta_C );
            }

        }

        
        Float recent_deviation = eta;
        
        /* Print information */
        
        bool print_condition = ( print_modulo > 0 and recent_iteration_count % print_modulo == 0 );
        
        if( print_condition and do_print_interim ) UNLIKELY {
            LOGPRINTF( "(%d/%d)   INTERIM: Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble) recent_deviation, (double)(safedouble)tolerance );
            LOGPRINTF( "(%d/%d)            Gamma: %.9le Eta: %.9le\n", recent_iteration_count, max_iteration_count, (double)(safedouble)gamma, (double)(safedouble)eta );
        }
        
        recent_iteration_count++;
        
    }
    
    /* HOW DID WE FINISH ? */
    
    Float recent_deviation = absolute( eta );
        
    if( do_print_finish ) 
        LOGPRINTF( "(%d/%d) %9s: "    "Residual norm is %.9le < %.9le\n", recent_iteration_count, max_iteration_count, recent_deviation < tolerance ? "SUCCESS" : "FAILED", (double)(safedouble)recent_deviation, (double)(safedouble)tolerance );

    return recent_iteration_count;
}
  


