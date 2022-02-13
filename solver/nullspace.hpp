
#include <ostream>
#include <fstream>
// #include <iomanip> // TODO: remove io manipulators 

#include "../basic.hpp"
#include "../operators/floatvector.hpp"
#include "../operators/linearoperator.hpp"
#include "../solver/sparsesolver.hpp"
#include "../solver/iterativesolver.hpp"


std::vector<FloatVector> computeNullspace(
    const LinearOperator& SystemMatrix,
    const LinearOperator& mass, // TODO: MassMatrix
    const int max_number_of_candidates,
    const int gram_schmidt_cycles,
    const int max_number_of_purifications,
    //
    const Float zero_vector_threshold
) {

    std::vector<FloatVector> nullvectorgallery;
    
    for( int no_candidate = 0; no_candidate < max_number_of_candidates; no_candidate++ )
    {
        
        FloatVector candidate( SystemMatrix.getdimin(), 0. ); 
        candidate.random(); 
        candidate.normalize(mass);
        
        {
            for( int s = 0; s < gram_schmidt_cycles; s++ )
            for( const auto& nullvector : nullvectorgallery ) {
                Float alpha = (mass*candidate*nullvector) / (mass*nullvector*nullvector);
                candidate = candidate - alpha * nullvector;
            }
            
            Float reduced_mass = candidate.norm(mass);
            LOG << "\t\t\t Preprocessed mass: " << reduced_mass << nl;
            
            if( reduced_mass < zero_vector_threshold ) {
                LOG << "**** The candidate already has very small mass" << nl;
//                                 continue;
            }
        }
        
        
        /* reduce the candidate to its nullspace component */
        {
            FloatVector rhs( SystemMatrix.getdimin(), 0. );
        
            FloatVector residual( rhs );
            
            for( int t = 0; t < max_number_of_purifications; t++ )
            {
                
                ConjugateResidualMethod solver( SystemMatrix );
                solver.print_modulo        = candidate.getdimension() / 20;
                solver.max_iteration_count = candidate.getdimension();
                solver.threshold           = sqrt( machine_epsilon );

                solver.solve( candidate, rhs );
                               
                
                residual = rhs - SystemMatrix * candidate;
                
                LOG << "\t\t\t (eucl) delta:     " << ( residual - rhs + SystemMatrix * candidate ).norm() << nl;
                LOG << "\t\t\t (mass) delta:     " << ( residual - rhs + SystemMatrix * candidate ).norm( mass ) << nl;
                LOG << "\t\t\t (eucl) res:       " << residual.norm() << nl;
                LOG << "\t\t\t (mass) res:       " << residual.norm( mass ) << nl;
                LOG << "\t\t\t (eucl) x:         " << candidate.norm() << nl;
                LOG << "\t\t\t (mass) x:         " << candidate.norm( mass ) << nl;
                LOG << "\t\t\t (eucl) Ax:        " << ( SystemMatrix * candidate ).norm() << nl;
                LOG << "\t\t\t (mass) Ax:        " << ( SystemMatrix * candidate ).norm( mass ) << nl;
                LOG << "\t\t\t (eucl) b - Ax:    " << ( SystemMatrix * candidate - rhs ).norm() << nl;
                LOG << "\t\t\t (mass) b - Ax:    " << ( SystemMatrix * candidate - rhs ).norm( mass ) << nl;
                
                
                candidate.normalize( mass );

                assert( candidate.isfinite() );
                
                LOG << "\t\t\t (norm eucl) x:         " << candidate.norm() << nl;
                LOG << "\t\t\t (norm mass) x:         " << candidate.norm( mass ) << nl;
                LOG << "\t\t\t (norm eucl) Ax:        " << ( SystemMatrix* candidate ).norm() << nl;
                LOG << "\t\t\t (norm mass) Ax:        " << ( SystemMatrix * candidate ).norm( mass ) << nl;
                
                
//                                 FloatVector zero( candidate.getdimension(), 0. );
//                                 FloatVector residual( candidate.getdimension(), 0. );
//                                 
//                                 ConjugateResidualSolverCSR( 
//                                     zero.getdimension(), 
//                                     zero.raw(), 
//                                     candidate.raw(), 
//                                     SystemMatrix.getA(), SystemMatrix.getC(), SystemMatrix.getV(),
//                                     residual.raw(),
//                                     1000 * machine_epsilon,
//                                     0
//                                 );
//                                 candidate = residual;
//                                 candidate.normalize( mass );
                
                
                
                
            }
        }
        
        
        /* Gram-Schmidt */
        
        for( int s = 0; s < gram_schmidt_cycles; s++ )
        for( const auto& nullvector : nullvectorgallery ) {
            Float alpha = (mass*candidate*nullvector) / (mass*nullvector*nullvector);
            candidate = candidate - alpha * nullvector;
        }
        
        Float reduced_mass = candidate.norm(mass);
        LOG << "\t\t\t Reduced mass: " << reduced_mass << nl;
        
        if( reduced_mass < zero_vector_threshold ) {
            LOG << "!!!!!!!!!!!!!Discard vector because mass is too small!" << nl;
            continue;
        }
        
        candidate.normalize(mass);
        
        Float residual_mass = ( SystemMatrix * candidate ).norm(mass);
        
        LOG << "\t\t\t Numerical residual: " << residual_mass << nl;
        
        if( residual_mass > 1e-6 ) {
            LOG << "!!!!!!!!!!!!!Discard vector because not nullspace enough!" << nl;
            continue;
        }
        
        assert( candidate.isfinite() );
        
        LOG << "Accept vector: " << nullvectorgallery.size() + 1 << nl;
    
        
        nullvectorgallery.push_back( candidate );
    }
    
    
    
    LOG << "How much nullspace are our vectors?" << nl;
    for( const auto& nullvector : nullvectorgallery ) {
        LOGPRINTF( "% 10.5e\t", ( SystemMatrix * nullvector ).norm(mass) );
        // LOG << std::showpos << std::scientific << std::setprecision(5) << std::setw(10) << ( SystemMatrix * nullvector ).norm(mass) << tab;
    }
    LOG << nl;
    
    LOG << "How orthonormal are our vectors?" << nl;
    for( const auto& nullvector1 : nullvectorgallery ) {
        for( const auto& nullvector2 : nullvectorgallery ) {
            LOGPRINTF( "% 10.5e\t", mass * nullvector1 * nullvector2 );
            // LOG << std::showpos << std::scientific << std::setprecision(5) << std::setw(10) << mass * nullvector1 * nullvector2 << tab;
        }
        LOG << nl;
    }

}

