
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
    const LinearOperator& MassMatrix, // TODO: MassMatrix
    const int max_number_of_candidates,
    const int gram_schmidt_cycles,
    const int max_number_of_purifications,
    //
    const Float threshold_purified
    const Float zero_vector_threshold
) {

    std::vector<FloatVector> nullvectorgallery;
    
    for( int no_candidate = 0; no_candidate < max_number_of_candidates; no_candidate++ )
    {
        
        /* 1. create a new candidate */
        
        LOG << "New candidate" << nl;

        FloatVector candidate( SystemMatrix.getdimin(), 0. ); 
        candidate.random(); 
        candidate.normalize( MassMatrix );
        
        /* 2. reduce the candidate to its nullspace component */
        /*    and orthogonalize against previous nullvectors  */
        
        for( int t = 0; t < max_number_of_purifications; t++ )
        {
            
            FloatVector rhs( SystemMatrix.getdimin(), 0. );
            FloatVector res( SystemMatrix.getdimin(), 0. );
        
            ConjugateResidualMethod solver( SystemMatrix );
            solver.print_modulo        = -1; //candidate.getdimension();
            solver.max_iteration_count = candidate.getdimension();
            solver.threshold           = threshold_residual;

            solver.solve( candidate, rhs );
                            
            
            res = rhs - SystemMatrix * candidate;
            
            LOG << "\t\t\t (eucl) res:       " << res.norm() << nl;
            LOG << "\t\t\t (mass) res:       " << res.norm( MassMatrix ) << nl;
            LOG << "\t\t\t (eucl) x:         " << candidate.norm() << nl;
            LOG << "\t\t\t (mass) x:         " << candidate.norm( MassMatrix ) << nl;
            LOG << "\t\t\t (eucl) Ax:        " << ( SystemMatrix * candidate ).norm() << nl;
            LOG << "\t\t\t (mass) Ax:        " << ( SystemMatrix * candidate ).norm( MassMatrix ) << nl;
            LOG << "\t\t\t (eucl) b - Ax:    " << ( SystemMatrix * candidate - rhs ).norm() << nl;
            LOG << "\t\t\t (mass) b - Ax:    " << ( SystemMatrix * candidate - rhs ).norm( MassMatrix ) << nl;
            
            /* check whether anything is left at all */
            
            Float purified_mass = candidate.norm( MassMatrix );
            LOG << "\t\t\t Filtered mass: " << purified_mass << nl;
            
            if( purified_mass < threshold_purified )
                break;
            
            
            candidate.normalize( MassMatrix );

            LOG << "\t\t\t (norm eucl) x:         " << candidate.norm() << nl;
            LOG << "\t\t\t (norm mass) x:         " << candidate.norm( MassMatrix ) << nl;
            LOG << "\t\t\t (norm eucl) Ax:        " << ( SystemMatrix* candidate ).norm() << nl;
            LOG << "\t\t\t (norm mass) Ax:        " << ( SystemMatrix * candidate ).norm( MassMatrix ) << nl;

        }

        
        /* 3. check whether anything is left at all */
        
        Float purified_mass = candidate.norm( MassMatrix );
        LOG << "\t\t\t Filtered mass: " << purified_mass << nl;
        
        if( purified_mass < threshold_purified ) {
            LOG << "\t\t\t Discard: purified mass too small!" << nl;
            continue;
        }
        


        
        /* 4. Gram-Schmidt the candidate against the previous nullvectors */
        
        for( int s = 0; s < gram_schmidt_cycles; s++ )
        for( const auto& nullvector : nullvectorgallery ) {
            Float alpha = ( MassMatrix * candidate * nullvector ) / ( MassMatrix * nullvector * nullvector );
            candidate = candidate - alpha * nullvector;
        }
        
        Float reduced_mass = candidate.norm( MassMatrix );
        LOG << "\t\t\t Reduced mass: " << reduced_mass << nl;
        
        if( reduced_mass < threshold_reduced ) {
            LOG << "\t\t\t Discard: reduced mass too small" << nl;
            continue;
        }
        
        /* 5. Serious nullspace contender. Normalize and see whether that persists */
        
        candidate.normalize( MassMatrix );
        
        Float residual_mass = ( SystemMatrix * candidate ).norm( MassMatrix );
        
        LOG << "\t\t\t Mass normalized, Numerical residual: " << residual_mass << nl;
        
        if( residual_mass > threshold_residual ) {
            LOG << "\t\t\t Discard: residual too large" << nl;
            continue;
        }
        
        assert( candidate.isfinite() );
        
        /* 6. If we reach here, accept as nullvector */
        
        LOG << "Accept vector: " << nullvectorgallery.size() + 1 << nl;

        nullvectorgallery.push_back( candidate );
    }
    
    
    
    LOG << "How much nullspace are our vectors?" << nl;
    for( const auto& nullvector : nullvectorgallery ) {
        LOGPRINTF( "% 10.5e\t", ( SystemMatrix * nullvector ).norm( MassMatrix ) );
        // LOG << std::showpos << std::scientific << std::setprecision(5) << std::setw(10) << ( SystemMatrix * nullvector ).norm( MassMatrix ) << tab;
    }
    LOG << nl;
    
    LOG << "How orthonormal are our vectors?" << nl;
    for( const auto& nullvector1 : nullvectorgallery ) {
        for( const auto& nullvector2 : nullvectorgallery ) {
            LOGPRINTF( "% 10.5e\t", MassMatrix * nullvector1 * nullvector2 );
            // LOG << std::showpos << std::scientific << std::setprecision(5) << std::setw(10) << MassMatrix * nullvector1 * nullvector2 << tab;
        }
        LOG << nl;
    }

}

