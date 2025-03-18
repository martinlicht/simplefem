
#include <vector>

#include "../base/include.hpp"
#include "../operators/floatvector.hpp"
#include "../operators/linearoperator.hpp"
#include "iterativesolver.hpp"

#include "nullspace.hpp"

std::vector<FloatVector> computeNullspace(
    const LinearOperator& SystemMatrix,
    const LinearOperator& MassMatrix,
    const LinearOperator& ResidualMassMatrix,
    const int max_number_of_candidates,
    //
    const Float mass_threshold_for_small_residual, 
    const Float mass_threshold_for_small_vectors,
    std::function<void(FloatVector&)> purifier
) { 

    const int max_number_of_purifications = 2;

    std::vector<FloatVector> nullvectorgallery;
    
    for( int candidate_number = 0; candidate_number < max_number_of_candidates; candidate_number++ )
    {
        
        /* Draw a random vector with unit mass norm as a candidate */

        FloatVector candidate( SystemMatrix.getdimin(), 0. ); 
        candidate.random(); 
        candidate.normalize(MassMatrix);
        assert( candidate.is_finite() );
        
        /* Orthogonalize that vector against previous nullspace vectors */

        {
            for( int s = 0; s < 2; s++ )
            for( const auto& nullvector : nullvectorgallery ) 
            {
                Float alpha = (MassMatrix*candidate*nullvector) / (MassMatrix*nullvector*nullvector);
                candidate = candidate - alpha * nullvector;
            }   
        }
        
        /* Re-normalize the candidate once again */

        candidate.normalize(MassMatrix);

        assert( candidate.is_finite() );
        
        /* Reduce the candidate to its nullspace component, orthogonalize, normalize */
        
        bool candidate_to_be_discarded = false;
        
        for( int t = 0; t < max_number_of_purifications and not candidate_to_be_discarded; t++ )
        {
            
            candidate.normalize(MassMatrix);

            assert( candidate.is_finite() );
            
            LOG << "\t\t\t (eucl) x:      " << candidate.norm() << nl;
            LOG << "\t\t\t (mass) x:      " << candidate.norm( MassMatrix ) << nl;
            LOG << "\t\t\t (eucl) Ax:     " << ( SystemMatrix * candidate ).norm() << nl;
            LOG << "\t\t\t (mass) Ax:     " << ( SystemMatrix * candidate ).norm( ResidualMassMatrix ) << nl;

            // PreconditionedConjugateResidualMethod PCRM(SystemMatrix,ApproxInverseSystem);
            // PCRM.verbosity = IterativeSolver::VerbosityLevel::startandfinish;

            // FloatVector solution( candidate.getdimension(), 0. );
            // PCRM.solve( solution, SystemMatrix * candidate );
            // candidate -= solution;

            purifier(candidate);
            
            /* Orthogonalize that candidate once again against nullspace vectors */
            /* Discard if nothing new is found, else re-normalize */
            
            for( int s = 0; s < 2; s++ )
            for( const auto& nullvector : nullvectorgallery ) {
                Float alpha = (MassMatrix*candidate*nullvector) / (MassMatrix*nullvector*nullvector);
                candidate = candidate - alpha * nullvector;
            }
            
            Float orthogonalized_candidate_mass = candidate.norm(MassMatrix);
            LOG << "\t\t\t Reduced candidate orthogonalized mass: " << orthogonalized_candidate_mass << nl;

            if( orthogonalized_candidate_mass < mass_threshold_for_small_vectors ) {
                candidate_to_be_discarded = true;
            }
            
            LOG << "\t\t\t (eucl) x:      " << candidate.norm() << nl;
            LOG << "\t\t\t (mass) x:      " << candidate.norm( MassMatrix ) << nl;
            LOG << "\t\t\t (eucl) Ax:     " << ( SystemMatrix * candidate ).norm() << nl;
            LOG << "\t\t\t (mass) Ax:     " << ( SystemMatrix * candidate ).norm( ResidualMassMatrix ) << nl;
            

        }
        
        /* Orthogonalize that candidate once again against nullspace vectors */
        
        LOG << "\t\t\t Reduced candidate non-orthogonalized mass: " << candidate.norm(MassMatrix) << nl;

        for( int s = 0; s < 2; s++ )
        for( const auto& nullvector : nullvectorgallery ) {
            Float alpha = ( MassMatrix * candidate * nullvector ) / ( MassMatrix * nullvector * nullvector );
            candidate = candidate - alpha * nullvector;
        }
        
        /* Discard if nothing new is found */

        const Float orthogonalized_candidate_mass = candidate.norm(MassMatrix);
        
        LOG << "\t\t\t Reduced candidate orthogonalized mass: " << orthogonalized_candidate_mass << nl;

        if( candidate_to_be_discarded or orthogonalized_candidate_mass < mass_threshold_for_small_vectors ) {
            LOGPRINTF("\n\t\t\t !!!!! Discard vector %d/%d because mass is too small!\n\n", candidate_number+1, max_number_of_candidates );
            continue;
        }
        
        /* Discard if insufficiently nullspace */ 
        
        Float residual_mass = ( SystemMatrix * candidate ).norm(ResidualMassMatrix);

        Float residual_ratio = residual_mass / orthogonalized_candidate_mass;
        
        LOG << "\t\t\t Numerical residual ratio: " << residual_ratio << nl;
        
        if( residual_ratio > mass_threshold_for_small_residual ) {
            LOGPRINTF("\n\t\t\t !!!!!Discard vector %d/%d because not nullspace enough!\n\n", candidate_number+1, max_number_of_candidates );
            continue;
        }
        
        /* Normalize and accept */

        candidate /= orthogonalized_candidate_mass;
        
        assert( candidate.is_finite() );
        
        LOG << "Accept vector: " << nullvectorgallery.size() + 1 << nl;
    
        
        nullvectorgallery.push_back( candidate );
    }
    
    
    
    LOG << "How much nullspace are our vectors?" << nl;
    for( const auto& nullvector : nullvectorgallery ) {
        Float residual_mass = ( SystemMatrix * nullvector ).norm( ResidualMassMatrix );
        
        Assert( residual_mass < mass_threshold_for_small_vectors, residual_mass, mass_threshold_for_small_vectors );
        
        LOGPRINTF( "% 10.5le\t", (double)(safedouble)residual_mass );
    }
    LOG << nl;
    
    LOG << "How orthonormal are our vectors?" << nl;
    for( int n1 = 0; n1 < nullvectorgallery.size(); n1++ ) {
        for( int n2 = 0; n2 < nullvectorgallery.size(); n2++ ) {
            auto nullvector1 = nullvectorgallery[n1];
            auto nullvector2 = nullvectorgallery[n2];
            
            Float mass_prod = MassMatrix * nullvector1 * nullvector2;
            
            LOGPRINTF( "% 10.5le\t", (double)(safedouble)mass_prod );
            
            if( n1 != n2 ) 
                assert( is_numerically_small( mass_prod ) );
            else
                assert( is_numerically_one( mass_prod ) );
        }
        LOG << nl;
    }
    
    
    
    return nullvectorgallery;

}

