#ifndef INCLUDEGUARD_SOLVER_NULLSPACE
#define INCLUDEGUARD_SOLVER_NULLSPACE

#include <vector>

#include "../base/include.hpp"
#include "../operators/floatvector.hpp"
#include "../operators/linearoperator.hpp"


std::vector<FloatVector> computeNullspace(
    const LinearOperator& SystemMatrix,
    const LinearOperator& MassMatrix,
    const LinearOperator& ResidualMassMatrix,
    const int max_number_of_candidates,
    //
    const Float mass_threshold_for_small_residual, 
    const Float mass_threshold_for_small_vectors,
    const std::function<void(FloatVector&)>& purifier
);

#endif
