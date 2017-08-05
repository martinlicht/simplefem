
#include "attest.core.hpp"

// Define macros that call the custom attest 
// and which are compiled unless NDEBUG is set.

#ifdef NDEBUG
# define attest(EX)     ((void)0)
# define warning(EX)    ((void)0)
# define shout(msg)     ((void)0)
# define raiserror(msg) ((void)0)
# define unreachable()  ((void)0)
#else
# define attest(EX)     enforce_attest(EX)
# define warning(EX)    enforce_warning(EX)
# define shout(msg)     enforce_shout(msg)
# define raiserror(msg) enforce_raiseerror(msg)
# define unreachable()  enforce_unreachable()
#endif


