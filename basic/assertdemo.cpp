
#include "debug.hpp"

/* A dummy function to make the backtrace more interesting. */
void
dummy_function (void)
{
  Assert(false);
}

int
main (void)
{
  dummy_function ();
  return 0;
}

