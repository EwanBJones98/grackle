#include <stdio.h>
#include "grackle_types.h"

// the following macros are auto-generated:
#define AUTO_VERSION "3.3.dev1"
#define AUTO_BRANCH "pencil_beam_test"
#define AUTO_REVISION "7cf6bce8adbc405fb73595b96754833b26839172"

// test that ensures that all macros were correctly defined:
#if !(defined(AUTO_VERSION) && defined(AUTO_BRANCH) && defined(AUTO_REVISION))
#error "Something went wrong while auto-generating macros"
#endif

grackle_version get_grackle_version(void) {
  grackle_version out;
  out.version = AUTO_VERSION;
  out.branch = AUTO_BRANCH;
  out.revision = AUTO_REVISION;
  return out;
}
