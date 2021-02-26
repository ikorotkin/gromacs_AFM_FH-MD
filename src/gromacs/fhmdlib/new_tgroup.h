#ifndef FHMD_NEW_TGROUP_H
#define FHMD_NEW_TGROUP_H

#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/fhmdlib/data_structures.h"

real sum_FHMD_ekin(t_grpopts *opts, gmx_ekindata_t *ekind, real *dekindlambda,
              gmx_bool bEkinFullStep, gmx_bool bScaleEkin, gmx_int64_t  step,FHMD *fh);
/* Sum the group ekins into total ekin and calc temp per group, calc temp per FHMD bin,
 * return total temperature.
 */

#endif
