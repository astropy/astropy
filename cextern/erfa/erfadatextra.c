#include "erfa.h"
#include "erfaextra.h"

static eraLEAPSECOND *changes;
static int NDAT = 0;


int eraGetLeapSeconds(eraLEAPSECOND **leapseconds)
{
    if (NDAT == 0) {
        double delat;
        int stat = eraDat(2000, 1, 1, 0., &delat);
        if (stat != 0 || NDAT == 0) {
            return -1;
        }
    }
    *leapseconds = changes;
    return NDAT;
}

void eraSetLeapSeconds(eraLEAPSECOND *leapseconds, int count)
{
    changes = leapseconds;
    NDAT = count;
}

/*
 * For internal use in dat.c
 */
int eraDatini(const eraLEAPSECOND *builtin, int n_builtin,
              eraLEAPSECOND **leapseconds)
{
    if (NDAT == 0) {
        eraSetLeapSeconds((eraLEAPSECOND *)builtin, n_builtin);
    }
    *leapseconds = changes;
    return NDAT;
}
