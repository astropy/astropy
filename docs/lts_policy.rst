*******************
LTS Backport Policy
*******************

Starting from astropy 5.0, backports to Long-Term Stable (LTS) releases
will only cover:

* Critical bug fixes.
* Critical security fixes.
* Special request from pipelines stuck to the LTS. For this use case,
  preferably the developers from affected pipelines would help with the
  backports.

This is because LTS lasts for about 2 years. During that time frame, as the
LTS branch diverges from the development branch, backports will become
increasingly difficult due to merge conflicts. When conflicts arise,
automation is not possible, therefore driving up the maintenance cost.
