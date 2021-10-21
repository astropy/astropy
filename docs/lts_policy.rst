*******************
LTS Backport Policy
*******************

Starting from astropy 5.0, backports to Long-Term Stable (LTS) releases
will only cover:

* Critical security fixes.
* Bug fixes where backporting is straightforward (e.g., no conflicts).
* More complex fixes are permissible only as resources allow, ideally with
  backports done by the original author, users, or institutions who are on LTS
  and need the fix.
* New "features" that are absolutely necessary to accomplish the above.
* Other critical additions absolutely necessary for users stuck to the LTS.
  In this case, the users themselves (e.g., developers from affected institutions)
  would perform the needed backports or implementation.

This is because LTS lasts for about 2 years. During that time frame, as the
LTS branch diverges from the development branch, backports will become
increasingly difficult due to merge conflicts. When conflicts arise,
automation is not possible, therefore driving up the maintenance cost.
