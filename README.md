# ARKODE-stepper-template

This repository contains a "template" for a new time stepping module to be built on top of ARKODE.  As we do not have a clear API for how ARKODE's time-stepping modules must interact with the underlying ARKODE infrastructure (more specifically, time-stepping modules currently "reach into" the main ARKODE memory structure to access relevant infrastructure components), this template repository can serve as a placeholder for users to build from when creating new ARKODE time stepping modules.  The directory layout of this repository matches that of the main SUNDIALS repository, allowing its contents to be more easily patched into SUNDIALS in the case that this will eventually be propagated upstream.

## Template contributions ##

Bug fixes or minor changes to this template are welcomed -- the preferred mechanism is via a pull request.

## Upstream contributions ##

New time-stepping modules for ARKODE cannot be submitted via pull requests to either this repository or the main [SUNDIALS GitHub repository](https://github.com/LLNL/sundials).  Anyone who desires to make these larger-scale upstream contributions should contact the SUNDIALS developers team directly.

## Authors ##

This template repository has been developed by [Daniel R. Reynolds](https://github.com/drreynolds).

## License ##

This template repository is released under the BSD 3-clause license. See the [LICENSE](./LICENSE) file for details. All new contributions must be made under the BSD 3-clause license.
