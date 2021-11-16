
<!-- NEWS.md is generated from NEWS.Rmd. Please edit NEWS.Rmd file -->

# R2BEAT 1.0.4

## Major changes

-   Added the function ‘eval\_2stage’ for the evaluation of the
    two-stage design by means of simulation

-   Added the functions ‘sensitivity\_min\_SSU’ and
    ‘sensitivity\_min\_SSU2’ to give hints about how to choose the value
    of ‘minimum’ parameter (minimum number of SSUs to be selected in
    each PSU). The first is to be used when a sampling frame of SSUs is
    available, the second when it is not available. These functions are
    in substitution of function ‘sensitivity’

-   Added the function ‘select\_PSU’ to select first stage units (in
    substitution of ‘StratSel’, eliminated)

-   Added the functions ‘prepareInputToAllocation1’ and
    ‘prepareInputToAllocation2’ to produce the inputs required by the
    allocation step (the first in case no previous round of the survey
    is available, the second when at least one round of the survey is
    available)

# R2BEAT 1.0.3

## Major changes

-   Added a function (‘select\_SSU’) for the selection of Secondary
    Stage Units

-   Added a function (‘sensitivity’) for the fine tuning of parameters

# R2BEAT 1.0.2

## Major changes

-   Added two functions (‘input\_to\_beat.2st\_1’ and
    ‘input\_to\_beat.2st\_2’) to prepare inputs for the allocation step

-   Added function ‘StratSel’ executing the PSU selection step

-   Fixed some bugs

# R2BEAT 1.0.1
