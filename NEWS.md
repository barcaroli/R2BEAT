
<!-- NEWS.md is generated from NEWS.Rmd. Please edit NEWS.Rmd file -->

# R2BEAT 1.0.5

## Major changes

- Added the function ‘select_PSU2’ to be used in alternative a
  ‘select_PSU’ when there is a quasi uniform distribution of the measure
  of size of the PSUs (this is the case of Enumeration Areas, in
  contrast to Municipalities): this function ensures a complete
  alignment of the number of selected PSUs and of the SSUs to be
  selected to corrispondent numbers of allocated PSUs and SSUs

- Added the function ‘build_dummy_variables’ to add dummy variables to
  the sampling frame (this is necessary if we want to differentiate
  precision constraints in the same domain level, for instance in the
  different regions)

# R2BEAT 1.0.4

## Major changes

- Added the function ‘eval_2stage’ for the evaluation of the two-stage
  design by means of simulation

- Added the functions ‘sensitivity_min_SSU’ and ‘sensitivity_min_SSU2’
  to give hints about how to choose the value of ‘minimum’ parameter
  (minimum number of SSUs to be selected in each PSU). The first is to
  be used when a sampling frame of SSUs is available, the second when it
  is not available. These functions are in substitution of function
  ‘sensitivity’

- Added the function ‘select_PSU’ to select first stage units (in
  substitution of ‘StratSel’, eliminated)

- Added the functions ‘prepareInputToAllocation1’ and
  ‘prepareInputToAllocation2’ to produce the inputs required by the
  allocation step (the first in case no previous round of the survey is
  available, the second when at least one round of the survey is
  available)

# R2BEAT 1.0.3

## Major changes

- Added a function (‘select_SSU’) for the selection of Secondary Stage
  Units

- Added a function (‘sensitivity’) for the fine tuning of parameters

# R2BEAT 1.0.2

## Major changes

- Added two functions (‘input_to_beat.2st_1’ and ‘input_to_beat.2st_2’)
  to prepare inputs for the allocation step

- Added function ‘StratSel’ executing the PSU selection step

- Fixed some bugs

# R2BEAT 1.0.1
