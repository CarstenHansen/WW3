! -------------------------------------------------------------------- !
! WAVEWATCH III ww3_multi.nml - multi-grid model                       !
! -------------------------------------------------------------------- !


! -------------------------------------------------------------------- !
! Define top-level model parameters via DOMAIN_NML namelist
! -------------------------------------------------------------------- !
&DOMAIN_NML
  DOMAIN%NRGRD	= 5
  DOMAIN%FLGHG1 = F
  DOMAIN%FLGHG2 = F
  DOMAIN%START  = '__SIMSTART__'
  DOMAIN%STOP   = '__SIMSTOP__'
/

! -------------------------------------------------------------------- !
! Define each input grid via the INPUT_GRID_NML namelist
! -------------------------------------------------------------------- !
&INPUT_GRID_NML
/

! -------------------------------------------------------------------- !
! Define each model grid via the MODEL_GRID_NML namelist
! -------------------------------------------------------------------- !
! TODO Profiling: Use switch MPRF, see https://\
!   polar.ncep.noaa.gov/waves/workshop/pdfs/WW3-workshop-exercises-day4-profiling.pdf
!
&MODEL_GRID_NML
  MODEL(1)%NAME                  = 'NUUP_5400'
  MODEL(1)%FORCING%WINDS         = 'native'
  ! MODEL(1)%FORCING%ICE_CONC      = 'native'
  ! MODEL(1)%FORCING%ICE_PARAM1    = 'native'
  ! MODEL(1)%FORCING%ICE_PARAM2    = 'H'
  ! MODEL(1)%FORCING%ICE_PARAM3    = 'H'
  ! MODEL(1)%FORCING%ICE_PARAM4    = 'H'
  MODEL(1)%RESOURCE%COMM_FRAC = 0.00,1.00
  MODEL(2)%NAME                  = 'NUUP_1800'
  MODEL(2)%FORCING%WINDS         = 'native'
  ! MODEL(2)%FORCING%ICE_CONC      = 'native'
  ! MODEL(2)%FORCING%ICE_PARAM1    = 'native'
  ! MODEL(2)%FORCING%ICE_PARAM2    = 'H'
  ! MODEL(2)%FORCING%ICE_PARAM3    = 'H'
  ! MODEL(2)%FORCING%ICE_PARAM4    = 'H'
  MODEL(2)%RESOURCE%COMM_FRAC = 0.00,1.00
  MODEL(3)%NAME                  = 'NUUP_600W'
  MODEL(3)%FORCING%WINDS         = 'native'
  ! MODEL(3)%FORCING%ICE_CONC      = 'native'
  ! MODEL(3)%FORCING%ICE_PARAM1    = 'native'
  ! MODEL(3)%FORCING%ICE_PARAM2    = 'H'
  ! MODEL(3)%FORCING%ICE_PARAM3    = 'H'
  ! MODEL(3)%FORCING%ICE_PARAM4    = 'H'
  MODEL(3)%RESOURCE%COMM_FRAC = 0.00,1.00
  MODEL(3)%RESOURCE%RANK_ID      = 3
  MODEL(4)%NAME                  = 'NUUP_600E'
  MODEL(4)%FORCING%WINDS         = 'native'
  ! MODEL(4)%FORCING%ICE_CONC      = 'native'
  ! MODEL(4)%FORCING%ICE_PARAM1    = 'native'
  ! MODEL(4)%FORCING%ICE_PARAM2    = 'H'
  ! MODEL(4)%FORCING%ICE_PARAM3    = 'H'
  ! MODEL(4)%FORCING%ICE_PARAM4    = 'H'
  MODEL(4)%RESOURCE%COMM_FRAC = 0.00,1.00
  MODEL(4)%RESOURCE%RANK_ID      = 3
  MODEL(5)%NAME                  = 'NUUP_200W'
  MODEL(5)%FORCING%WINDS         = 'native'
  ! MODEL(5)%FORCING%ICE_CONC      = 'native'
  ! MODEL(5)%FORCING%ICE_PARAM1    = 'native'
  ! MODEL(5)%FORCING%ICE_PARAM2    = 'H'
  ! MODEL(5)%FORCING%ICE_PARAM3    = 'H'
  ! MODEL(5)%FORCING%ICE_PARAM4    = 'H'
  MODEL(5)%RESOURCE%COMM_FRAC = 0.00,1.00
  MODEL(5)%RESOURCE%RANK_ID      = 4
  ! Defaults:
  !     MODEL(I)%RESOURCE%COMM_FRAC    = 0.00,1.00
  !     MODEL(I)%RESOURCE%RANK_ID      = I
  !     MODEL(I)%RESOURCE%GROUP_ID     = 1
/
  !
  ! &TIMESTEPS_NML                    New:
  ! 200:
  !   TIMESTEPS%DTMAX        =   20.  15. 
  !   TIMESTEPS%DTXY         =   5.   5.  
  !   TIMESTEPS%DTKTH        =   10.  5. 
  !   TIMESTEPS%DTMIN        =   5.   5.  
  !                                       
  ! 600:                                  
  !   TIMESTEPS%DTMAX        =   60.  60. 
  !   TIMESTEPS%DTXY         =   20.  15. 
  !   TIMESTEPS%DTKTH        =   30.  15. 
  !   TIMESTEPS%DTMIN        =   5.   5.  
  !                                       
  ! 1800:  No output!                         
  !   TIMESTEPS%DTMAX        =   200. 180.
  !   TIMESTEPS%DTXY         =   100. 60.
  !   TIMESTEPS%DTKTH        =   100. 60.
  !   TIMESTEPS%DTMIN        =   5.   5.  
  !                                       
  ! 5400:  No output!       
  !   TIMESTEPS%DTMAX        =   600. 540.
  !   TIMESTEPS%DTXY         =   150. 180.
  !   TIMESTEPS%DTKTH        =   150. 180.
  !   TIMESTEPS%DTMIN        =   5.   5.  
  !     
/


! -------------------------------------------------------------------- !
! Define the output types point parameters via OUTPUT_TYPE_NML namelist
! -------------------------------------------------------------------- !
&OUTPUT_TYPE_NML
  ALLTYPE%FIELD%LIST       = 'HS WND DPT DIR'
  ! Discarding ICE ICE1
  ITYPE(1)%POINT%NAME       = 'points_NUUP_5400'
  ITYPE(1)%POINT%FILE       = 'points_NUUP_5400.list'
/
! ALLTYPE%PARTITION        = 0 999 1 0 999 1 F

! -------------------------------------------------------------------- !
! Define output dates via OUTPUT_DATE_NML namelist
! -------------------------------------------------------------------- !
&OUTPUT_DATE_NML
  ALLDATE%FIELD          = '__SIMSTART__' '1200' '__SIMSTOP__'
  ALLDATE%POINT          = '__SIMSTART__' '1200' '__SIMSTOP__'
  ALLDATE%RESTART        = '__SIMSTART__' '__HOTDELTA__' '__SIMSTOP__'
/
! ALLDATE%PARTITION      = '20090525 000000' '3600' '20090526 000000'
! -------------------------------------------------------------------- !
! Define homogeneous input via HOMOG_COUNT_NML and HOMOG_INPUT_NML namelist
! -------------------------------------------------------------------- !

! Switch IC5, with namelist option IC5VEMOD=3:
! suggest eta=10.0 kg/m^3/s == parameter 'IC2' (IC3 and IC4 are read but not used!).
! Rogers et al., 2021, estimates a much higher coefficient,
!  HOMOG_INPUT(1)%VALUE1      = 23.4
! and also achieves a best fit with k_i multiplied by a frequency measure to power 1
! but they do not consider IS2 anelastic scattering that also dissipates energy
! at the highest frequencies ...

&HOMOG_COUNT_NML
  HOMOG_COUNT%N_MOV                =  3
/

!  HOMOG_COUNT%N_IC2                =  1 ! Make this work in ww3_multi !
!  HOMOG_COUNT%N_IC3                =  1
!  HOMOG_COUNT%N_IC4                =  1

&HOMOG_INPUT_NML
  HOMOG_INPUT(1)%NAME        = 'IC2'
  HOMOG_INPUT(1)%VALUE1      = 6.0

  HOMOG_INPUT(2)%NAME        = 'IC3'
  HOMOG_INPUT(2)%VALUE1      = 917.0

  HOMOG_INPUT(3)%NAME        = 'IC4'
  HOMOG_INPUT(3)%VALUE1      = 1.0
/

! -------------------------------------------------------------------- !
! WAVEWATCH III - end of namelist                                      !
! -------------------------------------------------------------------- !
