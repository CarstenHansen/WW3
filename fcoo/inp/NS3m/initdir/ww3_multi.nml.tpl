! -------------------------------------------------------------------- !
! WAVEWATCH III - ww3_multi.nml - multi-grid model                     !
! -------------------------------------------------------------------- !


! -------------------------------------------------------------------- !
! Define top-level model parameters via DOMAIN_NML namelist
!
! * IOSTYP defines the output server mode for parallel implementation.
!             0 : No data server processes, direct access output from
!                 each process (requires true parallel file system).
!             1 : No data server process. All output for each type
!                 performed by process that performs computations too.
!             2 : Last process is reserved for all output, and does no
!                 computing.
!             3 : Multiple dedicated output processes.
!
! * namelist must be terminated with /
! * definitions & defaults:
!     DOMAIN%NRINP  =  0  ! Number of grids defining input fields.
!     DOMAIN%NRGRD  =  1  ! Number of wave model grids.
!     DOMAIN%UNIPTS =  F  ! Flag for using unified point output file.
!     DOMAIN%IOSTYP =  1  ! Output server type
!     DOMAIN%UPPROC =  F  ! Flag for dedicated process for unified point output.
!     DOMAIN%PSHARE =  F  ! Flag for grids sharing dedicated output processes.
!     DOMAIN%FLGHG1 =  F  ! Flag for masking computation in two-way nesting
!     DOMAIN%FLGHG2 =  F  ! Flag for masking at printout time
!     DOMAIN%START  = '19680606 000000'  ! Start date for the entire model
!     DOMAIN%STOP   = '19680607 000000'  ! Stop date for the entire model
! -------------------------------------------------------------------- !
&DOMAIN_NML
  DOMAIN%NRINP  = 2 ! NS3 and BAL3, or 3-nmi NSBaltic extended forcing
  DOMAIN%NRGRD  = 2 ! NS3 and BAL3
  DOMAIN%UNIPTS = F ! Change to 'T' after tests
  DOMAIN%START  = '__SIMSTART__'
  DOMAIN%STOP   = '__SIMSTOP__'
/



! -------------------------------------------------------------------- !
! Define each input grid via the INPUT_GRID_NML namelist
!
! * index I must match indexes from 1 to DOMAIN%NRINP
! * INPUT(I)%NAME must be set for each active input grid I
!
! * namelist must be terminated with /
! * definitions & defaults:
!     INPUT(I)%NAME                  = 'unset'
!     INPUT(I)%FORCING%WATER_LEVELS  = F
!     INPUT(I)%FORCING%CURRENTS      = F
!     INPUT(I)%FORCING%WINDS         = F
!     INPUT(I)%FORCING%ATM_MOMENTUM  = F
!     INPUT(I)%FORCING%AIR_DENSITY   = F
!     INPUT(I)%FORCING%ICE_CONC      = F
!     INPUT(I)%FORCING%ICE_PARAM1    = F
!     INPUT(I)%FORCING%ICE_PARAM2    = F
!     INPUT(I)%FORCING%ICE_PARAM3    = F
!     INPUT(I)%FORCING%ICE_PARAM4    = F
!     INPUT(I)%FORCING%ICE_PARAM5    = F
!     INPUT(I)%FORCING%MUD_DENSITY   = F
!     INPUT(I)%FORCING%MUD_THICKNESS = F
!     INPUT(I)%FORCING%MUD_VISCOSITY = F
!     INPUT(I)%ASSIM%MEAN            = F
!     INPUT(I)%ASSIM%SPEC1D          = F
!     INPUT(I)%ASSIM%SPEC2D          = F
! -------------------------------------------------------------------- !
&INPUT_GRID_NML
  INPUT(1)%NAME                  = 'NS3'
  INPUT(1)%FORCING%WINDS         = T
  INPUT(2)%NAME                  = 'BAL3'
  INPUT(2)%FORCING%WINDS         = T
  INPUT(2)%FORCING%ICE_CONC      = T
  INPUT(2)%FORCING%ICE_PARAM1    = T
  INPUT(2)%FORCING%ICE_PARAM2    = T
  INPUT(2)%FORCING%ICE_PARAM3    = T
  INPUT(2)%FORCING%ICE_PARAM4    = T
  INPUT(2)%FORCING%ICE_PARAM5    = F
/


! -------------------------------------------------------------------- !
! Define each model grid via the MODEL_GRID_NML namelist
!
! * index I must match indexes from 1 to DOMAIN%NRGRD
! * MODEL(I)%NAME must be set for each active model grid I
! * FORCING can be set as :
!    - 'no'          : This input is not used.
!    - 'native'      : This grid has its own input files, e.g. grid
!                      grdX (mod_def.grdX) uses ice.grdX.
!    - 'INPUT%NAME'  : Take input from the grid identified by
!                      INPUT%NAME.
! * RESOURCE%RANK_ID : Rank number of grid (internally sorted and reassigned).
! * RESOURCE%GROUP_ID : Group number (internally reassigned so that different
!                                     ranks result in different group numbers).
! * RESOURCE%COMM_FRAC : Fraction of communicator (processes) used for this grid.
! * RESOURCE%BOUND_FLAG : Flag identifying dumping of boundary data used by this
!                         grid. If true, the file nest.MODID is generated.
!
! * Limitations relevant to irregular (curvilinear) grids:
!   1) Equal rank is not supported when one or more is an irregular
!       grid. Use non-equal rank instead. (see wmgridmd.ftn)
!   2) Non-native input grids: feature is not supported when either
!      an input grid or computational grids is irregular.
!      (see wmupdtmd.ftn)
!   3) Irregular grids with unified point output: This is supported
!      but the feature has not been verified for accuracy.
!      (see wmiopomd.ftn)
!
! * namelist must be terminated with /
! * definitions & defaults:
!     MODEL(I)%NAME                  = 'unset'
!     MODEL(I)%FORCING%WATER_LEVELS  = 'no'
!     MODEL(I)%FORCING%CURRENTS      = 'no'
!     MODEL(I)%FORCING%WINDS         = 'no'
!     MODEL(I)%FORCING%ATM_MOMENTUM  = 'no'
!     MODEL(I)%FORCING%AIR_DENSITY   = 'no'
!     MODEL(I)%FORCING%ICE_CONC      = 'no'
!     MODEL(I)%FORCING%ICE_PARAM1    = 'no'
!     MODEL(I)%FORCING%ICE_PARAM2    = 'no'
!     MODEL(I)%FORCING%ICE_PARAM3    = 'no'
!     MODEL(I)%FORCING%ICE_PARAM4    = 'no'
!     MODEL(I)%FORCING%ICE_PARAM5    = 'no'
!     MODEL(I)%FORCING%MUD_DENSITY   = 'no'
!     MODEL(I)%FORCING%MUD_THICKNESS = 'no'
!     MODEL(I)%FORCING%MUD_VISCOSITY = 'no'
!     MODEL(I)%ASSIM%MEAN            = 'no'
!     MODEL(I)%ASSIM%SPEC1d          = 'no'
!     MODEL(I)%ASSIM%SPEC2d          = 'no'
!     MODEL(I)%RESOURCE%RANK_ID      = I
!     MODEL(I)%RESOURCE%GROUP_ID     = 1
!     MODEL(I)%RESOURCE%COMM_FRAC    = 0.00,1.00
!     MODEL(I)%RESOURCE%BOUND_FLAG   = F
!
!     MODEL(4)%FORCING = 'no' 'no' 'no' 'no' 'no' 'no'
!
!     MODEL(2)%RESOURCE = 1 1 0.00 1.00 F
! -------------------------------------------------------------------- !
&MODEL_GRID_NML
  MODEL(1)%NAME                  = 'NS3'
  MODEL(1)%FORCING%WINDS         = 'NS3' ! Shall be the common 'NSB3' after tests
  MODEL(2)%NAME                  = 'BAL3'
  MODEL(2)%FORCING%WINDS         = 'BAL3'
  MODEL(2)%FORCING%ICE_CONC      = 'BAL3'
  MODEL(2)%FORCING%ICE_PARAM1    = 'BAL3'
  MODEL(2)%FORCING%ICE_PARAM2    = 'BAL3'
  MODEL(2)%FORCING%ICE_PARAM3    = 'BAL3'
  MODEL(2)%FORCING%ICE_PARAM4    = 'BAL3'
  MODEL(2)%FORCING%ICE_PARAM5    = 'no'
  MODEL(1)%RESOURCE = 3 3 0.00 0.50 T
  MODEL(2)%RESOURCE = 3 3 0.50 1.00 T
  MODEL(1)%RESOURCE%BOUND_FLAG   = T

/
! src/w3cspcmd.F90 src/w3gridmd.F90 w3initmd.F90 w3iogomd.F90 w3iopomd.F90
! w3servmd.F90 w3wavemd.F90
!
!
!

! -------------------------------------------------------------------- !
! Define the output types point parameters via OUTPUT_TYPE_NML namelist
!
! * index I must match indexes from 1 to DOMAIN%NRGRD
!
! * ALLTYPE will apply the output types for all the model grids
!
! * ITYPE(I) will apply the output types for the model grid number I
!
! * need DOMAIN%UNIPTS equal true to use a unified point output file
!
! * the point file is a space separated values per line :
!   longitude latitude 'name' (C*40)
!
! * the detailed list of field names is given in model/nml/ww3_shel.nml :
!  DPT CUR WND AST WLV ICE IBG TAU RHO D50 IC1 IC5
!  HS LM T02 T0M1 T01 FP DIR SPR DP HIG
!  EF TH1M STH1M TH2M STH2M WN
!  PHS PTP PLP PDIR PSPR PWS PDP PQP PPE PGW PSW PTM10 PT01 PT02 PEP TWS PNR
!  UST CHA CGE FAW TAW TWA WCC WCF WCH WCM FWS
!  SXY TWO BHD FOC TUS USS P2S USF P2L TWI FIC USP TOC
!  ABR UBR BED FBB TBB
!  MSS MSC WL02 AXT AYT AXY
!  DTD FC CFX CFD CFK
!  U1 U2
!
! * output track file formatted (T) or unformated (F)
!
! * namelist must be terminated with /
! * definitions & defaults:
!     ALLTYPE%FIELD%LIST         =  'unset'
!     ALLTYPE%POINT%NAME         =  'unset'
!     ALLTYPE%POINT%FILE         =  'points.list'
!     ALLTYPE%TRACK%FORMAT       =  T
!     ALLTYPE%PARTITION%X0       =  0
!     ALLTYPE%PARTITION%XN       =  0
!     ALLTYPE%PARTITION%NX       =  0
!     ALLTYPE%PARTITION%Y0       =  0
!     ALLTYPE%PARTITION%YN       =  0
!     ALLTYPE%PARTITION%NY       =  0
!     ALLTYPE%PARTITION%FORMAT   =  T
!
!     ITYPE(3)%TRACK%FORMAT      =  F
! -------------------------------------------------------------------- !
&OUTPUT_TYPE_NML
  ALLTYPE%FIELD%LIST     = 'HS DIR SPR WND T02'
  ITYPE(2)%FIELD%LIST    = 'HS DIR SPR WND ICE ICE1'

  ITYPE(1)%POINT%FILE    = 'points.NS3'
  ITYPE(2)%POINT%FILE    = 'points.BAL3'
/



! -------------------------------------------------------------------- !
! Define output dates via OUTPUT_DATE_NML namelist
!
! * index I must match indexes from 1 to DOMAIN%NRGRD
! * ALLDATE will apply the output dates for all the model grids
! * IDATE(I) will apply the output dates for the model grid number i
! * start and stop times are with format 'yyyymmdd hhmmss'
! * if time stride is equal '0', then output is disabled
! * time stride is given in seconds
! * it is possible to overwrite a global output date for a given grid
!
! * namelist must be terminated with /
! * definitions & defaults:
!     ALLDATE%FIELD%START         =  '19680606 000000'
!     ALLDATE%FIELD%STRIDE        =  '0'
!     ALLDATE%FIELD%STOP          =  '19680607 000000'
!     ALLDATE%POINT%START         =  '19680606 000000'
!     ALLDATE%POINT%STRIDE        =  '0'
!     ALLDATE%POINT%STOP          =  '19680607 000000'
!     ALLDATE%TRACK%START         =  '19680606 000000'
!     ALLDATE%TRACK%STRIDE        =  '0'
!     ALLDATE%TRACK%STOP          =  '19680607 000000'
!     ALLDATE%RESTART%START       =  '19680606 000000'
!     ALLDATE%RESTART%STRIDE      =  '0'
!     ALLDATE%RESTART%STOP        =  '19680607 000000'
!     ALLDATE%BOUNDARY%START      =  '19680606 000000'
!     ALLDATE%BOUNDARY%STRIDE     =  '0'
!     ALLDATE%BOUNDARY%STOP       =  '19680607 000000'
!     ALLDATE%PARTITION%START     =  '19680606 000000'
!     ALLDATE%PARTITION%STRIDE    =  '0'
!     ALLDATE%PARTITION%STOP      =  '19680607 000000'
!
!     ALLDATE%RESTART             =  '19680606 000000' '0' '19680607 000000'
!
!     IDATE(3)%PARTITION%START    =  '19680606 000000'
! -------------------------------------------------------------------- !
&OUTPUT_DATE_NML
  ALLDATE%FIELD%START         = '__SIMSTART__'
  ALLDATE%FIELD%STRIDE        = '1200'
  ALLDATE%FIELD%STOP          = '__SIMSTOP__'
  ALLDATE%POINT%START         =  '__SIMSTART__'
  ALLDATE%POINT%STRIDE        = '1200'
  ALLDATE%POINT%STOP          = '__SIMSTOP__'

  ALLDATE%RESTART             =  '__SIMSTART__' '__HOTDELTA__' '__SIMSTOP__'
/



! -------------------------------------------------------------------- !
! Define homogeneous input via HOMOG_COUNT_NML and HOMOG_INPUT_NML namelist
! 
! * the number of each homogeneous input is defined by HOMOG_COUNT
! * the total number of homogeneous input is automatically calculated
! * the homogeneous input must start from index 1 to N
! * if VALUE1 is equal 0, then the homogeneous input is desactivated
! * NAME can only be MOV
! * each homogeneous input is defined over a maximum of 3 values detailled below :
!     - MOV is defined by speed and direction
! 
! * namelist must be terminated with /
! * definitions & defaults:
!     HOMOG_COUNT%N_MOV             =  0
! 
!     HOMOG_INPUT(I)%NAME           =  'unset'
!     HOMOG_INPUT(I)%DATE           =  '19680606 000000'
!     HOMOG_INPUT(I)%VALUE1         =  0
!     HOMOG_INPUT(I)%VALUE2         =  0
!     HOMOG_INPUT(I)%VALUE3         =  0
! -------------------------------------------------------------------- !
! &HOMOG_COUNT_NML
!   HOMOG_COUNT%N_IC2                =  1
!   HOMOG_COUNT%N_IC3                =  1
!   HOMOG_COUNT%N_IC4                =  1
! /
! 
! &HOMOG_INPUT_NML
!   HOMOG_INPUT(1)%NAME        = 'IC2'
!   HOMOG_INPUT(1)%VALUE1      = 3.0
! 
!   HOMOG_INPUT(2)%NAME        = 'IC3'
!   HOMOG_INPUT(2)%VALUE1      = 917.0
!  
!   HOMOG_INPUT(3)%NAME        = 'IC4'
!   HOMOG_INPUT(3)%VALUE1      = 1.0
!   
! /
 
 
! -------------------------------------------------------------------- !
! WAVEWATCH III - end of namelist                                      !
! -------------------------------------------------------------------- !
