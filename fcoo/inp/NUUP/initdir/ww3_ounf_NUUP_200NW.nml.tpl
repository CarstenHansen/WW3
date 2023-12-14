! -------------------------------------------------------------------- !
! WAVEWATCH III - ww3_ounf.nml - Grid output post-processing           !
! -------------------------------------------------------------------- !


! -------------------------------------------------------------------- !
! Define the output fields to postprocess via FIELD_NML namelist
!
! * the detailed list of field names FIELD%LIST is given in ww3_shel.nml
!  DPT CUR WND AST WLV ICE TAU RHO IBG D50 IC1 IC5
!  HS LM T02 T0M1 T01 FP DIR SPR DP HIG WNM
!  EF TH1M STH1M TH2M STH2M WN
!  PHS PTP PLP PDIR PSPR PWS PDP PQP PPE PGW PSW PTM10 PT01 PT02 PEP TWS PNR
!  UST CHA CGE FAW TAW TWA WCC WCF WCH WCM FWS
!  SXY TWO BHD FOC TUS USS P2S USF P2L TWI FIC TOC XSP
!  ABR UBR BED FBB TBB
!  MSS MSC WL02 AXT AYT AXY
!  DTD FC CFX CFD CFK
!  U1 U2
!
! * namelist must be terminated with /
! * definitions & defaults:
!     FIELD%TIMESTART            = '19000101 000000'  ! Stop date for the output field
!     FIELD%TIMESTRIDE           = '0'                ! Time stride for the output field
!     FIELD%TIMECOUNT            = '1000000000'       ! Number of time steps
!     FIELD%TIMESPLIT            = 6                  ! [0(nodate),4(yearly),6(monthly),8(daily),10(hourly)]
!     FIELD%LIST                 = 'unset'            ! List of output fields
!     FIELD%PARTITION            = '0 1 2 3'          ! List of wave partitions ['0 1 2 3 4 5']
!     FIELD%SAMEFILE             = T                  ! All the variables in the same file
!     FIELD%VECTOR               = T                  ! Vector [T] or dir/magnitude [F] for
!                                                     !   directional fields
!     FIELD%TYPE                 = 3                  ! [2 = SHORT, 3 = it depends , 4 = REAL]
!     FIELD%FCVARS               = F                  ! Generate auxiliary forecast variables
!                                                     !   (forecast_period and forecast_reference_time)
!     FIELD%TIMEREF              = 'unset'            ! "Forecast reference time" for calculating
!                                                     !   forecast_period; defaults to TIMESTART
!     FIELD%TIMEVAR              = 'D'                ! Time var type ['D' = DOUBLE, 'I' = INT64]
!     FIELD%TIMEUNIT             = 'D'                ! Time units ['D' = days, 'I' = seconds]
!     FIELD%TIMEEPOCH            = '19900101 000000'  ! Epoch used for encoding of NC time variables
!     FIELD%NOVAL                = UNDEF              ! Value for wet cells that have an UNDEF value
!     FIELD%MAPSTA               = .TRUE.             ! Output MAPSTA field in file?
!
! Note: If FCVARS = T, the following auxiliary variables will be generated
! (see the manual entry for ww3_ounf for more information):
!
!    - forecast_reference_time: The time associated with the "analysis" of
!      the current forecast. Defaults to TIMESTART if TIMEREF not set.
!
!    - forecast_period: the time period elapsed (in seconds) since the
!      associated "forecast_reference_time".
! -------------------------------------------------------------------- !
&FIELD_NML
  FIELD%TIMESTART        =  '__SIMSTART__'
  FIELD%TIMESTRIDE       =  '1200'
  FIELD%TIMESPLIT        =  0
  FIELD%LIST             =  'DPT WND ICE HS T02 T0M1 DIR SPR TUS USS SVP MFIT DTD FC'
  FIELD%SAMEFILE         =  T
  FIELD%TYPE             =  4
  FIELD%FCVARS           =  F
/
! add 'IC1'

! -------------------------------------------------------------------- !
! Define the content of the output file via FILE_NML namelist
!
! * namelist must be terminated with /
! * definitions & defaults:
!     FILE%PREFIX        = 'ww3.'            ! Prefix for output file name
!     FILE%NETCDF        = 3                 ! Netcdf version [3|4]
!     FILE%IX0           = 1                 ! First X-axis or node index
!     FILE%IXN           = 1000000000        ! Last X-axis or node index
!     FILE%IY0           = 1                 ! First Y-axis index
!     FILE%IYN           = 1000000000        ! Last Y-axis index
! -------------------------------------------------------------------- !
&FILE_NML
  FILE%PREFIX        = 'NUUP_200NW.'
  FILE%NETCDF        = 3
/

! -------------------------------------------------------------------- !
! Define the content of the output file via SMC_NML namelist
!
! * For SMC grids, IX0, IXN, IY0 and IYN from FILE_NML are not used.
!   Two types of output are available:
! *   TYPE=1: Flat 1D "seapoint" array of grid cells.
! *   TYPE=2: Re-gridded regular grid with cell sizes being an integer
! *           multiple of the smallest SMC grid cells size.
!
! * Note that the first/last longitudes and latitudes will be adjusted
!  to snap to the underlying SMC grid edges. CELFAC is only used for
!  type 2 output and defines the output cell sizes as an integer
!  multiple of the smallest SMC Grid cell size. CELFAC should be a
!  power of 2, e.g: 1,2,4,8,16, etc...
!
! * namelist must be terminated with /
! * definitions & defaults:
!     SMC%TYPE          = 1              ! SMC Grid type (1 or 2)
!     SMC%SXO           = -999.9         ! First longitude
!     SMC%EXO           = -999.9         ! Last longitude
!     SMC%SYO           = -999.9         ! First latitude
!     SMC%EYO           = -999.9         ! Last latitude
!     SMC%CELFAC        = 1              ! Cell size factor (SMCTYPE=2 only)
!     SMC%NOVAL         = UNDEF          ! Fill value for wet cells with no data
! -------------------------------------------------------------------- !
&SMC_NML
/


! -------------------------------------------------------------------- !
! WAVEWATCH III - end of namelist                                      !
! -------------------------------------------------------------------- !
