!> @file w3stvpmd.F90
!> @brief Extended spectral tail for Stokes drift profile.
!> @author Carsten Hansen @date 21-Jan-2021
#include "w3macros.h"
!/ ------------------------------------------------------------------- /
      MODULE W3STVPMD
!/ over
!/                  +------------------------------------------+
!/                  | Stokes drift  profile with extended      |
!/                  | spectral tail                            |
!/                  |          Carsten Hansen                  |
!/                  | Danish defence Joint GEOMETOC Support    |
!/                  |                            FORTRAN 90    |
!/                  | Last update :             03-Jan-2023    |
!/                  +------------------------------------------+
!/
!/    02-Jan-2023 : Implement option USXT = 'Pf-5': Frequency^{-5} tail with no
!/                  further broadening as function of frequency
!-------------------
!>
!>  @brief Purpose, and output data.
!>
!>  @details
!>  Calculate the vertical profile of the Stokes drift, using an assumed
!>  power-law spectral tail (Phillips' frequency^{-5} for the diagnostic parts
!>  of the spectrum at high frequencies. It assumes that there is no spectral
!>  broadening as function of frequency, i.e. the so-called 'first circular
!>  moment' (Ewans, 1997) is assumed constant in the diagnostic tail. A namelist
!>  parameter USXT = 'Pf-5' must be set to add this tail.
!>  
!>  This replaces earlier versions with a described directional distribution
!>  following the empirical results of Ewans (1997) in conjunction with an
!>  assumed power-law of frequency^{-4}. This could be re-implemented and
!>  activated by a namelist parameter USXT = 'DoEw'. However the numerical
!>  logics are somewhat more complicated to comprehend, and the difference is
!>  small.
!>
!>  All calculations are performed in parallel and the Stokes drift is calculated
!>  at NDEPTH discrete depths given in an array Z_S(1:NDEPTH).
!>  The sum of prognostic and diagonostic (tail) contributions are stored in
!>  arrays U_S(1:NDEPTH), V_S(1:NDEPTH).
!>
!>  @brief Usage.
!>
!>  @details
!>  In ww3_grid.inp or namelists.nml for ww3_grid.nml, you must define three
!>  namelist parameters for the Stokes profile (&STVP), and two for the extended
!>  parametric spectral tail (&XSTP):
!>  @verbatim
!>  &STVP NDP = 11, DSC = 2.0, BP = 2.0 /
!>  &XSTP USXT = 'Pf-5' , USXF = 2.0 /
!>  @endverbatim
!>  - NDP: Number of depths (->NDEPTH)
!>  - DSC: Depth scale specifying the largest depth Z(NDP)=DSC/(2pi/T02)^2*GRAV
!>  - BP: Power of profile depth progression        
!>  - USXT: Tail type (C*4); 'Pf-5': Frequency^{-5} tail - no further broadening
!>  -                        'none': Prognostic spectrum frequencies, only
!>  - USXF: Tail truncation frequency (low-pass). Default 10 Hz (~infinity)
!>
!>  @brief Estimation of Stokes drift of the diagnostic tail.
!>
!>  @details
!>  The Stokes drift U_S is determined by the wave spectrum H(sig,theta) over
!>  (intrinsic) frequencies sig (= 2 pi f in the absense of currents) and
!>  directions theta. For wind waves over deep water the contribution
!>  in a frequency bin dsig to the directionally integrated Stokes drift at
!>  depth z below the surface is (as formulated by Ewans, 1997),
!> 
!>  @verbatim
!>    d U_S(sig,z) = 2/g sig^3 m1(sig) exp(-2kz) S(sig) d sig
!>  @endverbatim
!>  where k = sig^2/g. The directionally integrated spectrum is
!>  @verbatim
!>    S(sig) = int_-pi^pi H(sig,theta) d theta
!>  @endverbatim
!>  where H(sig,theta) is the 2D spectrum, and and m1(sig) is the so-called
!>  'first circular moment',
!>  @verbatim
!>    m1(sig) = int_-pi^pi cos(theta) H(sig,theta) d theta / S(sig)
!>  @endverbatim
!>
!>  With the option USXT = 'Pf-5' we assume that m1(sig) is constant for all
!>  sig > SIG(NK), and that the spectral tail has the form of Phillips (1985),
!>  @verbatim
!>    S(sig) = B g^2 sig^-5,
!>  @endverbatim
!>  Here, the factor B may vary with the integral wave properties, for example
!>  a dependance on wave age (U10/cP)*cos(theta_mean-theta_U). In the present
!>  usage, B is estimated based on a few of the highest frequency bins of
!>  the prognostic spectrum.
!>
!>  @brief Integrated pseudo-momentum over the diagnostic spectral tail.
!>
!>  @details
!>  The integrated pseudo-momentum [ per m^2 ] of the spectral tail is a vector
!>  rho g M_tail, 
!>  @verbatim
!>    M_tail = \int_ft^\fc int_-pi^pi kv/sig H(sig,theta) d theta d sig
!>  @endverbatim
!>  Here, ft it the tail start frequency, fc is a high-frequency cut-off, and
!>  kv is the spectral wave number vector, kv = k * (cos(theta), sin(theta)).
!>  Assume the tail components are all deep-water waves. Then k = sig^2/g.
!> 
!>  Assume a constant directional distribution,
!>  @verbatim
!>    H(sig,theta) = S(sig) H0(theta)
!>  @endverbatim
!>  with int_-pi^pi H0(theta) d theta = 1 and a constant mean direction theta0.
!>  The assumption implies that we have a constant value of m1(sig) = m1s, 
!>  @verbatim
!>    m1s = int_-pi^pi cos(theta-theta0) H0(theta) d theta
!>  @endverbatim
!>  
!>  For z = 0 at the very surface we find the contribution to the Stokes drift
!>  magnitude from the spectral tail between the frequencies ft and fc,
!>  @verbatim
!>    U_S_tail(z=0) = \int_ft^fc 2 m1s/g sig^3 B g^2 sig^-5 d sig
!>       = 2 m1s B g ( sig(ft)^-1 - sig(fc)^-1 )
!>        [ units: m s^-1 ]
!>  @endverbatim
!>  The tail contribution to the magnitude of the integral Stokes drift is
!>  @verbatim
!>    M_tail = \int_ft^fc m1s/g sig B g^2 sig^-5 d sig
!>       = 1/3 m1s B g ( sig(ft)^-3 - sig(fc)^-3 )
!>        [ units: m s ]
!>  @endverbatim
!>  At depth z>0, the Stokes drift magnitude is an 'incomplete-Gamma' integral
!>  @verbatim
!>    U_S_tail(z) = \int_ft^fc 2 m1s B g sig^-2 exp(-2 k z) d sig
!>       = \int_ft^fc m1s B g^1/2 k^-3/2 exp(-2 k z) d k
!>       = m1s B (2 g z)^1/2 \
!>         * ( GammaInc(-1/2, 2 k_t z) - GammaInc(-1/2, 2 k_c z) )
!>  @endverbatim
!>  where k_t = (2 pi ft)^2/g, k_c = (2 pi fc)^2/g .
!>
!>  In the actual numerical approximation we don't need to calculate the
!>  incomplete Gamma function. We can perform the discrete integration as:
!>  @verbatim
!>    U_S_tail(k_s z) = \sum_ft^fc m1s B g^1/2 k^-3/2 exp(-2 k z) D k
!>       =  \sum_ft^fc m1s B (g/k_s)^1/2 (k/k_s)^-3/2 exp(-2 z k_s k/k_s) D k/k_s
!>  @endverbatim
!>  This shows that we can apply a look-up function
!>  @verbatim
!>    NSIZ(IZK, IKt) = \sum_IKt^NP (k/k_s)^-3/2 exp(-2 z k_s k/k_s) D k/k_s
!>  @endverbatim
!>  where IZK is the index (1..NDP) of the user-defined discrete values of z k_s.
!>  The summation is over binned values of k/k_s over a subsection of the range
!>  ( XKH, XKH^2, ..., XKH^NP ) starting at bin number IKt where
!>  k_t/k_s ~ XKH^IKt. For the value of XKH we have chosen the same increments
!>  as for the model discrete wavenumber intervals, XKH = XFR**2, and we have
!>  set NP = NK. Similarly we subtract the highest frequencies beyond the
!>  cut-off k_c/k_s so that approximately,
!>  @verbatim
!>    U_S_tail(k_s z) = m1s B (g/k_s)^1/2 ( NSIZ(IZK, IKt) - NSIZ(IZK, IKc) )
!>  @endverbatim
!>  In the present implementation we *define* the wavenumber scale as
!>  k_s = SIG_S^2/g, where SIG_S = 2 pi/T02. This setting may be changed
!>  in the future, thus k_s is provided in the output as a netCDF variable
!>  'ksc'.
!>             
!/ ------------------------------------------------------------------- /
      USE CONSTANTS, ONLY: GRAV, TPI
!/
      USE W3ADATMD, ONLY: CG, WN, DW

! Integral pseudo-momentum vector and wave number scale K_S
      USE W3ADATMD, ONLY: M_X, M_Y, K_S
! Full Stokes profile vectors
      USE W3ADATMD, ONLY: U_S, V_S, ZK_S

      USE W3GDATMD, ONLY: SPND, SPDS, SPBP, NSEAL, DMIN, XFR, DTH,  &
                          SIG, DDEN, ECOS, ESIN, NK, NTH, USXT, USXF

      USE W3ODATMD, ONLY: NDST, NDSO, NDSE, IAPROC, NAPROC, NAPOUT, VERBOSENESS
      USE W3DISPMD, ONLY: DSIE, N1MAX, ECG1, EWN1


      PRIVATE

!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/

! Stokes drift for the diagnostic tail. Scalable array and depths factor
! DFH: Relative increment of wavenumber bin-to bin for the diagnostic Stokes
!      array
! NP:  Number of wavenumber bins in the diagnostic stokes array.

      REAL               :: XKH
      INTEGER            :: NP

! File id. for verbose output
      INTEGER            :: NDSV
!
      INTEGER            :: stvp_verbose = 1
      
! Look-up lists for Stokes spectral tail calculation
      REAL, allocatable  :: IKrange(:) ! (NP-1)
      REAL, allocatable  :: NSIZ(:,:) ! (NZ,NP)

! Power of the spectral tail decay
      REAL               :: Pn = 5.0

! k_t: Wavenumber at bin IKT; k_c: wavenumber at the cut-off frequency USXF
      REAL               :: k_t, k_c

! Profile arrays
      REAL, allocatable  :: U_S_tail(:)

      REAL, allocatable  :: A2S(:), TZ_S(:), ETKZ(:)

! Normalized integral momentum (might be an array if transition freq IKT varies)
      REAL               :: Im

! Tail transition indices
! IKST: Lower edge of tail estimation from prognostic spectrum
! IKT:  Upper edge of tail estimation.
      INTEGER            :: IKST, IKT

! For calculation of shallow water Stokes drift
      REAL, allocatable  :: fkd(:,:), fkdm(:,:) ! (1:NK,1:NSEAL)

! 2* (WN(IK,JSEA)-WN(IK-1,JSEA))
      REAL, allocatable  :: DWN(:,:)

      PUBLIC :: CALC_STVP
!/
      CONTAINS


!/ ------------------------------------------------------------------- /
!> @brief STVP_INIT(): Initialisation for CALC_STVP ().
!>        
!> @details        
!>  Common parameters and arrays are set based on the namelist settings.
!>  Lookup tables are calculated for the diagnostic tail.
!>
!/ ------------------------------------------------------------------- /
      SUBROUTINE STVP_INIT ()

!/            +------------------------------------------+
!/            | Extended spectral tail for Stokes drift  |
!/            |          Carsten Hansen                  |
!/            | Danish defence Joint GEOMETOC Support    |
!/            | Last update :             22-Dec-2017    |
!/            +------------------------------------------+
!/
!/    02-May-2011 : Origination. Development program w3stokes.ftn converted
!/                  from tested Python code
!/    27-Jun-2011 : First module version w3stvpmd.ftn
!/    22-Dec-2017 : Adoption to version WW3 v5.16
!/    15-Aug-2019 : Adoption to version WW3 v6.07
!/                  Called in first call CALC_XINIT()
!
      USE W3SERVMD, ONLY : EXTCDE

      IMPLICIT NONE

!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
! Location in spatial/spectral grid
      INTEGER             :: JSEA, ISEA, IK
! Control integers
      INTEGER             :: IERR

! Common Stokes drift depth distribution parameters
      REAL                :: DEPTH, XZK, kd, kd_t, kd_d, fkd2, fkdm2, e2kd, shkd4
      INTEGER             :: IZ
      logical             :: OPENED
! Frequency at the cut-off frequency USXF
      REAL                :: sig_c
! Building the normalized stokes tail look-up
      INTEGER             :: nt
      REAL, allocatable   :: DNSI0(:)
      REAL                :: tRk, XKHpm3h, XKHpm3nh, XKHpm1, XKHpp1

      
      stvp_verbose = VERBOSENESS%STVP
      NDSV = NDST
      ! Highly verbose output only with test output
      IF ( stvp_verbose .gt. 0 ) THEN
        INQUIRE (NDSV,OPENED=OPENED)
        IF ( .NOT. OPENED ) THEN
          stvp_verbose = 0
          IF ( IAPROC == NAPOUT) THEN
            NDSV = NDSO
            stvp_verbose = VERBOSENESS%STVP
            END IF
          END IF
        END IF

!/ ------------------------------------------------------------- /
! Allocate arrays for depths, prognostic stokes and tail stokes
!
      IF ( stvp_verbose .gt. 0 ) &
        WRITE (NDSV, 904),  'INITIALIZE Stokes drift for ', SPND, 'depths'

      ! Shallow water Stokes parameters
      allocate ( fkd(NK,NSEAL), stat=IERR )
      allocate ( fkdm(NK,NSEAL), stat=IERR )

      allocate( ETKZ(SPND), stat=IERR )

      ALLOCATE( A2S(SPND), stat=IERR )

      allocate( TZ_S(SPND), stat=IERR )

      allocate( U_S_tail(SPND), stat=IERR )

      ! For calculating exp(2 K Z) = exp(2 K0 Z) * exp( 2 (K-K0) Z )
      allocate( DWN(NK,NSEAL), stat=IERR )
      

      ! Construct a scalable array ZK_S of the stokes profile for each
      ! individual sea point:
      !   ZK_S  = K_S * Z(m),
      ! where K_S is a wave number scale and Z are the depths.

      ! Configured values:
      ! SPND: Number of depths for the Stokes profile
      ! SPDS: Depth scale so that the deepest point of the profile is
      ! SPBP: Power of profile depth progression
      !   Z_S(SPND) = SPDS / K_S
      ! Let K_S = (TPI/Tm02)**2 / GRAV
      ! Dimensionless depths will be constructed as
      !      ZK_S(IZ) = SPDS*(XZK**((IZ-1.)**SPBP) - 1.) for IZ = 1 .. SPND
      ! The exponential base XZK is defined so that the deepest point is
      ! ZK_S(SPND) = SPDS,

      XZK = 2.**((SPND*1.-1.)**(-SPBP))
      DO IZ = 1,SPND
        ZK_S(IZ) = XZK**((IZ*1.-1.)**SPBP)
      END DO
      ZK_S(:) = SPDS*(ZK_S(:) - 1.)

      ! Shallow water Stokes profile parameters fkd(1:NK,1:NSEAL),
      ! fkdm(1:NK,1:NSEAL) for calculation of the action-to-stokes factor
      ! A2S(:) in subroutine STOKES_PROFILE

      ! In the deep water limit,
      fkd = 1.
      fkdm = 0.
      
      ! Apply a linear interpolation in the range kd_t < kd <= kd_d of
      ! transitional to deep water. Our choise is:
      kd_t = 2. 
      kd_d = 3.
      
      e2kd = exp (2.*kd_t)
      shkd4 = 4. * (sinh( kd_t ))**2
      ! For use in calculation of fdk, fdkm:
      fkd2 = e2kd/shkd4
      fkdm2 = 1./ (e2kd * shkd4)

      DO JSEA = 1,NSEAL
#ifdef W3_DIST
        ISEA       = IAPROC + (JSEA-1)*NAPROC
#endif
#ifdef W3_SHRD
        ISEA       = JSEA
#endif

        ! For finite depth Stokes profile shape, parameters fkd, fkdm:
        DEPTH  = MAX ( DMIN, DW(ISEA) )

        DO IK=1,NK
          kd = WN(IK, ISEA) * DEPTH

          ! kd>2: pure deep water (then fkd(IK,JSEA)=1 and fkdm(IK,JSEA)=0)
          IF (kd .gt. kd_d) exit

          ! kd_t < kd <= kd_d: Linear interpolation of transitional to deep water
          IF (kd .gt. kd_t) THEN
            fkd(IK,JSEA) = fkd2 * (kd_d - kd) + (kd - kd_t)
            fkdm(IK,JSEA) = fkdm2 * (kd_d - kd)
            CYCLE
            END IF

          ! kd<kd_t: shallow water
          e2kd = exp (2.*kd)
          shkd4 = 4. * (sinh( kd ))**2
          fkd(IK,JSEA) = e2kd / shkd4
          fkdm(IK,JSEA) = 1./ (e2kd * shkd4)
          ! In subroutine STOKES_PROFILE the action-to-stokes factor is
          ! calculated as
          ! A2Sf[IK,iz] = fkd[IK,JSEA] * ETKZ + fkdm[IK,JSEA] / ETKZ, where
          ! ETKZ = EXP( - 2. * kz[iz] )
          ! Note: At the surface where kz[1] = 0 we have ETKZ = 1, and 
          !       A2Sf[IK,1] = 0.5 * COSH(2.*kd) / SINH(kd)**2
          !                  = fkd[IK,JSEA] + fkdm[IK,JSEA]
          
        END DO ! IK=1,NK

        ! For calculating exp(-2 K Z) = exp(-2 K0 Z) * exp( -2 (K-K0) Z ),
        ! DWN(IK) = K - K0 = WN(IK) - WN(IK-1)
        DWN(1,JSEA) = WN(1,ISEA)
        DWN(2:NK,JSEA) = WN(2:NK,ISEA) - WN(1:NK-1,ISEA)

      END DO ! JSEA = 1,NSEAL

      ! IKT: The upper edge of the prognostic range.
      ! For a truncated spectrum with no tail extension:
      IKT = NK
      IKST = NK+1
      NP = 0
      
      IF ( USXF*TPI <= SIG(NK) ) return
      
      
!/ ------------------------------------------------------------- /
! Initializations for calculations related to the diagnostic tail range

! cha 20121123: Shift fitting interval limits away from NK, from NK-4 to NK-1
      IF ( USXT == 'Pf-5' ) THEN
         ! Phillips saturation
         Pn = 5
         IKT = NK - 1
! cha 20110627: IKST: Lower edge of tail fitting to the prognostic spectrum
         IKST = NK - 4
      ELSE IF ( USXT == 'none' ) THEN
         RETURN
      ELSE
        WRITE (NDSV, 916), 'STVP_INIT: Unknown spectral tail type: ', USXT
        CALL EXTCDE ( 73 )
      END IF

! NP, XKH: Number of bins and relative wavenumber increment of the lookup
! array for the diagnostic part of the Stokes profile
      NP = NK
      XKH = XFR**2

      allocate ( DNSI0(NP), stat=IERR )
      allocate( NSIZ(SPND, NP), stat=IERR )
      allocate( IKrange(NP-1), stat=IERR )

      ! Wavenumber at bin IKT, approximated as deep water waves
      k_t = SIG(IKT)**2/GRAV
      
      ! The integrated pseudo-momentum of the diagnostic spectral tail is
      ! M_tail = m1Bg * Im        
      Im = 1./(3. * SIG(IKT)**3)

      if ( USXF < 10.0 ) then
        sig_c = TPI * USXF
        k_c = sig_c**2/GRAV
        Im = Im - 1./(3. * sig_c**3)
      end if
      
      ! Aassuming deep water waves overall, (g/k_s)^1/2 = g/SIG_S
      ! The tail Stokes drift profile is
      ! U_S_tail(Z) = m1s B g/SIG_S NSI(IZK, IKt)
      ! where
      ! NSI(IZK, IKt)= \sum_IKt^NP (k/k_s)^-3/2 exp (-2k/k_s ZK_S) D k/k_s
      
      ! k/k_s = k_t/k_s * ( XKH, XKH^2, ..., XKH^N, ... ), N=1..NP
      ! d k/k_s = k_t/k_s * (...,( XKH^(N+1) - XKH^(N-1) )/2,... ), N=1..NP

      ! DNSI0(1,2,...) = (k/k_s)^-3/2 D k/k_s
      !   = (...,XKH^-3N/2 * ( XKH^(N+1) - XKH^(N-1) )/2,... ), N=1..NP

      XKHpm3h = XKH**(-1.5)
      XKHpm3nh = 1
      XKHpm1 = 1/XKH * 0.5
      XKHpp1 = XKH * 0.5
      DO nt = 1,NP
        ! DNSI0(nt) = XKH**-(3*nt/2) * 0.5 * (XKH**(nt+1) - XKH**(nt-1))
        ! XKHp3nmh = (XKH**(-1.5))**nt
        ! XKHpp1 = XKH**(nt+1)
        ! XKHpm1 = XKH**(nt-1)
        XKHpm3nh = XKHpm3nh * XKHpm3h
        XKHpp1 = XKHpp1 * XKH
        XKHpm1 = XKHpm1 * XKH
        DNSI0(nt) = XKHpm3nh  * (XKHpp1 - XKHpm1)      
      END DO      
      
      ! NSIZ(KZ,nt): A 'normalized, look-up' vertical Stokes Drift profile
      ! shape for spectral bins of the diagnostic tail.
      ! Also, to help determine the look-up index nt, set IKrange(nt) = XKH**nt
      
      tRk = 2.
      DO nt=1,NP
        tRk = tRk*XKH ! = 2 * XKH**nt
        NSIZ(:,nt) = DNSI0(nt) * exp( -tRk * ZK_S(:) )
        if (nt == NP) exit
        IKrange(nt) = tRk
      END DO
      IKrange(:) = 0.5 * IKrange(:)
      ! Accumulated sum downwards from the high end of the spectral tail
      DO nt=NP-1,1,-1
         NSIZ(:,nt) = NSIZ(:,nt) + NSIZ(:,nt+1) 
      END DO

      deallocate ( DNSI0 )
      
      ! Application:
      ! ! 1) Determine nt so that k_t/k_s = XKH**nt, and nc so k_c/k_s = XKH**nc
      ! nt = index(IKrange,k_t/k_s); nc = index(IKrange*k_c,k_c/k_s)
      ! ! 2) Diagnostic tail contribution to the drift profile
      ! SIG_S = TPI/T02
      ! U_S_tail(ft)(:) = m1Bg / SIG_S * ( NSIZ(:,nt) - NSIZ(:,nc) )
      ! ! 3) Total Stokes drift
      ! M_tail = m1Bg * Im
  904 FORMAT ('  ', A,I5,A)
  916 FORMAT ('  ', A, A)

      END SUBROUTINE STVP_INIT


      SUBROUTINE CALC_STVP (A)
!/ ------------------------------------------------------------------- /
!> @brief CALC_STVP (A) - Calculate the profiles looping over all seapoints
!>      
!> @details
!>  At each WW3 output time step, the subroutine STOKES_PROFILE is called to
!>  calculate the Stokes drift profile and the pseudo-momentum for all sea points
!>  JSEA=1,NSEAL. The result is stored in a common array USVP(JSEA,1:3+2*SPND).
!>   
!/ ------------------------------------------------------------------- /

!/
!/            +------------------------------------------+
!/            | Calculate all the Stokes drift profiles  |
!/            |          Carsten Hansen                  |
!/            | Danish defence Joint GEOMETOC Support    |
!/            |                        FORTRAN 90        |
!/            | Last update :              02-Jan-2023   |
!/            +------------------------------------------+
!/
!/
!/    10-Oct-2011: Interface routine for Stokes drift profile calculation
!/                 To be called from W3OUTG() *after* calculation of both
!/                 integral parameters and swell partitions (W3CPRT)
!/                 The code structure is modelled after w3iosfmd.ftn
!/    22-Dec-2017 : Adoption to version WW3 v5.16
!/    02-Jan-2023 : Append the tail directly after bin IKT


!/    Remarks: This subroutine must be called after T02 has been calculated
!/             This subroutine can only be called under switch STVP.

! USVP: Full Stokes profile variables + M_X, M_Y, K_S, but not ZK_S
      USE W3ADATMD, ONLY: T02, USVP
      USE W3GDATMD, ONLY: NSEA, MAPSF, MAPSTA ! module scope:,NSEAL,NK,NTH,SIG
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(IN)        :: A(NTH,NK,0:NSEAL)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: ISEA, JSEA, IX, IY
      INTEGER, SAVE           :: OSTEP = 0


      IF ( .NOT. ALLOCATED(DWN) ) CALL STVP_INIT()
      IF ( stvp_verbose .GT. 0 ) THEN
           WRITE (NDSV, 912), OSTEP, 'Calculate Stokes drift ..'
           OSTEP = OSTEP + 1
      END IF

!
!
! -------------------------------------------------------------------- /
! 1.  Loop over sea points
!
      DO JSEA=1, NSEAL
!
! -------------------------------------------------------------------- /
! 2.  Check need for processing
!
#ifdef W3_DIST
        ISEA       = IAPROC + (JSEA-1)*NAPROC
#endif
#ifdef W3_SHRD
        ISEA       = JSEA
#endif
        IX     = MAPSF(ISEA,1)
        IY     = MAPSF(ISEA,2)
!
        IF ( MAPSTA(IY,IX) .LE. 0 ) CYCLE
!
! -------------------------------------------------------------------- /
! 3.  Perform Stokes drift profile calculation
        CALL STOKES_PROFILE ( A(:,:,JSEA), JSEA, T02(JSEA) )
!
! -------------------------------------------------------------------- /
! 4.  Save in output array
!
        USVP(JSEA,1) = K_S
        USVP(JSEA,2) = M_X
        USVP(JSEA,3) = M_Y
        USVP(JSEA,4:3+SPND) = U_S(1:SPND)
        USVP(JSEA,4+SPND:3+2*SPND) = V_S(1:SPND)
!
! -------------------------------------------------------------------- /
! 5.  End of loop over sea points
!
      END DO

  912 FORMAT (2X,I3,1X,A)

      END SUBROUTINE CALC_STVP

!/
!/ ------------------------------------------------------------------- /
!/

!/ ------------------------------------------------------------------- /
!> @brief STOKES_PROFILE (A, JSEA, T02) - Calculate the profile and pseudo-momentum
!>      
!> @details
!>  For the individual sea point, Stokes drift is calculated and
!>  provided in arrays
!> @verbatim
!>                 U_S(1:SPND), V_S(1:SPND)
!> @endverbatim
!>  for the vertical profile, and the pseudo-momentum is written to the variables
!> @verbatim
!>                 M_X, M_Y
!> @endverbatim
!>  An individual depth scale 1/K_S is applied based on an integral wave period
!>  measure (T02) for the time step, and the actual depths Z(1:SPND) can be
!>  determined from the global dimensionless dimension ZK_S(1:SPND) as
!> @verbatim
!>         Z(1:SPND) =  1/K_S * ZK_S(1:SPND)
!> @endverbatim 
!>  The Stokes profile is a sum of a prognostic part and a tail part
!>  U_S_tail(1:SPND)
!/ ------------------------------------------------------------------- /
      SUBROUTINE STOKES_PROFILE(A, JSEA, T02)
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/

      REAL, INTENT(IN)            :: A(NTH,NK)
      INTEGER, INTENT(IN)         :: JSEA
      REAL, INTENT(IN)            :: T02
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      REAL                  :: SIG_S, DEPTH, KDPT
      INTEGER               :: numDepths, IZ

      ! Factors for the prognostic calculation
      REAL                  :: A2S0, AX, AY, EDWTZ, DDWTZ

      INTEGER               :: ISEA, ITH, IK, NIK

! Tail spectral level
      REAL                  :: m1Bg
      INTEGER               :: nt, nc
      
! Integral pseudo-momentum parameters
      REAL                  :: A2M, M_tail

! Stokes drift wave number scale and parameters for 'inlined version of WAVNU1'
      REAL                  :: SIX, R1
      INTEGER               :: I1

      REAL                  :: CTH, STH, ksi

#ifdef W3_DIST
      ISEA       = IAPROC + (JSEA-1)*NAPROC
#endif
#ifdef W3_SHRD
      ISEA       = JSEA
#endif

      ! SIG_S is an integral wave frequency chosen as SIG_S == 2pi/T02
      SIG_S = TPI/T02

      CTH=0.
      STH=0.
      M_tail = 0.
      U_S_tail = 0.

      DEPTH  = MAX ( DMIN, DW(ISEA) )

      ! Inlined version of WAVNU1.
      SIX    = SIG_S * SQRT(DEPTH)
      I1     = INT(SIX/DSIE)
      IF (I1 .LE. N1MAX) THEN
        R1 = SIX/DSIE - REAL(I1)
        K_S = ( (1. - R1)*EWN1(I1) + R1*EWN1(I1 + 1) ) / DEPTH
      ELSE
        K_S = SIG_S**2 / GRAV
      END IF

      ! From here, we will use the dimensionless water depth
      KDPT = K_S * DEPTH

      ! numDepths: Number of depths of the discrete profile above the sea floor
      numDepths = SPND ! SPND is the configured value

      DO IZ = SPND,1,-1
        IF ( ZK_S(IZ) > KDPT ) CYCLE
        numDepths = IZ
        EXIT
      END DO

      IF ( NP == 0 ) THEN

        IF ( stvp_verbose .gt. 0  .and. JSEA .eq. 1 ) &
          WRITE (NDSV, 912), 'Truncate at NK-1. No diagnostic spectral extension'

      ELSE

        IF ( stvp_verbose .gt. 1  .and. JSEA .eq. 1 ) &
             WRITE (NDSV, 912), '    .. for the diagnostic part'

        ! At the prognostic range cut-off frequency, estimate the product
        ! m1Bg of the 1st circular moment m1t, and the spectral level Bg,
        ! and also estimate the mean *vector* orientation CTH, STH:
        ! (cos(theta0),sin(theta0)) m1 B g = (CTH, STH)*m1Bg

        call PROG_EDGE( A, Pn, m1Bg, CTH, STH, JSEA)

        ! The integrated pseudo-momentum of the diagnostic tail has the magnitude
        M_tail = m1Bg * Im        

        ! Determine the index nt of the Stokes look-up, so that k_t/k_s = XKH**nt
        ! ( IKrange = XKH**nt / k_t for nt=1,..shape(NSIZ)(2)-1 )
            ! This would yield same result: nt = minloc( abs(IKrange - k_t/k_s) )
        ksi=k_t/K_S
        do nt=2,size(IKrange)
          if (IKrange(nt) > ksi) exit
        end do
        ! Note, we append this at the prognostic bin IKT which has full bin width

        ! The diagnostic Stokes drift profile, at dimensionless depths Z K_S,
        ! has the magnitude
        U_S_tail(:) = m1Bg / SIG_S * NSIZ(:,nt)
        
        if ( USXF < 10.0 ) then
          ! If we assume a high-frequency cut-off
          ! Determine the index nc of the cut-off, so that k_c/k_s = XKH**nt
          ksi=k_c/K_S
          do nc=nt,size(IKrange)
            if (IKrange(nc) > ksi) exit
          end do
        
          U_S_tail(:) = U_S_tail(:) - m1Bg / SIG_S * NSIZ(:,nc)
        end if
        
        ! The Stokes drift values are zero below the sea floor
        U_S_tail(numDepths+1:SPND) = 0.

      END IF ! ( NP == 0 )

      ! Prognostic spectrum range IK=1, IKT

      IF ( stvp_verbose .gt. 1 ) THEN
        WRITE (NDSV, 912), '    .. for the prognostic part'
      END IF

      ! Stokes profile exponential shape for IK=1:
      ETKZ(:) = 1

      ! Initialize summation over the prognostic bins
      U_S = 0.
      V_S = 0.
      M_X = 0.
      M_Y = 0.

      TZ_S(:) = -2.0 * ZK_S(:)/K_S

      DO IK=1, IKT

        ! Conversion factor for Action spectrum to Stokes at the surface, A2S0
        ! SIG: angular frequency
        ! CG: Group velocity
        ! DDEN(IK): Spectrum 2D bin size DTH * DSII(IK)
        ! WN(IK,ISEA): Wave number of the spectral bin IK
        ! (Same calculation as in CALC_U3STOKES at Z=0 in w3iogomd)

        ! Action-to-momentum factor, times angular frequency.
        ! M = A k, also in shallow water.
        A2M = WN(IK,ISEA) * DDEN(IK) / CG(IK,ISEA)
        ! Action-to-Stokes factor at surface
        A2S0 = 2.0 * A2M * SIG(IK)
        ! Action-to-pseudo-momentum factor
        A2M = A2M / SIG(IK)

        ! Stokes profile exponential shape in deep water at ISEA for IK:
        ! ETKZ(:) = EXP( - 2. * WN * Z_S(:) )
        ! Remember, the values of Z_S depend on the actual value of K_S
        ! ETKZ(1) == 1. at the surface, always

        ! In order to calculate the exponential of a small argument (< order 1),
        ! it may be faster to increment from IK-1: Use DWN(IK)

        ETKZ(:) = ETKZ(:) * EXP( DWN(IK,JSEA) * TZ_S(:) )

        ! Action spectrum at frequency bin IK projected on x and y
        AX = 0.
        AY = 0.
        DO ITH=1, NTH
          AX = AX + A(ITH,IK) * ECOS(ITH)
          AY = AY + A(ITH,IK) * ESIN(ITH)
        END DO

        ! Stokes vertical profile for a finite water depth

        ! Action-to-Stokes-profile factor A2S
        IF (fkdm(IK,JSEA) .le. 0.) THEN
          ! Deep water
          A2S(:) = A2S0 * ETKZ(:)
        ELSE
          ! Finite water depth: Gridded parameters fkd, fkdm have been set in
          ! STVP_INIT with a linear interpolation in a transitional
          ! range, e.g. 2. < kd <= 3., from an exact shallow expression for
          !  kd <= 2. to assumed deep water for kd > 3. where fkd=1, fkdm=0 
          A2S(:) = A2S0 * ( fkd(IK,JSEA) * ETKZ(:) + fkdm(IK,JSEA) / ETKZ(:) )
        END IF

        ! Stokes drift projected on grid-Eastern and grid-Northern components
        U_S(:) = U_S(:) + AX * A2S(:)
        V_S(:) = V_S(:) + AY * A2S(:)

        ! Integral pseudo-momentum per unit surface area, divided by rho_w*GRAV
        M_X = M_X + A2M * AX
        M_Y = M_Y + A2M * AY

      END DO

      IF ( NP > 0 ) THEN
        ! Add the contributions from the tail, projected on x and y
        U_S(:) = U_S(:) + U_S_tail(:) * CTH
        V_S(:) = V_S(:) + U_S_tail(:) * STH
        M_X = M_X + M_tail * CTH
        M_Y = M_Y + M_tail * STH
      END IF

      ! Let the drift be zero where it may reach below the sea floor
      U_S(numDepths+1:SPND) = 0.
      V_S(numDepths+1:SPND) = 0.

      ! The integral pseudo-momentum per unit surface area, divided by density
      ! [ units m^2/s ]
      M_X = M_X * GRAV
      M_Y = M_Y * GRAV

! Formats
!
  907 FORMAT ('  ',3(A,I6))
  912 FORMAT ('  ',A)

      END SUBROUTINE STOKES_PROFILE


!/ ------------------------------------------------------------------- /
!> @brief PROG_EDGE (A, Pn, m1Bg, uX, uY, JSEA) - prognostic spectral level.
!>      
!>  @details
!>  Calculate the energy spectral level and circular moment from the average
!>  over an interval [IKST,IKT] of bands near NK.
!>  Assume the energy (variance) spectrum F(theta,sig) is symmetric around theta0
!>  @verbatim
!>     F(theta,sig) = H0(th)S(sig)
!>     m1 = \int_-pi/2^pi/2 H0(theta) cos(theta-theta0) d theta
!>     \int_-pi/2^pi/2 H0(theta) d theta = 1
!>  @endverbatim
!>  If Pn == 5, we assumea spectrum of shape
!>  @verbatim
!>     S(sig) = B g^2 * sig^{-5},
!>  @endverbatim
!>  We assume m1 B = constant and return the value ofm1Bg = mean (m1 B g)
!>  If Pn == 4, we assume a spectrum of shape
!>  @verbatim
!>     S(sig) = A g^2 /sig0 * sig^{-4},
!>  @endverbatim
!>  We assume m1 A = constant * sig^{-1} and return m1Bg = mean (m1 A g sig/sig0)
!>  The wave mean direction unit vector (direction of theta0) is returned
!>  as (uX, uY).
!>
!/ ------------------------------------------------------------------- /

      SUBROUTINE PROG_EDGE(A, Pn, m1Bg, uX, uY, JSEA)

!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      IMPLICIT NONE
      REAL, INTENT(IN)        :: A(:,:)   ! Wave action
      REAL, INTENT(IN)        :: Pn ! Negative power of spectral decay (4 or 5)
      REAL, INTENT(OUT)       :: m1Bg, uX, uY
      INTEGER, INTENT(IN)     :: JSEA
      
      REAL, allocatable, save :: DTHsigPn_CgG(:,:)

!     Integral band sums
      REAL                    :: ABX, ABY

      REAL                    :: SUM, SUMX, SUMY
      INTEGER                 :: KSEA, ISEA, ITH, IK, IERR

! On the first call, set the action_to_energy factor / GRAV for calculating
! the spectral level of the tail      
      IF ( .not. allocated(DTHsigPn_CgG) ) THEN
        allocate ( DTHsigPn_CgG(IKT-IKST+1,NSEAL), stat=IERR )

        DO KSEA=1,NSEAL
#ifdef W3_DIST
          ISEA       = IAPROC + (KSEA-1)*NAPROC
#endif
#ifdef W3_SHRD
          ISEA       = KSEA
#endif
          DO IK=IKST,IKT
            DTHsigPn_CgG(IK-IKST+1,KSEA)                             &
                 = DTH * SIG(IK)**(Pn+1) / CG(IK,ISEA) / GRAV
          END DO
        END DO
      END IF

! Integrate over wave energy in bands (as in subroutine CALC_U3STOKES, w3iogomd)
      SUMX=0.
      SUMY=0.
      
      ! Spectral level and circular moment from the average over an
      ! interval [IKST,IKT] of bands near NK
      
      DO IK = IKST,IKT
        ! Wave action in the band (each of the X- and Y-composants)
        ABX    = 0.
        ABY    = 0.
        DO ITH=1, NTH        
          ABX  = ABX + A(ITH,IK)*ECOS(ITH)
          ABY  = ABY + A(ITH,IK)*ESIN(ITH)
        END DO
        ! Spectral level in the band, assuming S(sig) = B g^2 * sig^{-5}
        ! m1Bg = m1s B g = action_to_energy * block_mean( AB ) /GRAV:
        SUMX = SUMX + ABX * DTHsigPn_CgG(IK-IKST+1,JSEA)
        SUMY = SUMY + ABY * DTHsigPn_CgG(IK-IKST+1,JSEA)
      END DO

      SUM = sqrt( SUMX**2 + SUMY**2 )
      ! Spectral level times circular moment
      m1Bg = SUM/(IKT-IKST+1)
      ! Wave mean direction unit vector
      uX = SUMX/SUM
      uY = SUMY/SUM
      
      END SUBROUTINE PROG_EDGE

!/
!/ End of W3STVPMD ----------------------------------------------------- /
!/
      END MODULE W3STVPMD

