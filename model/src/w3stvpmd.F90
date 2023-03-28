!> @file w3stvpmd.F90
!> @brief Vertical profile of the Stokes drift with extended spectral tail.

!> @author Carsten Hansen @date 21-Jan-2021

!> Under development: Remember set TSTOUT = .TRUE. in constants.F90

#include "w3macros.h"
!/ ------------------------------------------------------------------- /
      MODULE W3STVPMD
!/
!/                  +------------------------------------------+
!/                  | Stokes drift  profile with extended      |
!/                  | spectral tail                            |
!/                  |          Carsten Hansen                  |
!/                  | Danish defence Joint GEOMETOC Support    |
!/                  |                            FORTRAN 90    |
!/                  | Last update :             03-Feb-2023    |
!/                  +------------------------------------------+
!/
!/    02-Jan-2023 : Implement option USXT = 'Pf-5': Frequency^{-5} tail with no
!/                  further broadening as function of frequency
!-------------------
!>
!>  @brief Purpose, and output data.
!>
!>  @details
!>  Calculate the vertical profile of the Stokes drift at a few discrete depths,
!>  using an assumed power-law spectral tail (Phillips' frequency^{-5}) for the
!>  diagnostic parts of the spectrum at high frequencies. It assumes that there
!>  is no spectral broadening as function of frequency, i.e. the so-called
!>  'first circular moment' (Ewans, 1997) is assumed constant in the diagnostic
!>  tail. A namelist parameter USXT = 'Pf-5' must be set to add this tail.
!>  
!>  This replaces an earlier version that has a described directional
!>  distribution following the empirical results of Ewans (1997) in conjunction
!>  with an assumed power-law of frequency^{-4}. This could be re-implemented and
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
!>    M_tail = \int_sig_t^sig_c int_-pi^pi kv/sig H(sig,theta) d theta d sig
!>  @endverbatim
!>  Here, sig_t it the tail start frequency, sig_c is a user-specified
!>  high-frequency cut-off, and kv is the spectral wave number vector,
!>  kv = k * (cos(theta), sin(theta)). Assume the tail components are all
!>  deep-water waves. Then k = sig^2/g.
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
!>  magnitude from the spectral tail between the frequencies sig_t and sig_c,
!>  @verbatim
!>    U_S_tail(z=0) = \int_sig_t^sig_c 2 m1s/g sig^3 B g^2 sig^-5 d sig
!>       = 2 m1s B g ( sig(sig_t)^-1 - sig(sig_c)^-1 )
!>        [ units: m s^-1 ]
!>  @endverbatim
!>  The tail contribution to the magnitude of the integral Stokes drift is
!>  @verbatim
!>    M_tail = \int_sig_t^sig_c m1s/g sig B g^2 sig^-5 d sig
!>       = 1/3 m1s B g ( sig(sig_t)^-3 - sig(sig_c)^-3 )
!>        [ units: m s ]
!>  @endverbatim
!>  At depth z>0, the Stokes drift magnitude is an 'incomplete-Gamma' function
!>  defined as the integral
!>  @verbatim
!>    U_S_tail(z) = \int_sig_t^sig_c 2 m1s B g sig^-2 exp(-2 k z) d sig
!>       = \int_k_t^k_c m1s B g^1/2 k^-3/2 exp(-2 k z) d k
!>       = m1s B (2 g z)^1/2 (GammaInc(-1/2, 2 k_t z) - GammaInc(-1/2, 2 k_c z))
!>  @endverbatim
!>  where k_t = sig_t^2/g, k_c = sig_c^2/g .
!>
!>  In the actual numerical approximation we don't need to calculate the
!>  incomplete Gamma function. We can perform the discrete integration as:
!>  @verbatim
!>    U_S_tail(z) = \sum_i=it^ic m1s B g^1/2 k_i^-3/2 exp(-2 k_i z) D k_i
!>  @endverbatim
!>  Now, let the user configure a set of discrete values of dimensionless
!>  depths ZK = z k_s, which is set globally for all sea points at all timesteps,
!>  where k_s is an inverse depth scale. k_s controls the variation over space
!>  and time, and is an output field of the model run. In the present
!>  implementation, k_s is an integral wavenumber scale calculated from the
!>  prognostic spectrum. This means that any possible value of k_s is at
!>  least a half bin increment below k_t, that is k_t/k_s >= XFR.
!>
!>  Let us apply a look-up function NSIZ(IZ, it) where IZ is the index (1:NDP)
!>  of the user-defined discrete values of ZK, and it is the index in (1:NT) of
!>  a discrete set of values of the ratio k_t/k_s, that is where
!>  k_t/k_s = XFR * XKR^(it-1). At the surface (IZ=1) we define
!>  @verbatim
!>    NSIZ(1,it) = int_(k_t/k_s)^infty (k/k_s)^-3/2 d(k/k_s)
!>               = 2 (k_t/k_s)^-1/2
!>               = 2 XFT^-it XFR^1/2, it=1..NT
!>  @endverbatim
!>  If sig_c is finite, then NT should be large enough that NSIZ(IZ, ic) is
!>  defined even for k_s in the lowest prognostic bin.
!>  At finite depths below the surface, we apply the discrete approximation
!>  @verbatim
!>    NSIZ(IZ, it) = sum_i=it^NT (NSIZ(1,it+1)-NSIZ(1,it)) exp(-2ZK(IZ) k_i/k_s)
!>  @endverbatim
!>  where k_i/k_s = XFR * XKR^1/2 * XKR^(i-1), i=it:NT.
!>  If the cut-off frequency sig_c is very high and could be ignored, we have
!>  for the depths z in the set ZK/k_s, that
!>  @verbatim
!>    U_S_tail(z) = m1s B (g/k_s)^1/2 NSIZ(IZ, it) (infinite sig_c).
!>  @endverbatim
!>  If sig_c is finite, then
!>  @verbatim
!>    U_S_tail(k_s z) = m1s B (g/k_s)^1/2 ( NSIZ(IZK, it) - NSIZ(IZ, ic) )
!>  @endverbatim
!>
!>  We set the value of NT as the highest possible value reaching a
!>  user-specified high-frequency cutoff sig_c (sig_c/(2pi) defaults to 10 Hz),
!>  SIG(1)^2*XKR^NT ~ sig_c^2 rad/s, or NT = 2 LOG(sig_c/SIG(1))/LOG(XKR).
!>  In the present implementation we *define* the wavenumber scale k_s as based
!>  on the integral wave period measure T02, k_s = (2 pi/T02)^2/g.
!>  This setting may be changed in the future, thus k_s is provided in the
!>  output as a netCDF variable 'ksc'.
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

      USE W3ODATMD, ONLY: NDST, NDSO, NDSE, IAPROC, NAPROC, NAPOUT
      USE W3DISPMD, ONLY: DSIE, N1MAX, ECG1, EWN1


      PRIVATE

!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/

! DEBUG: File id. for verbose output
      INTEGER            :: NDSV, stvp_verbose = 1
      
! A global limiter to K_S
      REAL               :: K_S_max
      
! Profile array factors
      REAL, allocatable  :: A2S(:), TZ_S(:), ETKZ(:)

! Normalized integral momentum (might be an array if transition freq IKT varies)
      REAL               :: Im

! For calculation of shallow water Stokes drift
      REAL, allocatable  :: fkd(:,:), fkdm(:,:) ! (1:NK,1:NSEAL)
      
! For calculation of exp(2 K(IK) Z) = exp(2 K(IK-1) Z) * exp(2 (K(IK)-K(IK-1)) Z)
! K(IK)-K(IK-1) = WN(IK,JSEA) - WN(IK-1,JSEA)
      REAL, allocatable  :: DWN(:,:)
      
! Drift profile contribution from an extended spectral tail
      REAL, allocatable  :: U_S_tail(:)

! Stokes drift in the diagnostic tail (look-up table for the extended tail).

! IKT:  Prognostic band for tail estimation.
      INTEGER            :: IKT

! sig_c: Frequency at the cut-off frequency USXF
      REAL               :: sig_c
      
! XKR, XFT: Relative increment of wavenumber, frequency bin-to bin for the
!      diagnostic Stokes array (the extended tail)
! NT: Number of wavenumber bins in the full look-up array for the tail profile.
! NP:  Number of bins in the look-up arrays that matches the range 1, IKT
! rr: Number of tail bins per prognostic bin. rr=1 is reasonably precise.
      REAL               :: XKR, XFT
      INTEGER            :: NT
      INTEGER            :: NP
      REAL               :: rr

! Deriving m1Bg from U_S in band IKT by division with the norm:
      REAL               :: NSI0B

! Look-up lists for Stokes spectral tail calculation
      REAL, allocatable  :: NSIZ(:,:) ! (NZ,NT)
      REAL, allocatable  :: SIG_SI(:) ! (1,NP)
      
! Power of the spectral tail decay
      REAL               :: Pn = 5.0

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
      logical             :: OPENED

! Common Stokes drift depth distribution parameters
      REAL                :: DEPTH, XZK, kd, e2kd, shkd4
      REAL                :: kd_t, kd_d, dkd, fkd_t, fkdm_t
      INTEGER             :: IZ
! Upper cutoff frequency for Stokes drift
      REAL                :: sig_c
! Building the normalized stokes tail look-up
      INTEGER             :: ii
      REAL, allocatable   :: DNSI0(:)
      REAL                :: KRt
      
      NDSV = NDST
      ! Highly verbose output only with test output
      IF ( stvp_verbose .gt. 0 ) THEN
        INQUIRE (NDSV,OPENED=OPENED)
        IF ( .NOT. OPENED ) THEN
          NDSV = NDSO
          IF ( .NOT. IAPROC == NAPOUT )  stvp_verbose = 0
        END IF
      END IF

!/ ------------------------------------------------------------- /
! Allocate arrays for depths, prognostic Stokes, and tail Stokes
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

      ! For calculation of exp(2 K(IK) Z) ...
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
      ! The array of dimensionless depths will be constructed as
      !      ZK_S(IZ) = SPDS*(XZK**((IZ-1.)**SPBP) - 1.) for IZ = 1 .. SPND
      ! The exponential base XZK is here defined so that the deepest point is
      ! ZK_S(SPND) = SPDS,

      XZK = 2.**((SPND*1.-1.)**(-SPBP))
      DO IZ = 1,SPND
        ZK_S(IZ) = XZK**((IZ*1.-1.)**SPBP)
      END DO
      ZK_S(:) = SPDS*(ZK_S(:) - 1.)

      ! Calculation for finite water depth of the action-to-stokes factor
      ! A2S(:) in subroutine STOKES_PROFILE. There, the factor is:
      ! A2S[IK,iz] = fkd[IK,JSEA] * ETKZ + fkdm[IK,JSEA] / ETKZ, where
      ! ETKZ = EXP( - 2. * kz[iz] )
      ! The gridded parameters fkd(1:NK,1:NSEAL), fkdm(1:NK,1:NSEAL) are to be
      ! set here with a linear interpolation in a transitional range,
      ! e.g. 2 < kd <= 3, from an exact shallow expression for
      ! kd <= 2 to assumed deep water for kd > 3, where fkd=1, fkdm=0 
      ! Note: At the surface where kz[1] = 0 we have ETKZ = 1, and 
      !       A2Sf[IK,1] = 0.5 * COSH(2.*kd) / SINH(kd)**2
      !                  = fkd[IK,JSEA] + fkdm[IK,JSEA]
      !
      ! Start fill the arrays as if all points and bins are deep water waves,
      fkd = 1.
      fkdm = 0.
      
      ! Apply a linear interpolation in the transitional range kd_t < kd <= kd_d
      ! of shallow water waves to deep water. Our choice is:
      kd_t = 2. 
      kd_d = 3.
      dkd = kd_d-kd_t

      ! At the lower wavenumber kd_t of this range, fkd=fkd_t and fkdm=fkdm_t are
      ! calculated as:
      e2kd = exp (2.*kd_t)
      shkd4 = 4. * (sinh( kd_t ))**2
      fkd_t = e2kd/shkd4
      fkdm_t = 1./ (e2kd * shkd4)

      DO JSEA = 1,NSEAL
#ifdef W3_DIST
        ISEA       = IAPROC + (JSEA-1)*NAPROC
#endif
#ifdef W3_SHRD
        ISEA       = JSEA
#endif

        DEPTH  = MAX ( DMIN, DW(ISEA) )

        ! The parameters fkd, fkdm are:
        DO IK=1,NK
          kd = WN(IK, ISEA) * DEPTH

          ! kd>kd_d: pure deep water (then keep the filled fkd=1 and fkdm=0)
          IF (kd .gt. kd_d) exit

          ! kd_t < kd <= kd_d: Linear interpolation of transitional to deep water
          IF (kd .gt. kd_t) THEN
            ! Re-use the factors fkd_t, fkdm_t calculated for kd = kd_t
            fkd(IK,JSEA) = ( fkd_t * (kd_d - kd) + (kd - kd_t) ) / dkd
            fkdm(IK,JSEA) = ( fkdm_t * (kd_d - kd) ) / dkd
            CYCLE
          END IF

          ! kd<kd_t: shallow water
          e2kd = exp (2.*kd)
          shkd4 = 4. * (sinh( kd ))**2
          fkd(IK,JSEA) = e2kd / shkd4
          fkdm(IK,JSEA) = 1./ (e2kd * shkd4)
          
        END DO ! IK=1,NK

        ! For calculation of exp(-2 K Z)
        ! DWN(IK) = K(IK) - K(IK-1)
        DWN(1,JSEA) = WN(1,ISEA)
        DWN(2:NK,JSEA) = WN(2:NK,ISEA) - WN(1:NK-1,ISEA)

      END DO ! JSEA = 1,NSEAL
      
      ! Apply a global limit to K_S
      K_S_max = SIG(NK - 2)**2 / GRAV
      ! In the fittting of an extended tail at index IKT it will be required
      ! that K_S_max <= SIG(IKT)**2 / GRAV. Ensure also that this is so
      ! K_S_max = min( K_S_max, SIG(IKT)**2 / GRAV )


!/ ------------------------------------------------------------- /
! Initializations for calculations related to the diagnostic tail range

      ! IKT: The upper edge of the prognostic range.
      ! NP: Number of bins of the lookup array for the spectral tail extension
      !     that may match at the end of the prognostic range
      ! For a truncated spectrum with no tail extension:
      IKT = NK
      NP = 0 ! NP=0 indicates no addition of a spectral tail

      ! User-specified upper cutoff frequency for Stokes drift. At most 10 Hz
      sig_c = TPI * MIN(USXF, 10.)
      
      ! If the user-specified cutoff for Stokes drift is within SIG(1..NK),
      if ( sig_c < SIG(NK) ) then
        do IKT=NK-1,1,-1
          if ( sig_c > SIG(IKT) ) exit
        end do
        ! End the prognostic summation at IKT with no addition of a spectral tail
        K_S_max = min( K_S_max, SIG(IKT)**2 / GRAV )
        return
      end if      
      
! Fitting the tail at IKT=NK-1
      IF ( USXT == 'Pf-5' ) THEN
         ! Phillips saturation
         Pn = 5
         IKT = NK - 1
      ELSE IF ( USXT == 'none' ) THEN
         RETURN ! No tail extension. The value NP = 0 will indicate this.
      ELSE
        WRITE (NDSE, 916), 'STVP_INIT: Unknown spectral tail type: ', USXT
        CALL EXTCDE ( 73 )
      END IF

! Let XKR be the relative wavenumber increment of the lookup array for the
! extended tail part of the Stokes profile. In the present formulation, the
! choice of XKR must be an integer root r: XKR = (XFR**2)**(1/r). We choose
! r = 1.     
      rr=1
      XKR = (XFR**2)**(1./rr)
! We will later make use of the frequency increment factor for the look-up tail
      XFT = XFR**(1./rr)
      
! We should have in mind that when the depth is limited, the elements of the array SIG_SI
! are NOT (precisely) equal to the true angular frequencies SIG(1..NK)
      
! The length of the look-up table NT must be so that
! sig_c/SIG(1) <= XKR**(NT/2) (at least).
      NT = ceiling( 2. * LOG( sig_c/SIG(1) ) / LOG(XKR) ) + 1

      allocate( NSIZ(SPND, NT), stat=IERR )
      
! We let an array WN_SI(1:NP) be the possible start wavenumber value of the
! discrete look-up array for the extended tail. NP may be an uneven integer
! multiplum r of the prognostic bin arrays, and we require that
! WN_SI(1:NP:r) = WN(1:IKT) .
! To simplify the notations we will define an array (denoted 'frequencies')
! SIG_SI(1:NP:r) = sqrt(GRAV*WN_SI(1:NP)).

! We will further make the approximation, that 
! The highest possible look-up tail frequency start value (SIG_SI(NP)) is at
! the prognostic bin IKT, and we will further make the approximation, that
! the full tail, including SIG_SI(NP), have a dispersion relation as for deep
! water waves. This means that we may apply SIG_SI(NP) = SIG(IKT), and
! NP = (IKT-1)*rr + 1
      NP = IKT + (1 - rr)

      allocate( SIG_SI(NP), stat=IERR )

! The look-up array of possible tail frequency start values SIG_SI(1:NP) start
! at SIG(1)
      SIG_SI(1) = SIG(1)
      do ii=2,NP
         SIG_SI(ii) = SIG_SI(ii-1) * XFT
      end do
      
! The integrated pseudo-momentum of the diagnostic spectral tail is
! M_tail = m1Bg * Im. The tail is added from the upper edge of bin IKT.
      Im = 1./(3. * (SIG(IKT)*(SQRT(XFR)))**3)

! Subtract the range beyond the cut-off frequency
      Im = Im - 1./(3. * sig_c**3)
      
! Now we construct of NSIZ(IZ,it): A 'normalized, look-up' vertical Stokes
! drift profile for spectral bins of the diagnostic tail.
!
! Let the scale parameter K_S be defined in terms of one of the usual
! integral period measures. As input to the subroutine
! STOKES_PROFILE we choose the 2'nd moment upcrossing period T02,
! so that the phase speed term becomes (g/K_S)^1/2 = g T02 / (2 pi).

! The tail Stokes drift profile is estimated from a look-up table NSIZ(IZ, it),
!   U_S_tail(Z) = m1B (g/K_S)^1/2 NSIZ(IZ, it)
! where
!   NSIZ(IZ, it)= \sum_ii=it^NT DNSI0(ii) exp (-2k/K_S ZK_S)
!   DNSI0(ii) = (k/K_S)^-3/2 D k/K_S
!
! The look-up is constructed for k(ii) being elements in an array of discrete 
! bins. This means that we seek to discretize the ratio k/K_S into an array:
!   k(ii)/K_S ~= Kr(ii) = ( XKR^1, XKR^2, ..., XKR^ii, ... ), ii=1..NT
! This construction is to be matched at the highest prognostic band
! IKT, and applied to the tail bin at k_m = k_t/XFR = WN(IKT)*XKR**((rr-1)/2).
! We will determine the index it so that k(it) equals k_m * XKR.
! Let K_Sd be the nearest discrete value below or equal to K_S. Kr is taken at the
! index it where k_m * XKR = K_Sd * Kr(it). This means that the first element it=1 is
! Kr(1) = XKR**((rr+1)/2) = XFT**(rr+1)
! For example, if rr = 1 and there were energy at only one prognostic bin
! (i.e. at k_m only), then K_S = K_Sd = k_m precisely, and we would apply the
! look-up array at it = 1.
      
! The look-up table value at the surface, Z=0, NSIZ(1,it), is
! NSIZ(1,it) = int_(k_t/K_Sd)^inf (k/K_S)^-3/2 d(k/K_S)
!            = 2 (k_t/K_Sd)^-1/2
!            = 2 (k(it)/XFT/K_Sd)^-1/2
! Thus the first element needed is NSIZ(1,1) = 2 (Kr(1)/XFT)^-1/2
!    = 2 ( XFT**(rr+1)/XFT )^-1/2 = 2 XFT**(-(rr-1)/2)
      NSIZ(1,1) = 2.
      IF ( rr > 1 ) NSIZ(1,1) = 2. * XFT**(-(rr-1)/2)
      DO ii = 2,NT
          NSIZ(1,ii) = NSIZ(1,ii-1) / XFT
      END DO
      
! Each bin interval contributes to NSIZ at Z=0 as
! DNSI0(it) = int_(k_i/K_S/XFT)^(k_i/K_S*XFT) (k/K_S)^-3/2 d(k/K_S)
!           = 2 ( (k_i/K_S/XFT)^-1/2 - (k_i/K_S/XFT)^-1/2 / XFT )
!           = NSIZ(1,it) * (1-1/XFT), it=1..NT
      
      allocate( DNSI0(NT), stat=IERR )
      
      DNSI0(:) = NSIZ(1,:)* ( 1. - 1./XFT )

      ! NSIZ is a cumulated sum downwards from the high end of the spectral tail
      ! Any value beyond sig_c will be subtracted, thus we can safely ignore
      ! the rest of the tail. For the surface, we keep the values NSIZ(1,:)
      KRt = -2. * XKR**NT ! starting at ii = NT
      NSIZ(2:,NT) = 0.
      DO ii=NT-1,1,-1
         KRt = KRt / XKR ! KRt = -2 * XKR**ii
         NSIZ(2:,ii) = NSIZ(2:,ii+1) + DNSI0(ii) * exp( KRt * ZK_S(2:) )
      END DO      

! DNSI0 is not needed further
      deallocate ( DNSI0 )      

! The tail will be matched to the contribution to the Stokes drift, U_SB_t,
! from the highest prognostic bin with wavenumber WN(IKT).
! At the water surface (Z=0) the scale K_S is ambigious and may be chosen as
! K_S=WN(IKT). This prognostic bin is covered by a range just below NSIZ(1, 1),
!
! U_SB_t = m1 B g (g WN(IKT))^-1/2 ( NSIZ(1, 1)*XFR - NSIZ(1, 1) )
! We will keep the last factor for repeated use at all water points:
       NSI0B =  NSIZ(1, 1)*(XFR-1)
       
! During the model run, at all output time steps, the Stokes drift U_SB_t in
! the prognostic band IKT will be calculated at all water points. The Stokes
! spectral level is then derived as (see SUBROUTINE STOKES_PROFILE)
! m1Bg =  U_SB_t * (g WN(IKT))^1/2 / NSI0B
       
! Application procedure in the subroutine STOKES_PROFILE():
! 1) Determine the level m1 B g from the band IKT.
! 2) For a choice of SIG_S = TPI/T02 based on the integral period T02,
!    determine the nearest integer it so that sqrt(g k_m)/SIG_S~=XFT**(it-1),
!    and the tail cut-off index ic so that sig_c/SIG_S~=XFT**(ic-1)
! 3) The diagnostic tail contribution to the drift profile is, in case
!    SIG_S is one of the binned values, 
!        U_S_tail(:) = m1 B g / SIG_S * ( NSIZ(:,it) - NSIZ(:,ic) )
!    For a value of SIG_S between two bins, apply a linear interpolation
!    between the values of SIG_S for which step (2) would yield the neighbour
!    indices it-1, it.
! 4) The total Stokes drift is
!        M_tail = m1Bg * Im

  904 FORMAT ('  ', A,I5,A)
  916 FORMAT ('  ', A, A)

      END SUBROUTINE STVP_INIT


      SUBROUTINE CALC_STVP (A)
!/ ------------------------------------------------------------------- /
!> @brief CALC_STVP (A) - Calculate the profiles looping over all sea points
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
!> @brief STOKES_PROFILE (A, JSEA, PERIOD) - Calculate the profile and pseudo-momentum
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
!>  measure (PERIOD equal to i.e. T02) for the time step, and the actual depths
!>  Z(1:SPND) can be determined from the global dimensionless dimension
!>  ZK_S(1:SPND) as
!> @verbatim
!>         Z(1:SPND) =  1/K_S * ZK_S(1:SPND)
!> @endverbatim 
!>  The Stokes profile is a sum of a prognostic part and a tail part
!>  U_S_tail(1:SPND)
!/ ------------------------------------------------------------------- /
      SUBROUTINE STOKES_PROFILE(A, JSEA, PERIOD)
!/
      USE W3GDATMD, ONLY: MAPSF
      USE W3SERVMD, ONLY : EXTCDE
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/

      REAL, INTENT(IN)            :: A(NTH,NK)
      INTEGER, INTENT(IN)         :: JSEA
      REAL, INTENT(IN)            :: PERIOD
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      REAL                  :: SIG_S, DEPTH, KDPT
      INTEGER               :: numDepths, IZ

! Factors for the prognostic calculation
      REAL                  :: A2M, A2S0, AX, AY, EDWTZ, DDWTZ
      INTEGER               :: ISEA, ITH, IK, NIK

! Parameters for 'inlined version of WAVNU1'
      REAL                  :: SIX, R1
      INTEGER               :: I1

! Tail spectral level from prognostic band at IKT
      REAL                  :: m1Bg, U_SBX, U_SBY, U_SB
! Look-up tail index
      INTEGER               :: it ! Tail look-up index
      INTEGER               :: ip ! = NP - it
! Look-up tail and interpolation parameters
      REAL                  :: SIG_DS, wgt
      
! Tail integral pseudo-momentum and tail direction parameters
      REAL                  :: M_tail, CTH, STH
      
! Numerical consistency check output
      LOGICAL               :: INDTEST = .FALSE.
      CHARACTER(LEN=30)     :: rowfmt
      INTEGER               :: INDTESTX = 0, INDTESTY = 0
      ! To generate output to NDSE for verification of the central parameters
      ! for numerical consistency (correct indexing etc.), choose INDTESTX,
      ! INDTESTY as sea point values /=0. This will stop the program at the
      ! first instance.
      ! Set e.g. INTEGER               :: INDTESTX = 114, INDTESTY = 125

#ifdef W3_DIST
      ISEA       = IAPROC + (JSEA-1)*NAPROC
#endif
#ifdef W3_SHRD
      ISEA       = JSEA
#endif

! Numerical consistency check output
      if ( MAPSF(ISEA,1) == INDTESTX .and. MAPSF(ISEA,2) == INDTESTY ) &
        INDTEST = .TRUE.
      
! SIG_S is an integral wave frequency chosen as SIG_S == 2pi/PERIOD
      SIG_S = TPI/PERIOD

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

      ! If a tail is appended, a frequency scale SIG_DS for the tail look-up is
      IF ( NP > 0 ) THEN
        IF (I1 .LE. N1MAX) THEN
          SIG_DS = SQRT(GRAV * K_S)
        ELSE
          SIG_DS = SIG_S
        END IF
      END IF
      
! Apply the global limiter to K_S
      K_S = MIN(K_S_max,K_S)

      ! numDepths: Number of depths of the discrete profile above the sea floor

      ! SPND is the configured, maximum number of profile depths
      ! The dimensionless water depth is
      KDPT = K_S * DEPTH

      DO IZ = SPND,1,-1
        IF ( ZK_S(IZ) > KDPT ) CYCLE
        numDepths = IZ
        EXIT
      END DO

      ! Prognostic spectrum range IK=1, IKT

      IF ( stvp_verbose .gt. 1 ) THEN
        WRITE (NDSV, 912), '    .. for the prognostic part'
      END IF

      ! We represent the Stokes profile exponential shape exp(-2 K Z) by an
      ! array ETKZ(Z) which depends on the prognostic wavenumber bin IK. An
      ! initial value of K may be K=K0=0. For each of the prognostic wavenumbers
      ! K in WN(IK=1..NK), the exponential shape over Z(:) is calculated as 
      ! ETKZ(:)[IK=1] = exp( -2 WN(1) Z(:) )
      !               = exp( -2 K0 Z(:)) exp( -2 (WN(1)-K0) Z(:) ),
      ! ETKZ(:)[IK=2] = exp( -2 WN(2) Z(:) )
      !               = exp( -2 WN(1) Z(:) ) exp( -2 (WN(2)-WN(1)) Z(:) )
      !               = ETKZ(:)[IK=1] exp( -2 (WN(2)-WN(1)) Z(:) ),
      ! ETKZ(:)[IK=3] = exp( -2 WN(3) Z(:) )
      !               = ETKZ(:)[IK=2] exp( -2 (WN(3)-WN(2)) Z(:) ),
      ! etc.
      ! The increments DWN(IK,JSEA) = WN(IK+1,ISEA) - WN(IK,ISEA) are
      ! determined in STVP_INIT() for all local seapoints JSEA
      !
      ! The profile for the initial K0=0 is exp( -2 K0 Z(:) ) = exp( 0 )
      ETKZ(:) = 1

      ! Note that at the surface (Z=0) we have for all IK that ETKZ(1) == 1
      
      ! Initialize the summation over the prognostic bins
      U_S = 0.
      V_S = 0.
      M_X = 0.
      M_Y = 0.

      ! An array holding -2 Z(:) is
      TZ_S(:) = -2.0/K_S * ZK_S(:)

      ! Numerical consistency check output
      if ( INDTEST ) then
        write(NDSE,912), 'NT, SPND, IKT ='
        write(NDSE, FMT='(3(1X,I2))'), NT, SPND, IKT
        write(NDSE,912), 'IK, log10(D M_X/Y), log10(D U/V_S(IK,:)) ='
        WRITE(rowfmt,'(A,I4,A)') '(1X,I2,1X,F6.3,',size(U_S),'(1X,F6.3))'
      end if
      
      DO IK=1, IKT

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

        ETKZ(:) = ETKZ(:) * EXP( DWN(IK,JSEA) * TZ_S(:) )

        ! Omnidirectional action spectrum at frequency bin IK projected on x and y
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
          ! STVP_INIT() with a linear interpolation in a transitional
          ! range, e.g. 2 < kd <= 3, from an exact shallow expression for
          !  kd <= 2 to assumed deep water for kd > 3, where fkd=1, fkdm=0 
          A2S(:) = A2S0 * ( fkd(IK,JSEA) * ETKZ(:) + fkdm(IK,JSEA) / ETKZ(:) )
        END IF
        ! Stokes drift projected on grid-Eastern and grid-Northern components
        U_S(:) = U_S(:) + AX * A2S(:)
        V_S(:) = V_S(:) + AY * A2S(:)

        ! Integral pseudo-momentum per unit surface area, divided by rho_w*GRAV
        M_X = M_X + A2M * AX
        M_Y = M_Y + A2M * AY
        
        ! Numerical consistency check output
        if ( INDTEST ) then
           write(NDSE, FMT=rowfmt ) IK, log10(max(A2M * AX,1.e-09)), &
                ( log10(max(AX * A2S(IZ),1.e-09)), IZ=1,size(U_S) )
           write(NDSE, FMT=rowfmt ) IK, log10(max(A2M * AY,1.e-09)), &
                ( log10(max(AY * A2S(IZ),1.e-09)), IZ=1,size(U_S) )
        end if
      END DO

      IF ( NP == 0 ) THEN

        IF ( stvp_verbose .gt. 0  .and. JSEA .eq. 1 ) &
          WRITE (NDSV, 912), 'Truncate at NK-1. No diagnostic spectral extension'

      ELSE IF ( AX .gt. 0. .or. AY .gt. 0. ) THEN
        ! (If the wave action is zero at the highest prognostic frequency band,
        ! then there is a zero tail contribution)
         
        IF ( stvp_verbose .gt. 1  .and. JSEA .eq. 1 ) &
             WRITE (NDSV, 912), '    .. for the diagnostic part'

        ! We append the tail at frequencies above the prognostic bin IKT.   
        ! Here we estimate the product m1Bg of the 1st circular moment m1, and
        ! the spectral level Bg, and also estimate the band mean *vector*
        ! orientation CTH, STH:(cos(theta0),sin(theta0)) m1 B g = (CTH,STH)*m1Bg
        U_SBX = AX * A2S0 ! The contribution to U_S, VS from band IKT
        U_SBY = AY * A2S0 !
        U_SB = SQRT(U_SBX**2 + U_SBY**2)        
        
        ! At the highest prognostic frequency band, the Stokes spectral level
        ! is derived as explained in SUBROUTINE STVP_INIT:
        if (INT(I1 * SIG(IKT) / SIG_S) .LE. N1MAX) then
          ! For shallow depth
          m1Bg = U_SB * SQRT(GRAV * WN(IKT,ISEA)) / NSI0B
        else
          ! Deep water
          m1Bg = U_SB * SIG(IKT) / NSI0B
        end if
        
        ! Band mean direction unit vector
        CTH = U_SBX/U_SB
        STH = U_SBY/U_SB

        ! The integrated pseudo-momentum of the diagnostic tail has the magnitude
        M_tail = m1Bg * Im

        ! The diagnostic Stokes drift profile, at dimensionless depths Z K_S,
        ! has the magnitude U_S_tail(:).

        ! We introduced the approximation that above the prognostic bin IKT the
        ! dispersion relation is for deep water waves. We defined a frequency
        ! scale SIG_DS for the tail look-up, SIG_DS = SQRT(GRAV * K_S).
        ! 
        ! With this approximation the tail contribution to the surface Stokes
        ! drift is overestimated by a slight value.
        !
        ! For limited depths we have K_S = SIG_S**2 / GRAV / tanh(K_S D).
        ! Therefore SIG_DS is (slightly) larger than SIG_S, and thus always
        ! SIG_DS > SIG_S > SIG(1) as required for the look-up. For the look-up
        ! we also require that SIG(IKT)/SIG_DS >= 1. This is why we have set
        ! a limiter K_S_max = SIG(IKT)^2/GRAV above.
        
        ! Interpolate between the two discrete values of SIG_SI nearest
        ! SIG_DS where the tail matches most precisely.
        
        do ip=2,NP
          if (SIG_SI(ip) > SIG_DS) exit
        end do
        ip = ip - 1

        ! Perform linear interpolation of the tail Stokes drift between the
        ! two neighbour points ip and ip+1.

        ! The linear weights are
        wgt = ( SIG_SI(ip+1) - SIG_DS )/( SIG_SI(ip+1) - SIG_SI(ip) )

        ! At the neighbour points the look-up values of U_S_tail are
        !   m1Bg * NSIZ(:,it+1) / SIG_SI(ip)   ! if SIG_DS = SIG_SI(ip)
        ! and
        !   m1Bg * NSIZ(:,it) / SIG_SI(ip+1)   ! if SIG_DS = SIG_SI(ip+1)
        !   = m1Bg * NSIZ(:,it) / SIG_SI(ip) / XFT
        ! where
        it = NP - ip

        ! We further subtract the high-frequency cut-off beyond sig_c to
        ! a reasonable accuracy, and get
        !
        U_S_tail(:) = m1Bg * (   NSIZ(:,it+1)  * wgt         &
                               + NSIZ(:,it) / XFT * (1.-wgt) &
                               - NSIZ(:,NT-NP+it+1)          &
                             ) / SIG_SI(ip)

        ! Examples:
        !       if SIG_DS = SIG_SI(ip), then apply NSIZ(:,it+1), wgt = 1
        !       if SIG_DS = SIG_SI(ip = NP-1), then apply NSIZ(:, it = 1).
        !       if SIG_DS = SIG_SI(ip = 1), then apply NSIZ(:, it+1 = NP)
        !
        ! Note, at Z = 0 the result is rather simple, because
        ! NSIZ(1,ii) = 2. / XFT**(ii-1) and SIG_SI(ip) = SIG_SI(NP) / XFT**it.
        ! ( Remember that in the deep water limit, SIG_SI(NP) = SIG(IKT) )
        ! Thus, ignoring first the cut-off beyond sig_c, we have
        ! U_S_tail(1) = m1Bg *( NSIZ(1,it+1) )/ SIG_SI(ip) (minus cut-off)
        !             = 2 * m1Bg / XFT**it / ( SIG_SI(NP) / XFT**it ) (minus cut-off)
        !             = 2 * m1Bg / SIG_SI(NP) (minus cut-off)
        !
        ! and the cut-off beyond sig_c to be subtracted is
        ! cut-off = m1Bg * NSIZ(1,NT-NP+it+1))/ SIG_SI(ip)
        !         = 2 * m1Bg / XFT**(NT-NP) / SIG_SI(NP)
        !         

        ! Numerical consistency check output
        if ( INDTEST ) then
          write(NDSE,912), 'SIG_DS, m1Bg, SIG_SI(NP), it ='
          write(NDSE, FMT='(3(1X,F6.4),1X,I2)'), SIG_DS, m1Bg, SIG_SI(NP), it
          WRITE(rowfmt,'(A,I4,A)') '(1X,I2,',size(U_S),'(1X,F6.3))'
          write(NDSE,912), 'IK, log10(NSIZ(:,IK)) ='
          
          do it = 1, NT
            write(NDSE, FMT=rowfmt ) it, &
                 ( log10(max(NSIZ(IZ,it),1.e-06)), IZ=1,size(U_S) )
          end do
          WRITE(rowfmt,'(A,I4,A)') '(',size(U_S),'(1X,F6.3))'
          write(NDSE,912), 'log10(U_S_prog) ='
          write(NDSE,FMT=rowfmt) ( log10(max(U_S(IZ),1.e-09)), IZ=1,size(U_S) )
          write(NDSE,912), 'log10(V_S_prog) ='
          write(NDSE,FMT=rowfmt) ( log10(max(V_S(IZ),1.e-09)), IZ=1,size(U_S) )
          WRITE(rowfmt,'(A,I4,A)') '(',size(U_S_tail),'(1X,F6.3))'
          write(NDSE,912), 'log10(U_S_tail) ='
          write(NDSE,FMT=rowfmt) &
                ( log10(max(U_S_tail(IZ)* CTH,1.e-09)), IZ=1,size(U_S_tail) )
          write(NDSE,912), 'log10(V_S_tail) ='
          write(NDSE,FMT=rowfmt) &
                ( log10(max(U_S_tail(IZ)* STH,1.e-09)), IZ=1,size(U_S_tail) )
          write(NDSE,FMT='(A,F9.6)'), 'M_X = ', M_X
          write(NDSE,FMT='(A,F9.6)'),  'cos M_tail = ', M_tail * CTH
          write(NDSE,FMT='(A,F9.6)'), 'M_Y = ', M_Y
          write(NDSE,FMT='(A,F9.6)'),  'sin M_tail = ', M_tail * STH

          call extcde(988,MSG="Stop after point output for graphics", &
               FILE="w3stvpmd.F90")
        end if

        ! Add the contributions from the tail, projected on x and y
        U_S(:) = U_S(:) + U_S_tail(:) * CTH
        V_S(:) = V_S(:) + U_S_tail(:) * STH
        M_X = M_X + M_tail * CTH
        M_Y = M_Y + M_tail * STH
        
      END IF ! ( NP == 0 ) .. ELSE

      ! Let the drift be zero where it may reach below the sea floor
      U_S(numDepths+1:SPND) = 0.
      V_S(numDepths+1:SPND) = 0.

      ! The integral pseudo-momentum per unit surface area, divided by density
      ! [ units m^2/s ]
      M_X = M_X * GRAV
      M_Y = M_Y * GRAV      
      
! Formats
!
  912 FORMAT ('  ',A)

      END SUBROUTINE STOKES_PROFILE

!/
!/ End of W3STVPMD ----------------------------------------------------- /
!/
      END MODULE W3STVPMD

