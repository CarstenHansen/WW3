!> @file w3ounfmetamd.F90
!> @brief Extended spectral tail for Stokes drift profile.
!> @author Carsten Hansen @date 21-Jan-2021
#include "w3macros.h"
!/ ------------------------------------------------------------------- /
      MODULE W3STVPMD
!/
!/                  +------------------------------------------+
!/                  | Extended spectral tail for Stokes drift  |
!/                  |          Carsten Hansen                  |
!/                  | Danish defence Joint GEOMETOC Support    |
!/                  |                            FORTRAN 90    |
!/                  | Last update :             21-Jan-2021    |
!/                  +------------------------------------------+
!/
!/
!/    02-May-2011 : Origination. Development program w3stokes.ftn converted
!/                  from tested Python code
!/    27-Jun-2011 : First module version w3stvpmd.ftn. Output tested with
!/                  a stand-alone program w3stokes.ftn. Implementation plan.
!/    18-Aug-2011 : Calculate fit to J. Mattsson parametric profile (ver.0.12)
!/    10-Oct-2011 : Interface subroutine w3dstk() to ww3_shel (version 0.22)
!/    21-Oct-2011 : Detailed description of implementation into WW3 v.3.14
!/                  First compilation of ww3_shel in parallel mode
!/                  Test-scripts run_ww3test.sh,
!/                     $PERL5LIB/bin/ww3job_preproc (version 0.32)
!/    08-Nov-2011 : Full functioning module. Un-physical spatial discontinuities
!/                  identified: May need to match the Stokes profiles to
!/                  integral wave pseudo-momentum/depth-integrated Stokes drift.
!/    23-Nov-2011 : Ajust the deep profile to match the depth-integrated
!/                  Stokes drift. Function AJUST_FIT() (version 0.33).
!/    30-Nov-2011 : Re-formulated the look-up tables for the diagnostic tail.
!/                  Simplified tail extension (version 0.40)
!/    07-Dec-2011 : Refined ajustment of Mattsson profiles to integral
!/                  pseudo-momentum (version 0.41).
!/    16-Dec-2011 : Calculate integrated pseudo-momemtum directly (ver. 0.44).
!/                  For the prognostic spectrum, the momentum density is
!/                    m(K) = rho G K A(K), where K is the wave number vector
!/                  and A(K) is the action spectrum.
!/                  For the diagnostic tail, we apply a look-up array Im(n),
!/                    Im(n) = sum_fhn^fhcut m1t(fh) fh^-2 log(DFH)
!/                  where
!/                    fhn = ft/fP = DFH**(n-1), n = 1, ..., Np
!/                  m1t(fh) is the 1'st circular moment in Ewans' formulae,
!/                  and fh = f/fP.
!/                  The integrated pseudo-momentum (per m^2) from ft to fcut is
!/                    M(ft) = rho G alpha U sigP^-2 Im(n), where
!/                  n = log(ft/fP)/log(DFH) + 1,
!/                  sigP = 2 pi fP is the peak frequency in Ewans' formulae,
!/                  G is gravitational acceleration, alpha is the spectral
!/                  constant of Donelan et al., 1985.
!/    22-Dec-2017 : Adoption to version WW3 v5.16:
!/                  + Call CALC_STVP() from W3OUTG(), outside loop over JSEA
!/    09-Mar-2020 : Adoption to version WW3 v7.??
!/                  + Mattsson fitting is separated to module 'W3XSMFMD'
!/                  + Namelist group &STVP
!/                  + from W3OUTG pass the array A of spectra to CALC_STVP(A)
!/                  + Collect technical doc. in comments in diag_DE_init()
!/                  + For the extended tail to be representative, the local wind
!/                    sea peak must be within the prognostic range, fp < IKST
!/    21-Jan-2021 : Match the exact expressions in C. Hansen, 2020 (DRAFT):
!/                  "Mean Stokes drift profile from a wave model with a
!/                  parametric extension of the spectral tail".
!-------------------
!>
!>  @brief Purpose, and output data.
!>
!>  @details
!>  Calculate the vertical profile of the Stokes drift, using a
!>  directional distribution (Ewans, 1997) for the diagnostic parts of the
!>  spectrum at high frequencies in conjunction with an assumed power-law
!>  spectral tail.
!>
!>  All calculations are performed in parallel and the Stokes drift is calculated
!>  at NDEPTH discrete depths given in an array Z_S(1:NDEPTH).
!>  The sum of prognostic and diagonostic (tail) contributions are stored in
!>  arrays U_S(1:NDEPTH), V_S(1:NDEPTH).
!>
!>  @brief Usage.
!>
!>  @details
!>  In ww3_grid.inp/ww3_grid.nml, you may modify five namelist variables:
!>  @verbatim
!>  &STVP NDP = 11, DSC = 2.0, BP = 2.0 /
!>  &XSTP USXT = 'DoEw' , USXF = 2.0 /
!>  @endverbatim
!>  - NDP: Number of depths (->NDEPTH)
!>  - DSC: Depth scale specifying the largest depth Z(NDP)=DSC/(2pi/Tz)^2*GRAV
!>  - BP: Power of profile depth progression        
!>  - USXT: Tail type (C*4); 'DoEw': Donelan-Ewans,
!>  -                        'none': Prognostic spectrum frequencies, only
!>  - USXF: Tail truncation frequency(low-pass). Default 10 Hz (~infinity)
!>
!>  @brief Estimation of Stokes drift of the diagnostic tail.
!>
!>  @details
!>  The Stokes drift is determined by the wave spectrum H(sig,theta) over
!>  (intrinsic) frequencies sig (= 2 pi f in the absense of currents) and
!>  directions theta. For wind waves over deep water, and following
!>  Ewans (1997), the contribution in a frequency bin dsig to the
!>  directionally integrated Stokes drift is 
!>  @verbatim
!>    d Stokes(sig) = 2 g sig^3 m1(sig) S(sig) d sig
!>  where the directionally integrated spectrum is
!>    S(sig) = int_-pi^pi H(sig,theta) d theta
!>  @endverbatim
!>
!>  Following Ewans (1997) m1(sig) is denoted the 'first circular moment' of
!>  the spectrum H(sig,theta),
!>  @verbatim
!>    m1(sig) = int_-pi^pi cos(theta) H(sig,theta) d theta / S(sig)
!>  @endverbatim
!>
!>  Based on observational data (the 'Maui' experiment) Ewans (1997) suggests
!>  a mathematical model to fit the data, which results in the
!>  following expression for the first circular moment:
!>  @verbatim
!>    m1E(sig) = cos(theta_E(f/fp)) exp( -(Sigma(f/fp)**2)/2 ),
!>  @endverbatim
!>  where f/fp is the ratio of local spectral frequency f=sig/2pi over the
!>  peak frequency fp. Ewans provides expressions for theta_E, the reflection
!>  from the mean direction of two 'binodal' directional peaks, and Sigma,
!>  the standard deviation of Gaussian-shaped distributions around each of
!>  the two peaks.
!>
!>  In the present usage ('Donelan-Ewans': subroutine diag_DE_init(NP)),
!>  Ewans' m1 is approximated by a function,
!>  @verbatim
!>    m1t(f/fp) = m1E(f/fp=1) * exp( a sinh(phih * (f/fp - 1)) )
!>  @endverbatim
!>  where a = 1.25 and phih = -0.187 is a negative factor representing the
!>  increasing directional spread with increasing f/fp.
!>
!>  A formula for the spectral tail is assumed in the form of Donelan et. al
!>  (1985),
!>  @verbatim
!>    S(sig) = alpha U10 g sig^-4,
!>  @endverbatim
!>  where the factor alpha may vary with the integral wave properties, among
!>  which the dependance on wave age (U10/cP)*cos(theta_mean-theta_U) has been
!>  investigated. In the present usage, alpha is estimated based on a few of
!>  the highest frequency bins of the prognostic spectrum.
!>  We note that the resulting estimates of alpha is sensitive to the choice of
!>  formulation of wind input. It was first tested (and tuned) with the
!>  Tolman-Chalikov input terms formulation ( WW3 switch '/ST2'), and is
!>  also consistent with ST4-Romero or ST6, which are supposed with certain
!>  parameter combinations to have an sig^-4 spectral tail up to f/fp ~ 3-4.
!>
!>  With the spectrum of Donelan et. al and Ewans' formula for m1(sig), the
!>  integrated Stokes drift int_sigt^inf d Stokes(sig) will *not* converge,
!>  because the binodal direction theta_E(f/fp) converges to a value slightly
!>  larger that pi/2 as f/fp -> infinity. However, it is realistic that
!>  at values of f/fp in the range 3 to 6 there exist a further transition to
!>  a spectrum of satutated form, going asymptotically as sig^-5. Such
!>  spectral tail works as a low-pass filter, and the exact shape of m1 beyond
!>  frequencies of order f/fp ~ 6 does only influence the integrated Stokes
!>  drift by a small amount. The applied approximation m1t(sig) to Ewans'
!>  formula for m1(sig) provides a low-pass filter at an equivalent frequency
!>  scale. Thus it is suggested that there will be little effect of applying
!>  an explicit formula of an sig^-5 spectral tail.
!>
!>  @brief Integrated pseudo-momentum over the diagnostic spectral tail.
!>
!>  @details
!>  The integrated momentum of the spectral tail has a form
!>  @verbatim
!>    M(ft) = GRAV * Dsn fp^-2 Im(nt), nt = log(ft/fp)/log(DFH)
!>  @endverbatim
!>  where ft it the tail start frequency, Dsn is the prognostic spectral level
!>  at ft, and fp is a 'peak frequency' value that makes the tail extension
!>  match the prognostic spectrum at ft.
!>
!>  M is a first 'kv/sig'-moment - like the the first spectral moment
!>  (that determines the mean wave period) but considered a vector sum.
!>  kv is the wave number vector, kv = k * (cos(theta), sin(theta))
!>
!>  Assuming a Donelan-Ewans spectral tail form, Dsn is estimated from the 
!>  omnidirectional spectrum S(sig) as Dsn = S(sig)*sig^4.
!>          
!>  A common look-up table for Im(nt) is calculated in a subroutine specfic
!>  for the parametric tail (for Donelan-Ewans: subroutine diag_DE_init(NP)).
!>  The index nt corresponds to ft/fp rounded to nearest integer:   
!>  @verbatim
!>    ft/fp = DFH**ntf; nt=nint( ntf )
!>  @endverbatim
!>        
!>
!>  @brief Diagnostic tail until Ncut.
!>
!>  @details
!>  The diagnostic tail contribution is calculated in the subroutine
!>  CALC_STVP(A) using a common 2-D array SBDiaZ(:,:) for the
!>  normalized Stokes drift.
!>
!>  The look-up array start index nt, and the spectral level Dsn, are
!>  defined as above
!>
!>  @verbatim
!>  StkDiag(:) = SBDiaZ(:,nt) * (1. - ntf) ! Fraction at fd/fp
!>  do ntt = nt+1, min(ntr,NP)
!>    StkDiag(:) = StkDiag(:) + SBDiaZ(:,ntt)
!>  end do
!>  @endverbatim
!>
!>  Starting at the surface, SBDiaZ(0,:) = DStkOsn(:). This is calculated
!>  in a subroutine specfic for the parametric tail. This example is for
!>  Donelan-Ewans, subroutine diag_DE_init(NP):
!>
!>  @verbatim
!>  fh = fht
!>  do ntt = Np+1, Ncut
!>        DStkOsn(ntt) = exp( aa * sinh(phih * (fh - 1.)) )
!>        Sum = Sum + DStkOsn(ntt)
!>        fh = fh * DFH
!>  end do
!>  DStkOsn(Np+1)=Sum
!>  @endverbatim
!>
!>  The depth dependent array SBDiaZ(1:NZ,:) is calculated in the
!>  subroutine STVP_INIT:
!>
!>  Let XKH = DFH**2
!>
!>  Initial value of ratio of wavenumbers:
!>  @verbatim
!>  Rk = fOfp(JSEA)**2
!>  ! Loop over spectral bins
!>  do ntt = NP + 1, Ncut
!>    ! Loop over depth bins
!>    do iz = 1,NZ
!>      ! Assume no shallow water effect at high frequencies
!>      A2S = EXP( - 2 * Rk * ZK_S(iz) )
!>      SBDiaZ(iz,ntt) =  DStkOsn(ntt) * A2S
!>      end do
!>    ! Increment the ratio of wavenumbers
!>    Rk = Rk * XKH
!>    end do
!>  @endverbatim
!>
!>  @Brief Implementation in the WW3 code.
!>
!>  @details
!>  Code lines in WW3 cource files other than w3stvpmd.ftn have all been
!>  added under the compile switch 'STVP'
!>
!>  a) In w3iogomd.F90:
!>       CALC_STVP(A) is called from W3OUTG()
!>
!>  b) In w3gdatmd.F90, w3iogrmd.F90, w3gridmd.F90:
!>       Part of the GRIDS structure (W3GDATMD) and stored in mod_def
!>  @verbatim
!>       SPND: Number of depths for Stokes profile U_S(1:SPND)
!>       SPDS: Depth scale specifying the largest depth Z(SPND) for the profile
!>       SPBP: Power of profile depth progression
!>       USXT: Tail param. type. 'DoEw': Donelan-Ewans,
!>                               'none': Prognostic spectrum frequencies, only
!>  @endverbatim
!>
!>  c) In w3adatmd.F90:
!>       Declare and allocate arrays ZK_S, U_S, V_S together with the integral
!>       parameters M_X, M_Y, K_S, all to be written to the output packed in
!>       a 2D array USVP
!>
!>
!>  @Brief Subroutines and functions.
!>
!>  @details
!>  @verbatim
!>      Name      Type  Scope    Description
!>     ----------------------------------------------------------------
!>      CALC_STVP    Public  Interface routine for Stokes drift
!>                              profile calculation
!>     ----------------------------------------------------------------
!>  @endverbatim
!>

!/ ------------------------------------------------------------------- /
      USE CONSTANTS, ONLY: GRAV, TPI, RADE, TSTOUT
!/
      USE W3ADATMD, ONLY: CG, WN, DW

! Integral pseudo-momentum vector and wave number scale K_S
      USE W3ADATMD, ONLY: M_X, M_Y, K_S
! Full Stokes profile vectors
      USE W3ADATMD, ONLY: U_S, V_S, ZK_S

! Stokes drift parameters
! Part of the GRIDS structure (W3GDATMD) and stored in mod_def:
! SPND: Number of depths for Stokes profile U_S(1:SPND)
! SPDS: Depth scale specifying the largest depth Z(SPND) for the profile
! SPBP: Power of profile depth progression
! USXT: Tail parametric type. 'DoEw': Donelan-Ewans,
!                             'none': Prognostic spectrum frequencies, only
! USXF: Lowpass cutoff frequency for the diagnostic tail extension
      USE W3GDATMD, ONLY: SPND, SPDS, SPBP, NSEAL, DMIN, XFR, DTH,  &
                          SIG, DDEN, ECOS, ESIN, NK, NTH, USXT, USXF

! VERBOSENESS%STVP: Verbose level [0..4] of STVP output to NDSV. In WW3_shel.nml
      USE W3ODATMD, ONLY: NDST, NDSO, NDSE, IAPROC, NAPROC, NAPOUT, VERBOSENESS
      USE W3DISPMD, ONLY: DSIE, N1MAX, ECG1, EWN1

! ! Lowpass cutoff frequencies for specific output fields. These are namelist
! ! parameters for ww3_shel (and declared in W3ODATMD in order not to link
! ! w3nmlshelmd with the program ww3_ounf)
!       USE W3ODATMD, ONLY: OFCUT, OFCUT_COUNT

! TODO:
!   SPND, SPDS and SPBP may be set at run init from WW3_shel.nml

      PRIVATE

!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/

! Stokes drift for the diagnostic tail. Scalable array and depths factor
! DFH: Relative increment from frequency bin-to bin for the scalable
!      diagnostic Stokes array
! NP: Number of frequency bins in the scalable diagnostic stokes array.

      REAL               :: DFH, XKH
      INTEGER            :: NP

! File id. for verbose output
      INTEGER            :: NDSV
!
      INTEGER            :: stvp_verbose = 1

! Look-up lists for Stokes spectral tail calculation
      REAL, allocatable  :: DStkOsn(:), Im(:)
      REAL, allocatable  :: SBDiaZ(:,:) ! (NZ,NK)

      REAL, allocatable  :: fOfp(:)

! Spectral integral parameters in Ewans functions
      REAL               :: sa, sb, ta, tb, tc
      REAL               :: m1_max, m_deriv_min, fOfp_min, fOfpCut
! 1st circular moment fit to spectrum at diagnostic frequency
      REAL               :: m1E

! Spectral parameters for Stokes spectral tail calculation
      REAL               :: SBDia

! CHA 20200603  ftr: truncation frequency, default 2 Hz
!              ntre: relative truncation bin ntre = ntr - nt
      REAL               :: ftr = 2.0
      INTEGER            :: ntre

! Profile arrays
      REAL, allocatable  :: StkDiag(:)

      REAL, allocatable  :: A2S(:), TZ_S(:), ETKZ(:), DZK_S(:)

! Depth integrated Stokes drift. Equals the integral wave pseudo-momentum.
      REAL               :: MW_U, MW_V

! Tail transition indices
! IKST: Lower edge of tail estimation from prognostic spectrum
! IKT:  Upper edge of tail estimation.
      INTEGER            :: IKST, IKT

! For calculation of shallow water Stokes drift
      REAL               :: fkd2,fkdm2,kd,e2kd,shkd4
      REAL, allocatable  :: fkd(:,:), fkdm(:,:) ! (1:NK,1:NSEAL)

! 2* (WN(IK,JSEA)-WN(IK-1,JSEA))
      REAL, allocatable  :: DWN(:)

      PUBLIC :: CALC_STVP
!/
      CONTAINS


!/ ------------------------------------------------------------------- /
!> @brief STVP_INIT(): Initialisation.
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
!/    16-Aug-2011 :
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
      INTEGER                  :: JSEA, ISEA, IK
! Control integers
      INTEGER                  :: IERR

! Common Stokes drift depth distribution parameters
      REAL                     :: DEPTH, tRk, XZK, kd_t, kd_d
      INTEGER                  :: IZ, ntt, ikd
      logical                  :: OPENED

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
! Initialization for step b) CALC_STVP (). An interface to DSTOKES()

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

      allocate( StkDiag(SPND), stat=IERR )

      ! For calculating exp(2 K Z) = exp(2 K0 Z) * exp( 2 (K-K0) Z )
      allocate( DWN(NK), stat=IERR )


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

      ! Old hardcoded value: SPBP = 1.

      XZK = 2.**((SPND*1.-1.)**(-SPBP))
      DO IZ = 1,SPND
        ZK_S(IZ) = XZK**((IZ*1.-1.)**SPBP)
      END DO
      ZK_S(:) = SPDS*(ZK_S(:) - 1.)

      IF ( stvp_verbose .gt. 1 ) THEN
        ! DZK_S = K_S * approximate thicknesses (m) of each layer
        ! Used in control calculation of the integrated pseudo momentum

        allocate( DZK_S(SPND), stat=IERR )

        ! The uppermost layer at the surface has approximately half thickness.
        DZK_S(1) = ZK_S(1) * 0.5
        DZK_S(2:SPND-1) = (ZK_S(3:SPND) - ZK_S(1:SPND-2)) * (XZK-1./XZK) * 0.5
        ! A convenient thickness of the lowest layer is
        DZK_S(SPND) = (SPDS*(XZK**(SPND**SPBP)-1.)-ZK_S(SPND-1))    &
                                                          * (XZK-1./XZK)*0.5
      END IF

      ! Shallow water Stokes profile parameters fkd(1:NK,1:NSEAL),
      ! fkdm(1:NK,1:NSEAL) for calculation of the action-to-stokes factor
      ! A2S(:) in subroutine DSTOKES

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

          ! kd<kd_t: shallow water or transitional to deep water
          e2kd = exp (2.*kd)
          shkd4 = 4. * (sinh( kd ))**2
          fkd(IK,JSEA) = e2kd / shkd4
          fkdm(IK,JSEA) = 1./ (e2kd * shkd4)
          ! In subroutine DSTOKES the action-to-stokes factor is calculated as
          ! A2Sf[IK,iz] = fkd[IK,JSEA] * ETKZ + fkdm[IK,JSEA] / ETKZ, where
          ! ETKZ = EXP( - 2. * kz[iz] )
          ! Note: At the surface where kz[1] = 0 we have ETKZ = 1, and 
          !       A2Sf[IK,1] = 0.5 * COSH(2.*kd) / SINH(kd)**2
          !                  = fkd[IK,JSEA] + fkdm[IK,JSEA]
          
        END DO ! IK=1,NK

      END DO ! JSEA = 1,NSEAL

      ! If the spectral tail is truncated at a specifified frequency ftr,
      ! DO IK=1,OFCUT_COUNT%N_FIELD+1
      !   IF (IK == OFCUT_COUNT%N_FIELD+1) exit
      !   IF (OFCUT(IK)%FIELD == 'SVP') exit
      !   END DO
      ! IF (IK <= OFCUT_COUNT%N_FIELD) ftr = OFCUT(IK)%FREQ

      IF (USXF < 10.0) ftr = USXF ! Default 10.0 Hz understood as infinity
      
      ! IF ( stvp_verbose .gt. 0 ) &
      !     WRITE (NDSV, *), ' STVP_INIT: CUT-OFF AT',ftr,'Hz'

      ! IKT: The upper edge of the prognostic range.
      ! The default is a truncated spectrum with no tail extension
      IKT = NK-1
      IKST = NK
      ! IF ( ftr < SIG(NK)/TPI) THEN
      !     IF ( stvp_verbose .gt. 0 ) &
      !         WRITE (NDSV, *), ' CUT-OFF below NK - No tail extension'
      !     RETURN
      ! END IF
      
!/ ------------------------------------------------------------- /
! Initializations for calculations related to the diagnostic, tail range

! cha 20121123: Shift fitting interval limits away from NK, from NK-5 to NK-2
      IF ( USXT == 'DoEw' ) THEN
         ! Donelan-Ewans
         IKT = NK - MAX( 2, NINT( 0.2/log(XFR) ) )
! cha 20110627: IKST: Lower edge of tail fitting to the prognostic spectrum
         IKST = NK - MAX( 5, NINT( 0.5/log(XFR) ) )
      ELSE IF ( USXT == 'none' ) THEN
         RETURN
      ELSE
        WRITE (NDSV, 916), 'STVP_INIT: Unknown spectral tail type: ', USXT
        CALL EXTCDE ( 73 )
        END IF

      allocate ( fOfp(NSEAL), stat=IERR )
      ! A zero value will cause fd/fp to be initialized at the first timestep
      fOfp = 0.0

! DFH: Relative increment from frequency/wavenumber bin-to bin for the
!     scalable diagnostic stokes arrays
      DFH = XFR
      XKH = DFH**2
      NP = IKT

      ! The truncation frequency ftr corresponds to a bin number ntr, where
      ! ftr/fp=DFH**ntr, or ntr = log( ftr/fp ) / log(DFH).
      ! Relative to the transition frequency fd = SIG(IKT) / TPI we have
      ! ftr/fd = (ftr/fp) / (fd/fp). The tail truncation bin exceeds the
      ! transition bin at fd by ntre bins, where ftr/fd = DFH**ntre, or:
      ntre = log( ftr * TPI / SIG(IKT) ) / log(DFH)
      
      ! SBDiaZ(z,bin): A 'normalized, look-up' vertical Stokes Drift profile
      ! shape for spectral bins of the diagnostic tail.
      ! Factorized as DStkOsn(bin) * exp( -tRk(bin) * ZK_S(z) )
      !
      ! For the Diagnostic range the surface Stokes drift for every spectral
      ! bin is multiplied with a depth factor EXP( -2 K Z), assuming deep water
      ! dispersion at the high frequencies. We apply a look-up array
      !  SBDiaZ(1:NZ,1:NP) = DStkOsn(1:NP) * exp ( -2*Rk(1:NP)*ZK_S(1:NZ) )
      ! where Rk is a discrete array representing (f/fp)**2

      ! Let tRk = 2*Rk = 2 * XKH**ntt.

      ! NP  should span the full prognostic range, from F1 to
      ! F1 * DFH**NP == F1 * XFR**IKT or NP = ceiling(IKT * log(XFR)/log(DFH))
      ! If DFH==XFR, then NP=IKT
      ! Index shift to the next bin in DStkOsn: log(XFR) / log(DFH)
      ! NP = ceiling(IKT * log(XFR) / log(DFH))

      IF ( .not. allocated(SBDiaZ) ) allocate( SBDiaZ(SPND, NP), stat=IERR )

      allocate ( DStkOsn(NP), stat=IERR )

      ! Im(bin): A look-up table for integrated momemtum
      IF ( .not. allocated(Im) )  allocate ( Im(NP+1), stat=IERR )

      ! Stokes drift for the diagnostic tail. Calculate a normalized array
      ! DStkOsn(1:NP)
      ! Also calculate the look-up table Im(1:NP)
      IF ( stvp_verbose .gt. 1 ) &
        WRITE (NDSV, 904), 'call diag_DE_init(', NP, ')'
      IF ( USXT == 'DoEw' ) THEN
        ! Donelan-Ewans
        CALL diag_DE_init(NP)
        END IF

      ! Extend with a zero to be used for linear interpolation
      Im(NP+1) = 0.

      tRk = 2.
      DO ntt=1,NP
        tRk = tRk*XKH ! = 2 * XKH**ntt
        SBDiaZ(:,ntt) = DStkOsn(ntt) * exp( -tRk * ZK_S(:) )
        END DO

      IF ( stvp_verbose .gt. 2 ) WRITE (NDSV, *),                   &
              'NP =',NP,'2Rk(NP) =', tRk, 'ZK_S(:) =', ZK_S

      IF ( stvp_verbose .gt. 2 ) WRITE (NDSV, *),                   &
              '   SBDiaZ(1:SPND,int(NP/2)) =', SBDiaZ(1:SPND, int(NP/2))

      deallocate ( DStkOsn )

  904 FORMAT ('  ', A,I5,A)
  906 FORMAT ('  ', 6(A,G9.2))
  916 FORMAT ('  ', A, A)

      END SUBROUTINE STVP_INIT


!/ ------------------------------------------------------------------- /
!> @brief CALC_STVP (A) - an interface to DSTOKES().
!>      
!> @details
!>  At each WW3 output time step, calculation of a Stokes drift profile is
!>  provided in arrays
!> @verbatim
!>                 U_S(1:SPND), V_S(1:SPND)
!> @endverbatim
!>  The Stokes profile is a sum of a prognostic part and a tail part
!>  StkDiag(1:SPND), where the depths are given in the array Z_S(1:SPND)
!/ ------------------------------------------------------------------- /

      SUBROUTINE CALC_STVP (A)
!/
!/            +------------------------------------------+
!/            | Extended spectral tail for Stokes drift  |
!/            |          Carsten Hansen                  |
!/            | Danish defence Joint GEOMETOC Support    |
!/            |                        FORTRAN 90        |
!/            | Last update :              22-Dec-2017   |
!/            +------------------------------------------+
!/
!/
!/    10-Oct-2011: Interface routine for Stokes drift profile calculation
!/                 To be called from W3OUTG() *after* calculation of both
!/                 integral parameters and swell partitions (W3CPRT)
!/                 The code structure is modelled after w3iosfmd.ftn
!/    22-Dec-2017 : Adoption to version WW3 v5.16


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

!     Generate data arrays for Stokes drift calculation
!     including look-up tables for extended tail
!
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
! 3.  Perform Stokes drift calculation
        CALL DSTOKES ( A(:,:,JSEA), JSEA, T02(JSEA) )
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

      SUBROUTINE DSTOKES(A, JSEA, T02)

      USE W3ADATMD, ONLY: U10
      USE W3GDATMD, ONLY: MAPSF ! Present module scope: , NSEAL, NK, NTH, SIG

      USE W3SERVMD, ONLY : EXTCDE
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
      REAL                  :: Ts, DEPTH, KDPT
      ! 1st circular moment fit
      REAL                  :: m1t
      ! Tail mean direction and spectral level parameter
      REAL                  :: thetat, Dsn

      INTEGER               :: numDepths, IZ, nt, ntt, ntr

      ! Factors for the diagnostic (spectral tail) calculation
      REAL                  :: Rk, fOfp_ini, ntf

      ! Factors for the prognostic calculation
      REAL                  :: A2S0, AX, AY, EDWTZ, DDWTZ

! cha 20110627: Tail transition indices
      REAL                  :: FACDIA

      INTEGER               :: ISEA, ITH, IK, NIK
      REAL                  :: CTH, STH
      REAL                  :: bandDia

      ! Status integer
      INTEGER               :: II

! Integral pseudo-momentum parameters
      REAL                  :: A2M, Mdiag, MprogX, MprogY, MtransX, MtransY

! Stokes drift wave number scale and parameters for 'inlined version of WAVNU1'
      REAL                  :: SIX, R1
      INTEGER               :: I1

#ifdef W3_DIST
      ISEA       = IAPROC + (JSEA-1)*NAPROC
#endif
#ifdef W3_SHRD
      ISEA       = JSEA
#endif

      ! Ts is an integral wave period chosen as Ts==T02
      Ts = T02

      CTH=0.
      STH=0.
      MprogX = 0.
      MprogY = 0.
      MtransX = 0.
      MtransY = 0.
      Mdiag = 0.
      U_S = 0.
      V_S = 0.
      StkDiag = 0.


      IF ( stvp_verbose .gt. 0 .and. JSEA .eq. 1 ) &
          WRITE (NDSV, 907), '  Transition range IKST =', IKST, 'to IKT =', IKT

      DEPTH  = MAX ( DMIN, DW(ISEA) )

      ! Inlined version of WAVNU1.
      SIX    = TPI/Ts * SQRT(DEPTH)
      I1     = INT(SIX/DSIE)
      IF (I1 .LE. N1MAX) THEN
          R1 = SIX/DSIE - REAL(I1)
          K_S = ( (1. - R1)*EWN1(I1) + R1*EWN1(I1 + 1) ) / DEPTH
        ELSE
          K_S = (TPI/Ts)**2 / GRAV
        END IF

      ! From here, we will use the dimensionless water depth
      KDPT = K_S * DEPTH

      ! numDepths: Number of depths of the discrete profile above the sea floor
      numDepths = SPND ! SPND is the configured value

      DO IZ = SPND,1,-1
        IF ( ZK_S(IZ) > KDPT ) CYCLE
        numDepths = IZ
        exit
      END DO

      IF ( stvp_verbose .gt. 2 ) &
        WRITE (NDSV, 905), 'numDepths =', numDepths, 'XFR =', XFR, 'IKT =', IKT

      IF ( IKT < IKST ) THEN

        IF ( stvp_verbose .gt. 0  .and. JSEA .eq. 1 ) &
          WRITE (NDSV, 912), 'Truncate at NK-2. No diagnostic spectral extension'

      ELSE

        IF ( stvp_verbose .gt. 1  .and. JSEA .eq. 1 ) &
             WRITE (NDSV, 912), '    .. for the diagnostic part'

        ! At the prognostic range cut-off frequency, derive the mean direction
        ! thetat, the 1st circular moment m1t, and the spectral level Dsn
        call prog_edge_DE( A, thetat, m1t, Dsn, JSEA)

        CTH=cos(thetat)
        STH=sin(thetat)

        ! Initialize the module parameter fOfp, if needed
        IF ( fOfp(JSEA) .eq. 0.0 ) THEN
          ! Initial guess of equivalent peak wave period fp = 1.2* 1/Tz
          fOfp_ini = max( 1.2 * Ts * SIG(IKT) / TPI, fOfp_min )
          fOfp(JSEA) = fOfp_ini
          IF ( stvp_verbose .gt. 2 ) &
               WRITE (NDSV, 900), '-> initial estimate, fd/fp =', fOfp_ini
        ELSE
          fOfp_ini = fOfp(JSEA)
          IF ( stvp_verbose .gt. 2 ) WRITE (NDSV, 900),             &
               '. Use Last estimate as initial for fd/fp =', fOfp_ini
        END IF

        IF ( m1t .ne. 0 ) THEN
          ! Calculate the equivalent value of fd/fp in stokes drift tail
          ! ( using Ewans' formulas )
          call equivalent_tail(m1t, fOfp(JSEA), stat=II)

          IF ( II > 0 ) THEN
            ! New guess of equivalent peak wave period fd/fp

            fOfp(JSEA) = max( 1.2 * Ts * SIG(IKT) / TPI, fOfp_min )
            IF ( fOfp(JSEA) /= fOfp_ini ) THEN
               fOfp_ini = fOfp(JSEA)
               call equivalent_tail(m1t, fOfp(JSEA), stat=II)
            END IF
            IF ( II > 0 ) THEN
               fOfp(JSEA) = fOfp_ini + 0.0625 * ( 20.0 - fOfp_ini )
               IF ( stvp_verbose .gt. 0 ) THEN
                  WRITE (NDSV, 900),                                &
                       'Issue in w3stvpmd, equivalent_tail() for Ts = ', Ts
                  WRITE (NDSV, 908), 'status = ', II,               &
                       '  -> Use an enhanced estimate, fd/fp =', fOfp(JSEA)
                  WRITE (NDSV, 907), ' ISEA = ', ISEA,              &
                       'at I =', MAPSF(ISEA,1), ', J =',MAPSF(ISEA,2)
               END IF
               call equivalent_tail(m1t, fOfp(JSEA), stat=II)
            END IF
          END IF
          IF ( II > 0 ) THEN
            fOfp(JSEA) = fOfp_ini
            IF ( stvp_verbose .gt. 0 ) THEN
               WRITE (NDSV, 908), &
                    'Still a Stokes issue, status = ', II,          &
                    '. Use a parametric fd/fp =',  fOfp(JSEA)
            END IF
          END IF
        END IF ! ( m1t .ne. 0 )

        fOfp(JSEA) = min( max ( fOfp(JSEA), fOfp_min ), DFH**NP)

        ! Index nt in SBDiaZ(:,nt) corresponds to fd/fp rounded to nearest
        ! integer: fd/fp = DFH**ntf; nt=nint( ntf )
        ! The part of the pre-calculated tail that covers the present tail
        ! is the range ntt = nt,NP
        ntf = log( fOfp(JSEA) ) / log(DFH)
        ntf = max( min(ntf,real(NP)), 1. )
        nt = nint( ntf )

        ! The high-frequency tail truncation bin ntr exceeds the
        ! transition bin ntf at SIG(IKT) by ntre bins
        ntr = nint( ntf + ntre )
        ! There is a forced truncation at the highest pre-calculated tail
        ! bin NP
        ntr = min(ntr,NP)

        ! Now, let ntf be the fractional part of nt
        ntf = ntf - nt

        ! The integrated pseudo-momentum of the diagnostic spectral tail is
        ! Mdiag = Dsn * sigP^{-2} Im(nt), where nt = log(ft/fP)/log(DFH)
        ! and sigP = SIG(IKT) / fOfP.
        ! Linear interpolation to nt+ntf between bins nt and nt+1. Note
        ! that if nt == NP, then ntf == 0., and note that Im(NP+1) == 0.
        Mdiag = Dsn * ( fOfP(JSEA)/SIG(IKT) )**2                    &
                * ( Im(nt) * ( 1. - ntf ) + ntf * Im(nt+1) - Im(ntr+1) )

        ! Calculate a vertical profile for the Stokes Drift. The surface Stokes
        ! drift for every spectral bin is multiplied with a depth factor
        ! EXP( -2 K Z), assuming deep water dispersion relation at the high
        ! frequencies.

        ! We apply a look-up array at bin nt (see SUBROUTINE diag_DE_init()):
        !  SBDiaZ(iz,nt) = DStkOsn(nt) *  exp ( -2*Rk(nt) * ZK_S(iz) )
        ! where  fOfp(JSEA)**2 is in the interval [Rk(nt), Rk(nt+1)[
        ! The look-up array is initialized in the procedure STVP_INIT()

        ! The normalized ( tail Stokes drift profile ) / Dsn is

        ! -for the fractional interval nt to nt+ntf:
        StkDiag(:) = SBDiaZ(:,nt) * (1. - ntf)

        ! -for the rest of the precalculated tail
        DO ntt = nt+1, ntr
          StkDiag(:) = StkDiag(:) + SBDiaZ(:,ntt)
        END DO

        IF ( stvp_verbose .gt. 1 ) THEN
          WRITE (NDSV, 905), 'ISEA =', ISEA, ' Dsn/U10 =', Dsn/U10(JSEA), &
               ' U10 =', U10(JSEA), ' Ts =', Ts
          WRITE (NDSV, 905), &
               ' ft/fP =', fOfp(JSEA), ' StkDiag0 =', StkDiag(1)*Dsn, &
               ' Mdiag = ', Mdiag
        END IF

        ! Multiply with the normalization factor (spectral level * U10)
        StkDiag(:) = StkDiag(:) * Dsn

        ! The Stokes drift values are zero below the sea floor
        StkDiag(numDepths+1:SPND) = 0.

      END IF ! ( IKT < IKST )

      ! Prognostic spectrum range IK=1, IKST-1

      IF ( stvp_verbose .gt. 1 ) THEN
        WRITE (NDSV, 912), '    .. for the prognostic part'
      END IF

      ! Stokes profile exponential shape for IK=1:
      ETKZ(:) = 1

      MprogX = 0.
      MprogY = 0.

      ! From STVP_INIT:
      ! For calculating exp(-2 K Z) = exp(-2 K0 Z) * exp( -2 (K-K0) Z ),
      ! DWN(IK) = K - K0 = WN(IK) - WN(IK-1)
      DWN(1) = WN(1,ISEA)
      DWN(2:NK) = WN(2:NK,ISEA) - WN(1:NK-1,ISEA)

      TZ_S(:) = -2.0 * ZK_S(:)/K_S

      DO IK=1, IKST-1

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

        ETKZ(:) = ETKZ(:) * EXP( DWN(IK) * TZ_S(:) )

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
        MprogX = MprogX + A2M * AX
        MprogY = MprogY + A2M * AY

      END DO

      ! Smoothen the transition from the prognostic to the diagnostic range

      NIK = IKT - IKST + 1
      ! Existance of at least one transitional bin means there is a tail
      ! Note, the highest bin has width 0.5
      IF (stvp_verbose .gt. 1 .and. NIK > 0) &
        WRITE (NDSV, 907), '    .. for the transitional tail',IKST, 'to', IKT

      DO IK=IKST, IKT

        ! ( We silently ignore a possible case that the peak is above IKST
        ! This has little influence when IKST is large enough that the Stokes
        ! drift becomes small anyway.)

        ! Action-to-momentum factor, times angular frequency.
        ! M = A k, also in shallow water.
        A2M = WN(IK,ISEA) * DDEN(IK) / CG(IK,ISEA)
        ! Take only half of the highest bin
        IF ( IK == IKT ) A2M = A2M * ( sqrt(XFR) - 1 )/( XFR - 1 )

        ! Action-to-Stokes factor at surface
        A2S0 = 2.0 * A2M * SIG(IK)
        ! Action-to-momentum factor
        A2M = A2M / SIG(IK)

        ! Action in spectral band projected on x and y
        AX = 0.
        AY = 0.

        DO ITH=1, NTH
          AX = AX + A(ITH,IK) * ECOS(ITH)
          AY = AY + A(ITH,IK) * ESIN(ITH)
        END DO

        ! Surface diagnostic tail Stokes drift transformed to wave action
        bandDia = SBDiaZ( 1, nt - (IKT-IK) ) * Dsn / A2S0

        ! Apply a nearly linear interpolation between prognostic and
        ! diagnostic limits to damp influence of spurious peaks near IKT.
        ! Assuming surface Stokes drift in the band is nearly a constant as
        ! function of the frequency

        FACDIA = (float(IK - IKST) + 0.5)/NIK
        AX = AX * (1.-FACDIA) + bandDia * CTH * FACDIA
        AY = AY * (1.-FACDIA) + bandDia * STH * FACDIA

        ! Stokes profile exponential shape:
        ! ETKZ(:) = EXP( - 2. * WN * Z_S(:) )

        ETKZ(:) = ETKZ(:) * EXP( DWN(IK) * TZ_S(:) )

        ! Action-to-Stokes profile factor
        A2S(:) = A2S0 * ETKZ(:)

        ! Stokes drift projected on x and y
        U_S(:) = U_S(:) + AX * A2S(:)
        V_S(:) = V_S(:) + AY * A2S(:)

        ! Integrated pseudo-momentum
        MtransX = MtransX + A2M * AX
        MtransY = MtransY + A2M * AY

      END DO !IK=IKST, IKT

      IF ( IKT >= IKST ) THEN
        ! Add the contributions from the tail, projected on x and y
        U_S(:) = U_S(:) + StkDiag(:) * CTH
        V_S(:) = V_S(:) + StkDiag(:) * STH
        MprogX = MprogX + MtransX + Mdiag * CTH
        MprogY = MprogY + MtransY + Mdiag * STH
      END IF

      ! Let the drift be zero below the sea floor
      U_S(numDepths+1:SPND) = 0.
      V_S(numDepths+1:SPND) = 0.

      ! The integral pseudo-momentum per unit surface area, divided by density
      M_X = MprogX * GRAV
      M_Y = MprogY * GRAV

      ! Output if high level of verboseness

      IF ( stvp_verbose .le. 1 ) RETURN

      WRITE (NDSV, 912), '    .. done'

      ! Compare M_X, M_Y with the depth-integrated Stokes drift.
      ! Depth-integration is approximated as a sum over Z_S(1:numDepths).
      !
      MW_U = SUM ( U_S(1:numDepths) * DZK_S(1:numDepths)/K_S )
      MW_V = SUM ( V_S(1:numDepths) * DZK_S(1:numDepths)/K_S )

      ! Depth of lower edge of the layer at numDepths
      IF ( ZK_S(SPND) > KDPT ) THEN
         MW_U = MW_U - SUM ( U_S(1:numDepths) * ( DZK_S(numDepths)  &
                             - (0.5*(ZK_S(IZ-1) + ZK_S(IZ)) - KDPT) )/K_S )
         MW_V = MW_V - SUM ( V_S(1:numDepths) * ( DZK_S(numDepths)  &
                             - (0.5*(ZK_S(IZ-1) + ZK_S(IZ)) - KDPT) )/K_S )
      END IF
      ! Linear interpolation in the uppermost layer
      MW_U = MW_U - 0.25 * ( U_S(1) - U_S(2) ) * DZK_S(1)/K_S
      MW_V = MW_V - 0.25 * ( V_S(1) - V_S(2) ) * DZK_S(1)/K_S

      WRITE (NDSV, 914), ' MW_U, MW_V  = ', MW_U,  MW_V
      WRITE (NDSV, 914), ' M_X, M_Y    = ', M_X,  M_Y
      IF ( stvp_verbose .gt. 2 ) THEN
        Mdiag = Mdiag  * GRAV
        WRITE (NDSV, 914), ' Mt_X, Mt_Y  = ', MtransX * GRAV, MtransY * GRAV
        WRITE (NDSV, 914), ' Md_X, Md_Y  = ', Mdiag * CTH, Mdiag * STH
        WRITE (NDSV, 914), 'MDi = ', &
             SUM ( StkDiag(1:numDepths) * DZK_S(1:numDepths)/K_S )  &
             - 0.25 * ( StkDiag(1) - StkDiag(2) ) * DZK_S(1)/K_S
      END IF

! Formats
!
  900 FORMAT ('  ',A,G9.2)
  902 FORMAT ('  ',A,I5)
  905 FORMAT ('  ',A,I6,6(A,G9.2))
  906 FORMAT ('  ',6(A,G9.2))
  907 FORMAT ('  ',3(A,I6))
  908 FORMAT ('  ',A,I2,A,G9.2)
  909 FORMAT ('  ',2(A,G9.2))
  910 FORMAT ('  ',A,100G9.2)
  912 FORMAT ('  ',A)
  914 FORMAT ('  ',A,4(F8.4,', ',F8.4))
  915 FORMAT ('  ',A,9F8.4)

      END SUBROUTINE DSTOKES


      SUBROUTINE diag_DE_init(NP)

      !  Stokes drift for the diagnostic tail of Donelan_Ewans spectrum.
      !  Set parameters in Ewans' formulas
      !  Pre-calculate for STVP_INIT and CALC_STVP:
      !   - A normalized array DStkOsn(1:NP). The array should be approximately
      !     as long as the prognostic tail. We have chosen NP=NK-2
      !   - A table Im(1:NP) of normalized integrated momemtum
      !
      !  The surface Stokes drift for every spectral bin is multiplied with
      !  the depth factor EXP( -2 K Z), assuming deep water dispersion relation
      !  at the high frequencies.

      !  Calculation of the tail directional spreading:
      !
      !  Given a spectrum of shape
      !  S(sig) = H(sig,theta) * alpha * U * g * sig^{-4}
      !
      !  Let fh = f/fP denote the ratio of spectral bin frequency to the peak
      !
      !  At f/fP = 1., the contribution for the logarithmic bin of width log(df)
      !  to the surface Stokes drift is alpha * U times the factor
      !    SOsn = 2 * ( first circular moment ) * d f / f
      !  The first circular moment is defined as
      !    m1 = int_{-pi}^pi cos(theta) H(sig,theta) d theta
      !  Ewans (1997) provides an expression given as a function
      !  m1E(f/fP), based on at fit to observations (the 'Maui'
      !  experiment).
      !
      !  The binned contribution to the stokes drift is influenced by the first
      !  circular moment m1. A freq^{-4} spectral tail is assumed in the form
      !     d Stokes(f) = 2 alpha U m1 d f / f
      !
      !  With Ewans' expression for the directional spread:
      !     sigma = sa - sb * (fP/f)**2
      !     theta_m1 = ta * exp( tb * fP/f )
      !     m1E = cos(theta_m1) * exp( -(sigma**2)/2 )
      !
      !  In the vicinity of fh=f/fP ~ fh0 ~ 2.5, this may be approximated as
      !     m1t(fh) = m1E(fh=1) * exp( aa sinh(phih * (fh - 1.)) )
      !  where
      !     phih = (1/(aa m1) * {d m1 / d fh}
      !  m1t(fh) converges to zero at high frequences.
      !
      !  Choices of these parameters were compared graphically in a Python
      !  program and we end up with a good overall fit choosing
      !
      !  phih=-0.187
      !  aa = 1.25
      !

!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      IMPLICIT NONE

      INTEGER, intent(in)        :: NP

      ! local variables
      REAL                       :: aa,phih,logdf
      REAL                       :: fht, fh
      REAL                       :: dfhi
      REAL                       :: fhPm2, idf2, df2, ImS

      INTEGER                    :: iz, N, nt, ntt, Npcut

      ! The part of the pre-calculated tail that covers the present tail
      ! is the set of bins ntt >= nt

      ! Set parameters in Ewans' formulas
      sa = 32.13/RADE
      sb = 15.393/RADE
      ta = exp(5.453)/RADE / 2.
      tb = -2.750
      tc = 0.0845

      ! fOfpCut is an order of magnitude cut-off at the upper limit of
      ! validity for the directional distribution formula. The associated bin
      ! number is Npcut = NP + int(log(fOfpCut)/log(DFH))
      fOfpCut = 10.0

      ! Smallest fOfp valid for our approximations to Ewans' formulas
      fOfp_min = sqrt(sa/sb)

      ! A practical asymptotic value for small fOfp in m1_Ew98():
      m_deriv_min = m1deriv_Ew98( fOfp_min )

      ! First circular moment at fOfp_min
      m1_max = m1_Ew98( fOfp_min )

      logdf = log(DFH) ! Log of the increment frequency factor
                       ! Be aware: DFH may *differ* from the prognostic bin
                       ! increment factor XFR

      ! The bin number associated with fOfpCut is
      Npcut = NP + int( log(fOfpCut)/logdf )

      ! Inverted integration frequency increment
      dfhi = 1./DFH

      ! Suggested approximation, parameter values
      aa = 1.25
      phih = -0.187

      fh = 1.
      DO nt = 1,NP
        DStkOsn(nt) = exp( aa * sinh(phih * (fh - 1.)) )
        fh = fh * DFH
      END DO

      ! A discrete array representing d Stokes(f) / (alpha U) then becomes
      ! (with d f / f = logdf)

      DStkOsn = m1_Ew98(1.) * DStkOsn * 2. * logdf

      ! Lookup table Im(1:NP) (See '4. Integrated pseudo-momentun')

      fht = DFH**NP ! Value of fh just beyond the prognostic tail

      fh = fht
      fhPm2 = 1./fht**2
      idf2 = 1./DFH**2
      ImS = 0.
      DO ntt = NP + 1, Npcut
        fhPm2 = fhPm2 * idf2 != 1./(fh*DFH)**2 = fh ^ -2
        fh = fh * DFH
        ImS = ImS + exp( aa * sinh(phih * (fh - 1.)) ) * fhPm2
      END DO

      Im(NP) = ImS

      ! Step downwards the spectral frequencies of possible values of ft/fP
      fh = fht
      fhPm2 = 1./fht**2
      df2 = DFH**2
      DO nt = NP-1,1,-1
        Im(nt) = Im(nt+1) + exp( aa * sinh(phih * (fh - 1.)) ) * fhPm2
        fhPm2 = fhPm2 * df2 != 1./(fh/DFH)**2 = fh ^ -2
        fh = fh * dfhi
      END DO

      Im(NP+1) = 0. ! When used in extended linear interpolation
      Im(:) = Im(:) * ( m1_Ew98(1.) * logdf )

      IF ( stvp_verbose .gt. 2 ) WRITE (NDSV, *), 'Im(:) =', Im(:)

      IF ( stvp_verbose .gt. 2 ) &
        WRITE (NDSV, 906), 'Ewans parameters: sa =', sa,            &
               'sb =', sb, 'ta =', ta, 'tb =', tb,                  &
               'fOfp_min =', fOfp_min, 'm_deriv_min =', m_deriv_min

  904 FORMAT ('  ', A,I5,A)

  906 FORMAT ('  ', 6(A,G9.2))

      END SUBROUTINE diag_DE_init



!/
!/ ------------------------------------------------------------------- /
!/

      SUBROUTINE prog_edge_DE(A, theta, m1, alphaU, JSEA)

! Donelan-Ewans spectral tail
! At the prognostic edge where the tail is to be matched, estimate:
!        theta:   The mean direction
!        m1:      The 1st circular moment m1
!        alphaU: The spectral level (alpha * U)
!

!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      IMPLICIT NONE
      REAL, INTENT(IN)        :: A(:,:)   ! Wave action
      REAL, INTENT(OUT)       :: theta, m1, alphaU
      INTEGER, INTENT(IN)     :: JSEA

      REAL, allocatable, save :: DTHsig5_CgG(:,:)
      REAL, allocatable       :: DTHsig5_G(:)

!     Integral band sums
      REAL                    :: ABst, ABsX, ABsY

      REAL                    :: a1, b1
      REAL                    :: SUM, SUMX, SUMY
      INTEGER                 :: KSEA, ISEA, ITH, IK, IERR

! Spectral level of the tail, assuming S(sig) = alphaU * GRAV * sig^-4:
! alphaU = band_mean( \sum A(IK) * FACTOR*SIG(IK)*4/GRAV ), where
! FACTOR = SIG(IK) * DTH / Cg(IK,JSEA). This is similar to the subroutine
! CALC_U3STOKES in w3iogomd, but here we assume the frequency band contain
! deep water waves, thus WN = SIG^2/GRAV.

! Spectral factor for calculating the spectral level of the tail
      IF ( .not. allocated(DTHsig5_CgG) ) THEN
        allocate ( DTHsig5_CgG(NK,NSEAL), stat=IERR )
        allocate ( DTHsig5_G(NK), stat=IERR )

        DO IK=1,NK
          DTHsig5_G(IK) = DTH * SIG(IK)**5/GRAV
        END DO

        DO KSEA=1,NSEAL
#ifdef W3_DIST
          ISEA       = IAPROC + (KSEA-1)*NAPROC
#endif
#ifdef W3_SHRD
          ISEA       = KSEA
#endif
          DO IK=1,NK
            DTHsig5_CgG(IK,KSEA) = DTHsig5_G(IK) / CG(IK,ISEA)
          END DO

        END DO

          deallocate ( DTHsig5_G )
      END IF

! Integrate over wave action in bands (As in subroutine CALC_U3STOKES,
! in w3iogomd)
      SUM=0.
      SUMX=0.
      SUMY=0.
      
      ! Spectral level and circular moment from the average over an
      ! interval [IKST,IKT] of bands near NK
      
      DO IK = IKST,IKT
        ! Initialize wave action in the band
        ABst    = 0.
        ABsX    = 0.
        ABsY    = 0.
        ! Integration over wave action in band
        DO ITH=1, NTH        
          ABst  = ABst + A(ITH,IK)
          ABsX  = ABsX + A(ITH,IK)*ECOS(ITH)
          ABsY  = ABsY + A(ITH,IK)*ESIN(ITH)
        END DO
        ! Spectral level of the tail, assuming S(sig) = alphaU * GRAV * sig^{-4}
        ! alphaU = block_mean( ABst ):
        SUM  = SUM  + ABst * DTHsig5_CgG(IK,JSEA)
        SUMX = SUMX + ABsX * DTHsig5_CgG(IK,JSEA)
        SUMY = SUMY + ABsY * DTHsig5_CgG(IK,JSEA)
      END DO
        
      alphaU = SUM / (IKT-IKST+1)      

      theta = ATAN2( SUMY, SUMX )
      
      ! Spectral circular moment
      a1 =  SUMX/SUM
      b1 =  SUMY/SUM

      m1 = sqrt( a1**2 + b1**2 )

! cha first notes:
!    Wave-number and frequency energy spectrum are related as
!    (User manual for WW3, MMAB_276.pdf, Eq. 2.4):
!    F(k,th) dth dk = S(\sigma,th) dth d\sigma
!    Thus, wave-number action spectrum, A(k) = F(k)/\sigma, is:
!    A(k,th) dth dk = S(\sigma,th) / \sigma dth d\sigma
!
!    Komen et al., 1994 (the WAM Book) has a different definition
!    of spectral density in the wavenumber domain (Eq. 1.137):
!    F(k,th) k dth dk = S(\sigma,th) dth d\sigma
!
  906 FORMAT ('  ',6(A,G9.2))

      END SUBROUTINE prog_edge_DE

!/
!/ ------------------------------------------------------------------- /
!/

      SUBROUTINE equivalent_tail(m1,fOfpE,fOfp0,fOfp1,xacc,stat)
      ! Given a first circular moment m1, determine an equivalent value of
      ! the spectral peak, f/fpE, that reproduces the Maui spectral shape
      ! (see function m1_Ew98()).

      IMPLICIT NONE

      REAL, intent(in)           :: m1
      REAL, intent(inout)        :: fOfpE
      INTEGER, intent(inout), optional  :: stat
      REAL, intent(in), optional :: fOfp0, fOfp1, xacc

      REAL                       ::  fOfp0_, fOfp1_, xacc_, fOfpE_
      INTEGER                    ::  status_

      ! Default values
      fOfp0_=1.0
      ! fOfp1_=10.5751
      fOfp1_= 30.
      xacc_=0.01
      ! Substitute with optionals
      IF ( present(fOfp0)) fOfp0_=fOfp0
      IF ( present(fOfp1)) fOfp1_=fOfp1
      IF ( present(xacc)) xacc_=xacc

      ! Determine the value of the equivalent peak (fp) from a root fOfpE of
      ! m1Ed_Ew98(fOfp) = ( m1_Ew98( fOfp ) - m1 ). First set a global
      ! representation of m1 as a limited value m1E:

      IF ( m1 > m1_max ) THEN
        m1E = m1_max
      ELSE
        m1E = m1
      END IF
        
      fOfpE_=fOfpE
      ! Determine the root fOfpE
      fOfpE = rtnewt(m1Ed_Ew98, fOfpE_, 0.0, fOfp1_, xacc_, status_)

      ! Validate the result
      IF ( (fOfp0_-fOfpE)*(fOfpE-fOfp1_).lt.0.) THEN
        IF (stvp_verbose .gt. 0) &
          WRITE (NDSV, 908), 'WARNING: fOfpE = rtnewt() =', fOfpE,  &
            '. Result beyond valid range [',fOfp0_,',',fOfp1_,'] . m1=', m1
        status_ = 3
        fOfpE = fOfpE_
      END IF

      IF (present(stat)) stat = status_

  908 FORMAT ('  ',8(A,G9.2))

      END SUBROUTINE equivalent_tail

!/ ------------------------------------------------------------------- /
!/
!     Application of Numerical Recipes Par. 9.4:
!     Root of function by Newton-Raphson Method Using Derivative
!
      REAL FUNCTION rtnewt(funcd,rtnewtv,x1,x2,xacc,status)

      IMPLICIT NONE

      INTEGER                 :: JMAX
      REAL, INTENT(IN)        :: rtnewtv,x1,x2,xacc
      INTEGER, INTENT(OUT)    :: status
      EXTERNAL funcd
      PARAMETER (JMAX=20) ! Set to maximum number of iterations.
      ! Using the Newton-Raphson method, find the root of a function known
      ! to lie in the interval [x1, x2]. The root rtnewt will be refined
      ! until its accuracy is known within +/-xacc. funcd is a user-supplied
      ! subroutine that returns both the function value and the first
      ! derivative of the function at the point x.
      INTEGER                 :: j
      REAL                    :: df,dx,f
      status = 0

      IF (stvp_verbose .gt. 2) &
        WRITE (NDSV, 906), '  m1E = ', m1E, 'fOfpf =', rtnewtv

      rtnewt = rtnewtv !Initial guess.

      DO j = 1, JMAX
        call funcd(rtnewt,f,df)
        dx=f/df
        rtnewt=rtnewt-dx

        ! Remedy initial over-/undershoot
        IF (rtnewt < x1 .and. j < 3) rtnewt = x1 + 0.1*(x2-x1) * j
        IF (rtnewt > x2 .and. j < 3) rtnewt = x2 - 0.1*(x2-x1) * j
        
        IF ((x1-rtnewt)*(rtnewt-x2).lt.0.) THEN
          IF (stvp_verbose .gt. 1) &
            WRITE (NDSV, 914), '  WARNING in w3stvpmd.ftn: rtnewt =', rtnewt, &
                               'jumped out of brackets'
          IF (stvp_verbose .gt. 2) THEN
            WRITE (NDSV, 906), '  [',x1,',',x2,'] . f=', f, 'df=',df
            WRITE (NDSV, 906), '  j=', j, 'rtnewtv=', rtnewtv, 'm1E=', m1E
            END IF
          status=2
          RETURN
        END IF

        IF ( abs(dx) .lt. xacc) RETURN ! Convergence.
      END DO

      status=1
      rtnewt = rtnewtv
        IF (stvp_verbose .gt. 1) &
          WRITE (NDSV, 912), '  WARNING in w3stvpmd.ftn: rtnewt exceeded', &
                             ' maximum iterations'


  906 FORMAT ('  ',6(A,G9.2))
  912 FORMAT ('  ',6A)
  914 FORMAT ('  ',A,G9.2,A)

      END FUNCTION rtnewt

!/ ------------------------------------------------------------------- /
!/
!     Subroutine to determine with rtnewt() at which value of f/fp the first
!     directional moment m1 of Ewans 1998 equals the value at the tail cut-off
!
      SUBROUTINE m1Ed_Ew98(fOfp, f, df)

      IMPLICIT NONE
      REAL, INTENT(IN)        :: fOfp
      REAL, INTENT(OUT)       :: f, df

      f = m1_Ew98(fOfp) - m1E
      df = m1deriv_Ew98(fOfp)

      END SUBROUTINE m1Ed_Ew98


!/ ------------------------------------------------------------------- /
!/
      REAL FUNCTION m1_Ew98(fOfp)

      REAL, intent(in)        :: fOfp
      REAL                    :: m1, fpOf, sigma, theta_m1, fpD

      IF (fOfp .lt. fOfp_min) THEN
        fpOf = 1./fOfp_min
        fpD = fOfp - fOfp_min
      ELSE
        fpOf=1./fOfp
        fpD = 0.0
      END IF

      sigma = sa - sb * fpOf**2

      theta_m1 = ta * exp( tb * fpOf )

      m1_Ew98 = cos(theta_m1) * exp( -(sigma**2)/2 )

      ! Define a practical asymptotic value for small fOfp:
      IF ( fpD .lt. 0.0 ) THEN
        m1_Ew98 = m1_Ew98 + fpD * m_deriv_min
      END IF

      END FUNCTION m1_Ew98


!/ ------------------------------------------------------------------- /
!/
      REAL FUNCTION m1deriv_Ew98(fOfp)

      REAL, intent(in)        :: fOfp
      REAL                    :: fpOf, sigma, theta_m1

      IF (fOfp .lt. fOfp_min) THEN
        m1deriv_Ew98 = m_deriv_min
        RETURN
      END IF

      fpOf = 1./fOfp
      sigma = sa - sb * fpOf**2

      theta_m1 = ta * exp( tb * fpOf )

      m1deriv_Ew98 = ( sin(theta_m1) * theta_m1 * tb                &
                        - cos(theta_m1) * sigma * 2. * fpOf * sb    &
                     ) * exp( -(sigma**2)/2 ) * fpOf**2


      END FUNCTION m1deriv_Ew98


!/
!/ End of W3STVPMD ----------------------------------------------------- /
!/
      END MODULE W3STVPMD


