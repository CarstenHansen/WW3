!> @file
!> @brief Contains module W3XSMFMD.
!> 
!> @author C. Hansen  @date 28-Apr-2021
!>
#include "w3macros.h"
!/ ------------------------------------------------------------------- /
      MODULE W3XSMFMD
!/
!/                  +---------------------------------------------+
!/                  | Mattsson-Hansen fit to Stokes drift profile |
!/                  |          Carsten Hansen                     |
!/                  | Joint GEOMETOC Support Centre Denmark       |
!/                  |                          FORTRAN 90         |
!/                  | Last update :           09-OCT-2019         |
!/                  +---------------------------------------------+
!/
!/    02-May-2011 : Origination. Development program w3stokes.ftn converted
!/                  from tested Python code
!/    27-Jun-2011 : First module version w3dstkmd.ftn. Output tested with
!/                  a stand-alone program w3stokes.ftn.Implementation plan.
!/    18-Aug-2011 : Calculate fit to Mattsson parametric profile (ver.0.12)
!/    10-Oct-2011 : Interface subroutine w3dstk() to ww3_shel (version 0.22)
!/    21-Oct-2011 : Detailed description of implementation into WW3 v.3.14
!/                  First compilation of ww3_shel in parallel mode
!/                  Test-scripts run_ww3test.sh,
!/                     $PERL5LIB/bin/ww3job_preproc (version 0.32)
!/    08-Nov-2011 : Full functioning module. Checked against patches from
!/                  the WW3 errata page, applied: DEPTH=max(DW,DMIN).
!/                  Un-physical spatial discontinuities identified: May
!/                  need to match the Stokes profiles to integral wave
!/                  pseudo-momentum/depth-integrated Stokes drift (ver. 0.32).
!/    23-Nov-2011 : Ajust the deep profile to match the depth-integrated
!/                  Stokes drift. Function AJUST_FIT() (version 0.33).
!/    07-Dec-2011 : Refined ajustment of Mattsson-Hansen profiles to integral
!/                  pseudo-momentum (version 0.41).
!/    20-Feb-2012 : The integrated normalized pseudo momentum of the JM-profiles
!/                  is now (ver. 0.48) calculated from the surface to the
!/                  actual depth, and based on a mathematical expression for
!/                  the integral (not a numerical integral).
!/    02-Aug-2012 : Fitting of deep parameter PMd skipped (ver. 0.52)
!/
!/    06-Aug-2012 : (Ver. 0.55). Introduced a new procedure in profile fitting,
!/                  subroutine MH_fit(), which fully replaces the former
!/                  two-layers fit procedure.
!/                  MH_fit() performs, in a few iterative steps, a
!/                  parametric fit alternating with ajustment to the
!/                  integrated pseudo-momentum.
!/
!/    22-Dec-2017 : Adoption to WW3 version 5.16
!/                  + Call STVP() from W3OUTG(), outside the loop over JSEA
!/
!/    19-Aug-2019 : Adoption to WW3 version 7.XX
!/                  + Mattsson-Hansen fitting is separated to module 'W3XSMFMD'
!/                  + Namelist group &STVP
!/                  + Clearifying the in-code documentation
!/                  + Dimensionless parameters in calculations, e.g. M -> Kz M
!/                  + Replace Gradshteyn and Ryzhik with Numerical Recipes
!/                    in function gi()
!/                  + Quality count output field
!/                  +
!/
!/    11-Sep-2020 : Apply a lower limiter to abs(M)
!/                  Disabled: AM_tune
!/                  Disabled: swap if 1. - IPp/IPd < 0                    
!/                          
!/    30-Sep-2020 : Lookup table for the reverse of IP_calc when PM==PMds
!/                          
!/    16-Oct-2020 : Limit Ud projected on the opposite direction of Up
!/
!/    28-Apr-2021 : List of dimnsionless depths focus stronger to the surface.
!/                  Re-introduce fit with PMd and avoid use of the reverse
!/                  of IP_calc - just match Ud.
!/
!  1. Purpose, and output data
!
!  A parametric representation of the Stokes drift profile vector U_S(Z) is
!  provided as a sum of two expressions (suggested by J. Mattsson and C.Hansen),
!        U_S(Z) = Up exp( -AM (Z Kz)^PM ) + Ud exp( -AMd (Z Kz)^PMd ) ,
!  where Up + Ud = U0 is the Stokes drift vector at the surface (Z=0) and
!  Z is pointing downwards.
!  A wave number Kz is related to the zero-upcrossing frequency by the
!  linear-wave dispersion relation, (2pi/Tz)^2 = g Kz tanh( Kz D). The
!  local depth D should also be provided in the output field for the bathymetry.
! 
!  The vertical integral of U_S(Z) equals the integral wave (pseudo) momentum
!  vector M. Further details on the MHfit parameters AM, PM, AMd, PMd are
!  found in Sect, 3, 'More on MHfit and suggestions for postprocessing'
!
!  U_S(Z) is determined as a best fit to the total model profile U_mod(Z)
!  constructed as a sum of a prognostic and a diagonostic (tail) contribution.
!  These profiles have been derived from the model spectra + spectral tail
!  using the module STVPMD and subroutine CALL_STVP(). NDEPTH discrete depths
!  have been given in an array Z_S(1:NDEPTH).
!
!  Also calculated are a scatter index of the profile fit, and a measure of
!  the convergence of the iteration procedure.
!
!  All calculations are performed in parallel and the results are gathered
!  to the process that writes the output  variable XSMH.
! 
!  From the array USVP (see w3stvpmd.ftn) at least the following 5 values are
!  stored to out_grd.ww3 to be used in conversion to NetCDF (ww3_ounf) .
!
!    USVP(JSEA,1)      = Kz     ! g Kz tanh( Kz D) = (2pi/T02)^2 [ 1/m ]
!
!    USVP(JSEA,2)      = M_X    ! (M_X, M_Y) is the integral pseudo-momentum
!    USVP(JSEA,3)      = M_Y    ! per unit surface area, divided by density.
!                               ! [ m^2/s ]        
!    USVP(JSEA,4)      = U0_X   ! Surface Stokes drift [ m/s ]
!    USVP(JSEA,4+SPND) = U0_Y 
!
!  Via the array USEROXSMH we store other seven fitting parameters
!
!    XSMH(JSEA,1)     = Ud_X   ! Surface drift of 'deep' partition
!    XSMH(JSEA,2)     = Ud_Y   ! Primary partition drift is Up = U0 - Ud
!    XSMH(JSEA,3)     = AM     ! Fit parameter, primary partition
!    XSMH(JSEA,4)     = PM     ! Fit parameter, primary partition
!    XSMH(JSEA,5)     = AMd    ! Fit parameter, deep partition
!    XSMH(JSEA,6)     = PMd    ! Fit parameter, deep part.
!
!    XSMH(JSEA,7)     = SIC  ! Quality count, see Stokes_MHfit()
!    !    IT  : Number of iterations in MH_fit
!    !    SI  : Scatter index  std(|U_S - U_mod| * dz)/|M|
!    !    C   : Convergence of iteration
!        
!  2. Usage
!
!  a) In ww3_grid.inp/namelists.nml, you may modify four namelist variables:
!     $ NDP: Number of depths (->NDEPTH)
!     $ DSC: Depth scale specifying ZK, the largest depth Z(NDP) * Kz
!     $ TYP: Tail type (C*4);
!     $      TYP='DoEw': Donelan-Ewans extension to infinite freq.,
!     $      TYP='DE20': Donelan-Ewans truncated at 2.0 Hz,
!     $      TYP='None': Prognostic spectrum truncated with no extension
!     &STVP NDP = 31, DSC = 3.0, TYP = 'DE20' /
!
!  b) In w3initmd.ftn: Call STVP_INIT() (module W3STVPMD):
!     Allocation of StkProg(1:NDEPTH,I), StkDiag(1:NDEPTH,I), Z_Stk(1:NDEPTH)
!
!  c) In w3iogomd.ftn: Call the subroutine STVP() (module W3STVPMD): At each
!     WW3 output time step, calculation of the Stokes drift profile
!
!  d) Call the subroutine Stokes_MHfit () (present module W3XSMFMD):
!     Approximation of the profile with a parametric expression
!        U_S(Z) = Up exp( -AM (Z Kz)^PM ) + Ud exp( -AMd (Z Kz)^PMd ), 
!     where Kz is a characteristic wavenumber corresponding to T02, and the sum
!     of vectors Up + Ud is the Stokes drift at the surface.
!     Further details on the MHfit parameters Up, Ud, AM, PM, AMd, PMd are found
!     in Sect. 3.
!
!  3. More on Stokes_MHfit and suggestions for postprocessing
!
!  In the fitting procedure it is assured that the vertical integral of the
!  parametric drift profile equals the pseudo-momentum (or total Stokes
!  transport vector). Two dimensionless profile shape parameters are calculated,
!     IPp = int_0_D Kz exp( -AM (Z Kz)^PM ) d Z
!     IPd = int_0_D Kz exp( -AMd (Z Kz)^PMd ) d Z
!  where D is the water depth, and Kz a scaling wavenumber related to the
!  zero-crossing period T02 by the dispersion relation for low amplitude waves,
!  (2*pi/Tz)**2 = GRAV * Kz * tanh( D * Kz). (GRAV = 9.81).
!  IP is calculated in a function IP_calc.
!
!  Inclusion of the pseudo-momentum vector in the output allows for a possible
!  post-processing with temporal-spatial smoothing of the MHfit parameters
!  AM, PM, AMd, PMd. After smooting, the surface drift should be partioned
!  into a sum of a 'primary' partition and a 'deep' partition, U0 = Up + Ud,
!  following the operations below (Eqs. 1-7). After smoothing it must be
!  assured that the 'deep' profile decays slower than the 'primary' profile,
!  so that the difference of integrated profile shapes IPd - IPp is
!  significantly positive. This can be achieved by ajusting AM towards
!  a value of IPd - IPp which differ only slightly from an original, smoothed
!  value.
!        
!  Output data to use for reconstruction of the Stokes profile are:
!  Seafloor depth, Kz, M = (MX,MY), U0 = U_S(0) = (U0X,U0Y), Ud = (UdX,UdY),
!  AM, PM, AMd, and PMd.
!
!  The pseudo-momentum vector, normalized by its depth scale, is
!        
!  (1)     MK = (MX,MY) * Kz
!        
!  The 2-dimensional surface Stokes drift vector
!        
!  (2)     U0 = (U0X, U0Y)
!
!  Factorize the normalized pseudo-momentum vector MK in an 'primary part' and
!  a 'deep part' as
!          MK = Up * IPp + Ud * IPd
!  where the surface Stokes drift is
!          U0 = Up + Ud
!
!  Here, IPp and IPd are dimensionless integrated profiles calculated by a
!  numerical function IP_calc():
!
!  (3)     IPp = IP_calc(Kz*D, AM, PM)
!  (4)     IPd = IP_calc(Kz*D, AMd, PMd)
!
!  Substitution yields a vector expression in the 'primary part' of
!  the surface drift Up, MK = (U0-Up)* IPd + Up * IPp, or
!        
!  (5)     Up = ( U0 * IPd - MK ) / (IPd - IPp)
!        
!  The 'deep part' of the surface drift is the remainder
!
!  (6)     Ud = U0 - Up
!
!  In case IPd ≃ IPp we have MK ≃ U0 * IPp and the values of Up, Ud become
!  ambigous. In any case, a new estimate of the parameter pair (AM, PM) will
!  be calculated by a best-fit method. Let proj_p(Ud) and proj_p(U_mod(z)) be
!  projections on the direction of Up. A 1-dimensional residual profile
!
!  (7)     Ur(Z) = proj_p( U_mod(Z) ) - proj_p(Ud) * exp ( -AMd (Z Kz)^PMd )
!        
!  will be fitted by an shape Up exp( -AM (Z Kz)^PM ). Here the surface 
!  value Ur(Z=0) is the value of Up (Eq. 5). Thus we can  fit the values of
!  log( log( Ur(Z)/Up ) ) (taking the log twice) with a linear expression
!  log(-AM) + PM log(Z Kz). Thereby the values of AM, PM, have been altered.
!
!  An iterative procedure is constructed as follows:
!
!  Calculate the new value of IPp using equation (3), and repeat calculation
!  corresponding to equations (5,6,7), but with the 'primary' and 'deep' parts
!  reverted. Thereby a new value of AMd is obtained by fitting a shape
!  Ud exp( -AMd (Z Kz)^PMd ).
!  Then calculate the new value of IPd using Eq. (4), and repeat so on.
!  End by determining the parameters AM, PM that provides a best to fit the
!  residual Ur(Z) of Eq. (7) with a shape Up exp( -AM (Z Kz)^PM ).
!     
!  Convergence at iteration step number ii is measured as
!  C =  | 1 - ( Up * IPp(ii-1) + Ud * IPd(ii-1) ) | / | MK |.
!
!  4. Implementation in WW3:
!
!  Code lines are added to WW3 under the two compile switches 'W3_STVP' and
!  'W3_MFIT'. First, follow the implementations as described in w3stvpmd.ftn.
!
!  a) w3odatmd.ftn
!
!     Real arrays USVP(JSEA,1:4), USVP(JSEA,4+SPND) and XSMH(JSEA,1:6) 
!     contain the 11 Stokes fitting parameters for the sea point
!     JSEA. USVP and XSMH are structures for later communication to the job
!     that performs field output and belong to the module w3adatmd.ftn.
!
!     In subroutine W3NOUT(), specify the MHfit parameter names
!           NOGE(10) = NOEXTR+1
!     #ifdef W3_MFIT
!         IDOUT(10, NOEXTR+1)  = 'Stokes profile fit  '
!
!  b) w3iogomd.ftn
!
!     b1) SUBROUTINE W3READFLGRD, and W3FLGRDFLAG
!        
!     b2) In subroutine W3FLGRDUPDT
!     #ifdef W3_STVP
!        ! Set NZO for use in W3DIMA, W3DIMX
!         NZO = 0
!         IF ( FLGRD(ISVP,JSVP) ) NZO = SPND
!     #ifdef W3_MFIT
!         ! For FLGRD(10,NOEXTR+1), if FLGRD( ISVP, JSVP) is False, only
!         ! the surface and integral variables are in USVP(1:NSEA,1:5)
!         IF ( FLGRD(10,NOEXTR+1) .AND. .NOT. FLGRD(ISVP, JSVP) ) NZO = 1
!        
!     b3) In subroutine W3OUTG:
!     Declare use of CALC_MFIT()
!     #ifdef W3_MFIT
!         USE W3XSMFMD, ONLY: CALC_MFIT
!     (...)
!       !
!       !/STVP! USVP(JSEA,1:3+2*SPND) is set by CALC_STVP(A). Thereafter,
!       !/MFIT! XSMH(JSEA,1:7) is determined by CALC_MFIT(A)
!       !/MFIT! The combined CALC_STVP(A); CALC_MFIT() are called
!       !/MFIT! outside and after the present loop
!
!     b4) Later in  W3OUTG() just following CALL CALC_STVP():
!     !/STVP ! Stokes drift with extended tail
!     !/STVP       IF ( FLOLOC(ISVP,JSVP)        &
!     !/MFIT            .OR.  FLOLOC(10,NOEXTR+1)  &
!     !/STVP                         ) THEN
!     !/STVP           CALL CALC_STVP(A)
!     !/STVP       END IF
!     !/MFIT       IF ( FLOLOC(10,NOEXTR+1) ) CALL CALC_MFIT()
!     !/STVP!
!        
!  c) In ww3_ounf.ftn:
!     Define the relevant NetCDF variable names etc.
!
!  d) Modify the makefile scripts (make_makefile.sh, w3_new) as
!     described in the WW3 manual chpt. 5.5
!
!     d1) make_makefile.sh:
!       for type in ...  stvp xsmf; do
!       (...)        
!       case $type in ...
!         (...)        
!     #sort:stvp:
!         stvp ) TY='upto1'
!                ID='Stokes drift vertical profile'
!                TS='STVP'
!                OK='STVP' ;;
!     #sort:xsmf:
!         xsmf ) TY='upto1'
!                ID='Stokes drift MHfit'
!                TS='MFIT'
!                OK='MFIT' ;;
!
!       (...)
!       case $prog in
!         (...)
!         ww3_shel) (...)
!                  IO="$IO w3stvpmd w3xsmfmd"
!         (...)
!         ww3_multi|ww3_multi_esmf)
!                   (...)
!                   IO="$IO w3stvpmd w3xsmfmd"
!       (...)        
!       case $mod in
!          (...)
!          'W3STVPMD'     ) modtest=w3stvpmd.o ;;
!          'W3XSMFMD'     ) modtest=w3xsmfmd.o ;;
!
!     d2) w3_new:
!
!       for key in $keys
!       do
!         case $key in
!           (...)
!           'stvp' ) cd $main_dir/ftn ; touch w3stvpmd.ftn
!                                       touch w3initmd.ftn
!                                       touch w3odatmd.ftn
!                                       touch w3gdatmd.ftn
!                                       touch w3iogomd.ftn
!                                       touch w3iogrmd.ftn
!                                       touch ww3_grid.ftn ;;
!           'mfit' ) cd $main_dir/ftn ; touch w3stvpmd.ftn
!                                       touch w3xsmfmd.ftn
!                                       touch w3odatmd.ftn
!                                       touch w3gdatmd.ftn
!                                       touch w3iogomd.ftn
!                                       touch w3iogrmd.ftn
!                                       touch ww3_ounf.ftn ;;
!
!  5. Variables and types :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!                C*10  Private
!
!     ----------------------------------------------------------------
!
!  6. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!
!     ----------------------------------------------------------------
!
!  7. Subroutines and functions used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!
!     ----------------------------------------------------------------
!
!  8. Remarks :
!
!
!  9. Switches :
!
!       !/SHRD  Switch for shared / distributed memory architecture.
!       !/DIST  Id.
!
!       !/ST2   Source term set 2 (Tolman and Chalikov)
!               ST2 is required for consistent similarity.
!

!  10. Source code :
!
      USE W3SERVMD, ONLY: EXTCDE
      USE W3ADATMD, ONLY: XSMH

      USE W3ODATMD, ONLY: NDSO, NDSE, NDST, IAPROC, NAPROC, NAPOUT
      USE W3GDATMD, ONLY: SPND, SPDS, SPBP, NSEAL, NK, NTH, SIG, DMIN
      ! Stokes drift profile data      
      USE W3ADATMD, ONLY: DW, USVP, ZK_S
      USE W3ODATMD, ONLY: NZO

      PRIVATE
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/

      ! Common Stokes drift depth distribution parameters to be set
      ! in MHfit_init()      
      real                        :: M_X, M_Y, K_S
      real, pointer               :: U_S(:), V_S(:)
      
      real                        :: KsD, AMs, PMs, AMds, PMds, AMIN, AMAX
      real                        :: PMIN, PMINd, PMAX, PMAXd, AMINd, AMAXd

      ! Look-up parameters to determine profile coefficients AM, PM

      real, allocatable           :: DZK_S(:), log_KZ(:)
      integer                     :: Zdi
      real                        :: delAMlu

      ! Local number of profile depths above the sea floor
      integer                     :: NPD

      ! Profile levels that cannot be fitted
      integer, allocatable        :: Zmask(:)
      
      ! DEBUG: Test output file id, Verboseness of test output
      integer                     :: NDSV, xsmf_verbose = 1

      PUBLIC :: CALC_MFIT
!/
      
      CONTAINS
!/ ------------------------------------------------------------------- /
!>
!> @brief Output a parametric representation of the Stokes drift profile vector
!>
!> @author C. Hansen  @date 28-Apr-2021
!>
!> @details 1. Purpose, and output data
!>
!> A parametric representation of the Stokes drift profile vector U_S(Z) is
!> provided as a sum of two expressions (suggested by J. Mattsson and C.Hansen),
!>       U_S(Z) = Up exp( -AM (Z Kz)^PM ) + Ud exp( -AMd (Z Kz)^PMd ) ,
!> where Up + Ud = U0 is the Stokes drift vector at the surface (Z=0) and
!> Z is pointing downwards.
!> A wave number Kz is related to the zero-upcrossing frequency by the
!> linear-wave dispersion relation, (2pi/Tz)^2 = g Kz tanh( Kz D). The
!> local depth D should also be provided in the output field for the bathymetry.
!>
!> (More @details to be converted to doxygen foemat from comments at file head)
!>
      SUBROUTINE CALC_MFIT ()      
!/
!/                  +-------------------------------------+
!/                  | FCOO Stokes profile parametric fit  |
!/                  |           Carsten Hansen            |
!/                  |                        FORTRAN 90   |
!/                  | Last update :         03-Oct-2019   |
!/                  +-------------------------------------+
!/
!/      
!/    Remarks: This subroutine can be called after CALC_STVP
!/             Call CALC_MFIT only under switch MFIT
!/
      USE W3GDATMD, ONLY: MAPSF, MAPSTA ! module scope:,NSEAL
!/
      IMPLICIT NONE

!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: ISEA, JSEA, IX, IY

! Initialization
! Generate common arrays to determine profile coefficients
! AM, PM that reproduce a given normalized profile IP0

      IF ( .NOT. ALLOCATED(log_KZ) ) THEN         
        call MHfit_init()
        IF ( xsmf_verbose .GT. 0 ) &
          WRITE (NDSV, *) '    Stokes profile fit initialized for',NSEAL,'points'
      END IF
     
      IF ( xsmf_verbose .GT. 0 ) WRITE (NDSV, 912), 'Stokes profile fit ..'
!
! -------------------------------------------------------------------- /
! 1.  Loop over sea points
!     
      DO JSEA=1, NSEAL
!
! -------------------------------------------------------------------- /
! 2.  Process only for water points
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
        IF ( MAPSTA(IY,IX) .GT. 0 ) &
          call Stokes_MHfit(JSEA)
          ! If under MPI only USVP(JSEA,1:5) is communicated to XUSVP for output
        if ( NZO < SPND ) &
          USVP(JSEA,4+NZO:3+2*NZO) = USVP(JSEA,4+SPND:3+SPND+NZO)
      END DO
      
      IF ( xsmf_verbose .GT. 0 ) WRITE (NDSV, 912), 'Stokes fits done'
      
  912 FORMAT ('  ',A)

      END SUBROUTINE CALC_MFIT    
      
      
      SUBROUTINE Stokes_MHfit(JSEA)

      IMPLICIT NONE
      INTEGER, INTENT(IN)         :: JSEA
      
! Parametric profile fitting
      integer                     :: ISEA, IZ
      real                        :: UsD, VsD, AM, PM, AMd, PMd
      real                        :: SIC
      real                        :: rel_ms, UN, fac

!/ ------------------------------------------------------------------- /

!  Calculate fit to Mattsson-Hansen parametric profile
! --------------------
!
!     MHfit profile fit. Find the surface vectors Up, Ud and
!     the parameters AM, PM, AMd, PMd so that with a minimum e(:),
!     (U_S(i),V_S(i)) = e(i) + Up*exp(-AM*ZK(i)**PM) + Ud*exp(-AMd*ZK(i)**PMd)
!     where
!     ZK(i) = Z(i) * Kz, Kz  = ( TPI / Tz)**2 / GRAV
!     J. Mattsson first suggestion for the parameters:
!     AM=2.5, PM=0.45, Ud=0

      if ( xsmf_verbose > 1 .and. JSEA == 1 ) &
           WRITE (NDSV, *), JSEA, '  Fit the MHfit profile'
      

      K_S = USVP(JSEA,1)
      M_X = USVP(JSEA,2)
      M_Y = USVP(JSEA,3)
      U_S => USVP(JSEA,4:3+SPND)
      V_S => USVP(JSEA,4+SPND:3+SPND*2)
      
      ! Feed with number of fitting cycles in MH_fit
      ! A minimum of 1 cycle + a half cycle is required:
      ! Cycle 1: fit AM,PM -> update Ud,Up -> fit AMd,PMd -> update Up,Ud 
      ! Half cycle: fit AM,PM -> update Ud,Up

      SIC = 1.
      
      call MH_fit(JSEA, AM, PM, AMd, PMd, UsD, VsD, SIC)

      USVP(JSEA,2)  = M_X
      USVP(JSEA,3)  = M_Y

      XSMH(JSEA,1)  = UsD     ! Surface drift of 'deep' part.
      XSMH(JSEA,2)  = VsD     ! Primary part drift is UsP=U_S-UsD.

      XSMH(JSEA,3)  = AM      ! MHfit parameter, primary part
      XSMH(JSEA,4)  = PM      ! MHfit parameter, primary part

      XSMH(JSEA,5)  = AMd     ! Deep part MHfit parameter.
      XSMH(JSEA,6)  = PMd     ! Deep part MHfit parameter. May be constant PMds

      XSMH(JSEA,7)  = SIC     ! Quality integer
      ! SIC = si*100 + min(int(qiUd+3*niter),9)*10 + min(int(limAP+3*uopp),9)
      ! si : Scatter index % std(|U_S - U_mod| * dz)/std(|U_mod|)
      !
      ! qiUd  =  0:none, 1:Udcap, 2:Ud=0
      ! limAP =  0:none, 1:lim A or P, 2: lim Ad or Pd only
      ! uopp  =  0:false, 1:true, 2: IP2_min

      ! Reconstruction:
      ! (...)

      if ( xsmf_verbose > 1  .and. JSEA == NSEAL )  &
         WRITE (NDSV, *), JSEA, '  MHfit done. SIC =', SIC

      END SUBROUTINE Stokes_MHfit

      
      SUBROUTINE MHfit_init()

      IMPLICIT NONE
      integer                     :: iret, IZ
      real                        :: XKZ
      logical                     :: OPENED
      
      NDSV = NDST
      ! Highly verbose output only with test output
      if ( xsmf_verbose .gt. 0 ) then
        INQUIRE (NDSV,OPENED=OPENED)
        IF ( .NOT. OPENED ) THEN
          NDSV = NDSO
          IF ( .NOT. IAPROC == NAPOUT )  xsmf_verbose = 0
        END IF
      end if

      ! Initializing values (based on tests 2020 with ST4-Romero)
      AMds = 0.6; PMds = 1.5
      AMs = 2.5; PMs = 0.57

      ! Overall min/max values in MH_fit
      PMIN = 0.1  ! PM
      PMAX = 2.0  ! PM
      AMIN = 0.2  ! AM
      AMAX = 15.0  ! AM        
      PMINd = 0.5  ! PMd
      PMAXd = 2.4  ! PMd
      AMINd = 0.4  ! AMd        
      AMAXd = 5.0  ! AMd        

      allocate( DZK_S(SPND), stat=iret )

      ! DZK_S(:) = K_S * approximate thicknesses (m) of each layer
      ! In w3stvpmd, the dimensionless profile depths have been defined as
      !      ZK_S(IZ) = SPDS*(XZK**((IZ-1.)**SPBP) - 1.), where 
      ! XZK = 2.**((SPND*1.-1.)**(-SPBP))
      ! Let us untroduce a set of half-way dimlensionless depths
      ! ZK_SH(IZ) = SPDS*(XZKH**((IZ-1.)**SPBP) - 1.), IZ=1,SPND*2-1, where
      ! XZKH = 2.**((SPND*2.-2.)**(-SPBP).
      ! ZK_S(1:SPND) is every second of these depths.
      ! The uppermost layer thickness DZK_S(1) is the first of these half-way
      ! intervals
      ! DZK_S(1) = ZK_SH(2) = SPDS*(XZKH**(((2-1.)**SPBP) - 1.)
      !          = SPDS*(2.**((SPND*2.-2.)**(-SPBP) - 1.)
      !          = SPDS*(2.**(0.5/(SPND-1.)**SPBP) - 1)
      ! and each of the rest of the thicknesses DZK_S(2:) span over two of
      ! these level intervals
      ! DZK_S(IZ) = ZK_SH(IZ*2) - ZK_SH((IZ-1)*2), IZ = 2, SPND
      ! For IZ*2 = 2, 4, 6, 8, ...
      ! ZK_SH(IZ*2) = SPDS*(XZKH**((IZ*2-1.)**SPBP) - 1.)
      !             = SPDS*(2.**[(SPND*2.-2.)**(-SPBP) * (IZ*2-1)**SPBP] - 1 )
      !             = SPDS*(2.**[( 0.5*(IH-1.)/(SPND-1.) )**SPBP ] - 1 )      
      ! For  (IZ-1)*2 = 2, 4, 6, ... :
      ! ZK_SH((IZ-1)*2) = SPDS*(XZKH**(((IZ-1)*2-1.)**SPBP) - 1.)
      !            = SPDS*(2.**[( 0.5*(IH0-1.)/(SPND-1.) )**SPBP ] - 1 )
      ! Let the array DZK_S(IZ) hold ZK_SH(IZ*2), then subtract ZK_SH((IZ-1)*2)
      DO IZ = 1, SPND
         DZK_S(IZ) = SPDS*(2.**( ( 0.5*(IZ*2-1.)/(SPND-1.) )**SPBP ) - 1 )
      END DO
      DO IZ = 2, SPND
         DZK_S(IZ) = DZK_S(IZ) -  DZK_S(IZ-1)
      END DO
      
      ! For fitting the profiles, we use
      allocate( log_KZ(SPND), stat=iret )
      log_KZ(2:SPND) = log(ZK_S(2:SPND))
      ! Dummy value at surface
      log_KZ(1) = 0.

      ! Max depth index for initial deep fit is for KZ = 1.0
      ! Zdi = min(int((log(2.)/log(XKZ))**(1./SPBP)) + 1, SPND - 3)
      Zdi = SPND - 3

      allocate( Zmask(SPND), stat=iret )
      Zmask = 0

      end SUBROUTINE MHfit_init


      SUBROUTINE MH_fit(JSEA, AM, PM, AMd, PMd, UdX, UdY, SIC)
        
      implicit none

!     Mattsson-Hansen profile fit. Find the surface vectors Up, Ud and
!     the parameters AM, PM, AMd, PMd so that with a minimum rms error |e(:)|,
!     (U_S(i),V_S(i)) = e(i) + Up*exp(-AM*ZK(i)**PM) + Ud*exp(-AMd*ZK(i)**PMd)
!     where
!     ZK(i) = Z(i) * Kz, Kz  = ( TPI / Tz)**2 / GRAV


      integer, intent(in)      :: JSEA
      real, intent(out)        :: AM, PM, AMd, PMd, UdX, UdY
      real, intent(inout)      :: SIC
      
      complex                  :: MK, Ud, U0, Up, unit, Us, Ut
        
      real                     :: IPp, IPd, tmp, Udl, Upr, Mpr, IP2_min, mlU
      real, allocatable, save  :: Pd(:), Pp(:), Umr(:)
      integer                  :: iter, numPos, ii, im, nneg
      integer                  :: ISEA, IZ, Zt
      complex                  :: Upi(5),Udi(5)
      real                     :: IPpi(5),AMi(5),PMi(5),AMdi(5),PMdi(5),IPdi(5)
      real                     :: msi(5)
      logical                  :: wi(5), warning_issued, capped
      real                     :: fi, fiU, Umr_min, Umd_min, dIP, dIP_min
      real                     :: si, cit, sm, qcr, qcd, ms_U
      real                     :: IPm, IP
      integer                  :: NA, niter, qiUd, uopp, limAP
      character(len=5)         :: fit_kind

      if ( .not. allocated(Umr) )  allocate( Umr(SPND), stat=ii )
      if ( .not. allocated(Pd) )   allocate( Pd(SPND), stat=ii )
      if ( .not. allocated(Pp) )   allocate( Pp(SPND), stat=ii )
      
      ! Initial quality index parameter
      qiUd=0
      uopp=0
      limAP=0
      
      niter = min(max(int( SIC+0.1 ), 1 ),3)

      dIP_min = 0.01

!     Initially, assume standard parameters of the profiles
      AMd = AMds; PMd = PMds; AM = AMs; PM = PMs
      
      ! Convert input vectors to complex numbers, and normalize M
      MK = cmplx( M_X, M_Y ) * K_S
      U0 = cmplx( U_S(1), V_S(1) )

#ifdef W3_DIST
      ISEA       = IAPROC + (JSEA-1)*NAPROC
#endif         
#ifdef W3_SHRD
      ISEA       = JSEA
#endif         
      KsD  = K_S * MAX ( DMIN, DW(ISEA) )

      ! NPD: Number of depths of the discrete profile
      NPD = SPND ! Configured value. To be reduced if below the sea bed     
      do IZ = SPND,1,-1
        if ( ZK_S(IZ) > KsD ) cycle
        NPD = IZ
        exit
      end do
      

      ! If there is a turning depth, estimate a deep profile fit
      
      ! Longest vector deviation from surface value and its depth index
      Zt = NPD
      call S_span(Us,Zt)
      if ( Zt > Zdi .or. NPD - Zt <= 3 ) Zt = NPD

      ! Given PMd = PMds, estimate Ud and AMd
      if ( Zt /= NPD ) then

        ! Deep profile top:
        Ut = cmplx(U_S(Zt), V_S(Zt))
        ! Projection on the direction of Ut
        unit = Ut / abs(Ut)
        Umr(Zt:NPD) = U_S(Zt:NPD) * real(unit) + V_S(Zt:NPD) * aimag(unit)
        Umr(Zt:NPD) = Umr(Zt:NPD)/abs(Ut)
        Umr_min = 0.01
        do ii = Zt,NPD
          if  (Umr(ii) < Umr_min ) exit
        end do
        if (ii == NPD) then
          numPos = ii
        else
          numPos = ii - 1
        end if
        Zmask=0
        call lin_reg(-log(Umr(Zt:numPos)), ZK_S(Zt:numPos)**PMd, &
                     Zmask(Zt:numPos), AMd, mlU, AMINd, AMAXd)
        Ud = Ut * exp(-mlU)
        
!       The 'upper part' of the surface drift is the remainder
        Up = U0 - Ud
        uopp = 1
      end if

      ! Integrated drift of the deep profile
      IPd = IP_calc(KsD, AMd, PMd)

      if ( Zt == NPD ) then
!       Factorize the normalized M*K in an 'primary part' and a 'deep part' as
!       MK = Ud * IPd + Up * IPp, where
!       U0 = Ud + Up, and where IPd, IPp are normalized profile shapes
        IPp = IP_calc(KsD, AM, PM)

!       Substitution yields a vector expression in one unknown part Up,
!       MK = (U0-Up)* IPd + Up * IPp,  or
        Up = ( U0 * IPd - MK ) / (IPd - IPp)
!       The 'deep part' of the surface drift is the remainder
        Ud = U0 - Up
      end if
      
      ! Handling of small drift
      if ( abs(U0) < 0.003 .and. abs(MK) < 0.003 ) then
        ! Neglect the details when very weak Stokes drift
        SIC = -3276.9
        return
      end if
      
      ! Smallest drift values allowed in the improved estimates
      ! of (AM, PM, AMd) in the iteration steps below

      Umr_min = abs(Up) * exp( - AMs * ZK_S(NPD)**(PMs*2) )
      Umd_min = abs(Ud) * exp( - AMd * ZK_S(NPD)**(PMds*2) )

      ! The two parametric profile components
      Pp = exp ( - AM * ZK_S ** PM )
      Pd = exp ( - AMd * ZK_S ** PMd )

      warning_issued = .false.
      
! Cycles:

      ! The number (niter) of fitting cycles must be at least 1 + a half cycle:
      ! Cycle 1: fit AM,PM -> update Ud,Up -> fit AMd,PMd -> update Up,Ud.
      ! A half cycle (niter+1) is added to update Ud,Up again.
      ! If Ud is nearly zero, then an extra cycle (up to niter+2) is taken.
      
      do iter=1,niter+2

!  Register the primary profile parameters
      
        AMi(iter)  = AM
        PMi(iter)  = PM
        AMdi(iter) = AMd
        PMdi(iter) = PMd
        Upi(iter)  = Up
        Udi(iter)  = Ud
        IPpi(iter) = IPp
        IPdi(iter) = IPd
        wi(iter)   = warning_issued

        ! Mean square deviation from predicted profile
        msi(iter) = calc_msd( Up, Ud, Pp, Pd, NPD)                 

        warning_issued = .false.

!
!  step 1: Apply the model profile Um(z ) to yield improved estimates
!       of (AM, PM) determined as the best fit to the residual
!       Umr(z) =  proj(Um(z)) - proj(Ud)*Pd(z),
!       where proj is the projection on the direction of Up given from
!       the previous iteration, and Pd(z) = exp ( - AMd * zK ** PMd ) is
!       given from the parametric deep parameters of the previous iteration.
!       Note that Pd(z=0) = 1 and Umr(z=0) = Up

        ! Projection Ud * unit of vector Ud on direction of Up
        unit = Up/abs(Up)
        Umr = U_S * real(unit) + V_S * aimag(unit)
        
        if (abs(Ud)>0.) then
          tmp = real(Ud) * real(unit) + aimag(Ud) * aimag(unit)
          Umr = Umr - tmp * Pd
        end if
          
        fit_kind='upper'
        call fit_proj(AM,PM,log_KZ,Umr,Umr_min,'upper',im)
        if (im/=0) limAP = 1
        if (PM > 2.*PMs .or. 3.*PM < PMs+PMIN ) warning_issued = .true.
        
!  step 2: Update the profile parts, given the old IPd and new AM, PM

        if (xsmf_verbose > 2 ) &
                write (NDSV,*) ' ISEA=', ISEA, '  IP_calc( ', AM, PM,')'
        IPp = IP_calc( KsD, AM, PM)

        Pp = exp ( - AM * ZK_S ** PM )
        
!       Make a match to MK of the surface vectors for the deep and primary parts
        dIP = IPd/IPp - 1. 
        if ( abs(dIP) > dIP_min ) then
          ! Finalize the loop if iter >= niter + 1 half step
          if ( iter > niter .and. .not. capped ) exit
        
          ! MK = (U0-Ud)* IPp + Ud * IPd, which means that
          Ud = ( MK/IPp - U0 ) / dIP
          Up = U0 - Ud
          ! Limiter to the magnitude of Ud. Require |Up| >= 2 * |Ud|
          ! This is with the purpose of preserving positive values to
          ! some depth of the 1-d profile Umr below
          call Ud_cap(Ud,Up,U0,Us,capped=capped)
        
        else
          ! If dIP is of small magnitude, We may better neglect the deep part
          Ud = 0.
          Up = U0
          if ( iter > niter+1) exit
          capped = .False.
          qiUd = 2
          cycle ! Repeat steps 1, 2 (even for iter=niter+2)
        end if
          
!  step 3: Determine AMd as the best fit of the secondary, 'deep',
!       profile given its surface drift vector Ud

        if (xsmf_verbose > 2 ) write (NDSV,*), '  Fit AMd'
            
        ! Find  Umr(z)= proj_d(Um) - proj_d(Up)*Pp
        ! projected on the direction of Ud
        
        Udl = abs(Ud)
        unit = Ud /  Udl
        tmp = real(Up) * real(unit) + aimag(Up) * aimag(unit)
        
        Umr = U_S * real(unit) + V_S * aimag(unit) - tmp * Pp
        ! 
        AMd = AMds ! For a ref. profile within fit_proj
        PMd = PMds ! - do. -
        fit_kind='deep '
        
        ! Fit for two parameters AMd and PMd
        call fit_proj(AMd,PMd,log_KZ,Umr,Umd_min,fit_kind,im)
        
        ! If AMd or PMd is an outlier, fix PMd to a standard values
        if (.not. (PMd>PM .and. PMd<PMAXd .and. AMd>AMINd .and. AMd<AM)) then
          if ( im/=0 .and. limAP == 0 ) limAP = 2
          AMd = AMds
          if (PMd < PM) then
             PMd = PM
             if (PMd > PMds) PMd = PMds
          else if (PMd > PMAXd) then
             PMd = PMAXd
          end if
          call fit_proj(AMd,PMd,log_KZ,Umr,Umd_min,fit_kind,im,PMds)
        end if
          
        if ( iter > niter ) then ! if capped
           Pd = exp ( - AMd * ZK_S ** PMd )
           exit
        end if
        
        ! Apply a low limiter to AMd
        AMd = max(AMd, AMINd)
        if ( AMd ==AMINd ) limAP=2
        
          
!  step 4: Update the profile partitions based on the old IPp and a new AMd, PMd
        if (xsmf_verbose > 2 ) write (NDSV,*), '  IP_calc(KsD', AMd, PMd,')'
        
        IPd = IP_calc( KsD, AMd, PMd)
        
!       Make a match to MK of the surface vectors for the profile partitions
        dIP = 1. - IPp/IPd
        
        if ( dIP > dIP_min ) then
          ! Matching MK = (U0-Ud)* IPp + Ud * IPd:
          Up = ( U0 - MK/IPd ) / dIP
          Ud = U0 - Up
          Pd = exp ( - AMd * ZK_S ** PMd )
          ! Limiting cap to the magnitude of Ud. Require |Up| >= 2 * |Ud|
          ! This is with the purpose of preserving positive values to
          ! some depth of the 1-d profile Umr in the refinement

          if (iter <= niter ) call Ud_cap(Ud,Up,U0,Us,capped=capped)
          if (capped) qiUd=1
        else
          Up = U0
          Ud = 0.
          Pd = 0.
          qiUd = 2
          if ( IPd < -dIP_min ) then
            capped=.True.
            qiUd=1
          end if
        end if
          
!     Repeat steps 1-2 with the new AMd to get a more precise fit of AM, PM
        
      end do ! iter = 1,niter+2

      if ( AMd <= AMINd ) warning_issued = .true.
      
!     Scatter index of the final fit
      sm  = calc_msd( Up, Ud, Pp, Pd, NPD)

      if ( xsmf_verbose > 2 &
           .or. ( xsmf_verbose  > 1 .and. warning_issued ) ) then
        write (NDSV,*), '  U0,     Ud,     AM,  PM,  IPp,  AMd,  IPd, si = '
        write (NDSV,*), U0, Ud, AM, PM, IPp, AMd, IPd, si

        write (NDSV,*), '    MK =  ', MK
        write (NDSV,*), '    Upi = ', Upi(1:iter) , Up
        write (NDSV,*), '   IPpi = ', IPpi(1:iter), IPp
        write (NDSV,*), '    AMi = ', AMi(1:iter) , AM
        write (NDSV,*), '    PMi = ', PMi(1:iter) , PM
        write (NDSV,*), '   AMdi = ', AMdi(1:iter), AM
        write (NDSV,*), '    si = ', sqrt( msi(1:iter) )
      end if
        
      im = iter + 1
      
      if ( niter > 0 ) then
       
        do ii=1,niter
          if ( msi(ii) > sm - 0.01 ) cycle ! Ignore insignificant difference
          if ( xsmf_verbose > 0 ) &
            write (NDSV,*) '  WARNING in w3xsmfmd, ISEA=', ISEA, &
                ': msi(',im,')/msi(',ii,')=',sm,'/',msi(ii)
          if (ii==1 ) cycle
          im = ii
          sm = msi(ii)
        end do
        
        if ( im <= niter ) then
          AM  = AMi(im)
          PM  = PMi(im)
          AMd = AMdi(im)
          PMd = PMdi(im)
          Up  = Upi(im)
          Ud  = Udi(im)
          IPp = IPpi(im)
          IPd = IPdi(im)
          warning_issued = wi(im)
        
        end if          
      end if
      
      niter = max(im - 2, 0)
      
      if ( AMd >= AMAX .and. xsmf_verbose > 0 ) warning_issued = .true.
          
      UdX = real(Ud)
      UdY = aimag(Ud)
      
      ! sqrt(|U_S^2|)
      ms_U = sum((U_S(1:NPD)*U_S(1:NPD) + V_S(1:NPD)*V_S(1:NPD))* DZK_S(1:NPD))/ZK_S(NPD)
      
      ! Check of the match to MK of the surface vectors:
      ! cit = abs ( MK - ( Up * IPp + Ud * IPd ) ) / ms_U
      
      ! Scatter index of the fit; si = sqrt(mean(U_S - U_fit)^2)/mean(U_S^2))
      si = sqrt( sm / ms_U)

      niter = max(niter-1,0)
      SIC = min(int(qiUd+3*niter),9)*10 + min(int(limAP+3*uopp),9)
      
      si = 1000.*si
      if (si > 326.) si = 0.5*(326.-si)
      si =  max(si,-325.)
      if ( si >= 0 ) then
        SIC = nint(si)*100 + SIC
      else
        SIC = nint(si)*100 - SIC
      end if

      if ( xsmf_verbose == 0 .and. .not. warning_issued ) return

      ! Ignore deviations warnings in calm sea
      if ( abs(U0) + abs(MK) < 0.02 ) return
      
      if ( .not. ( AM > AMIN .and. 3.*PM >= PMs+PMIN .and. AMd >= AMINd ) ) then
        write (NDSV,*) '  WARNING in w3xsmfmd: Very small AM, PM, or AMd'
        warning_issued = .true.
      end if

      if ( .not. ( AM < AMAX .and. AMd < AMAX .and. PM <= 2.*PMs ) ) then
        write (NDSV,*) '  WARNING in w3xsmfmd: Very large AM, PM, or AMd'
        warning_issued = .true.
      end if

      if ( .not. abs(Ud) < abs(U0) + abs(MK) ) then
        write (NDSV,*), '  WARNING in w3xsmfmd:', &
                ' Strong secondary Stokes drift |Ud| > |U0| + |MK|'
        warning_issued = .true.
      end if

      if ( xsmf_verbose > 1 .or. warning_issued ) then
#ifdef W3_DIST
        ISEA       = IAPROC + (JSEA-1)*NAPROC
#endif         
#ifdef W3_SHRD
        ISEA       = JSEA
#endif         
        write (NDSV,*) '  ISEA=', ISEA, ', niter=',niter
        write (NDSV,*) '    KsD, AM, PM, AMd = ', KsD, AM, PM, AMd
        write (NDSV,*) '    |U0|,   |Ud|,   |MK|,   IPd,   IPp ='
        write (NDSV,*) abs(U0), abs(Ud), abs(MK), IPd, IPp
      end if

      end SUBROUTINE MH_fit


      SUBROUTINE fit_expr(numD,log_KZ,Um,AM,PM,Um_min,fit_kind,P0,A0)

      !  Determine AM, PM as a best fit of an expression
      !  Ufit = Um(1) exp(-AM * KZ**PM) to a profile Um(KZ).
      implicit none

      integer, intent(inout)      :: numD
      real, intent(in)            :: log_KZ(:), Um(:)
      real, intent(in)            :: Um_min
      real, intent(inout)         :: AM, PM
      character(len=5), intent(in):: fit_kind
      real, intent(in), optional  :: P0
      real, intent(in), optional  :: A0

      real, allocatable, save     :: loglog_ref(:), log_Um(:)
      real                        :: log_Up, log_AM, log_Um_min, logAM0, amx, amn
      integer                     :: ii, num_unmasked

      if ( .not. allocated(loglog_ref) ) allocate( loglog_ref(SPND), stat=ii )
      if ( .not. allocated(log_Um) )     allocate( log_Um(SPND), stat=ii )

      ! Zmask = 1 for the depths that are ignored in the fit
      Zmask = 0
      
      ! Profile model:
      ! Ufit = Um(1) * exp (-AM * KZ**PM)

      ! Regression shape valid for KZ > 0:
      ! log_Um = -log( Um/Um(1) ) ~ AM * KZ**PM

      ! Reference shape:
      ! log_ref(2:numD) = AM * KZ(2:numD)**PM
      ! The log of log_Um is linear in log_KZ:
      if ( present(A0) ) AM = A0 ! Best fit of PM, only, given AM=A0
      if ( present(P0) ) PM = P0 ! Best fit of log_AM, only, given PM=P0
      
      logAM0 = log(AM)
      loglog_ref(2:SPND) = logAM0 + PM * log_KZ(2:SPND)
      
      log_Um_min = -log(Um_min/Um(1))

      where ( Um > Um_min )
        log_Um = -log( Um/Um(1) )
      elsewhere
        ! Avoid depths of near negative speed:      
        log_Um = log_Um_min
        ! Register the depths that are ignored
        Zmask = 1
      end where
      
      ! Set dummy values at the surface and below the sea floor (numD+1:SPND)
      ! (not to be used in lin_reg)
      loglog_ref(1) = loglog_ref(2) - 1.
      log_Um(1) = exp(loglog_ref(1))
      log_Um(numD+1:SPND) = exp(loglog_ref(numD))
      
      ! First, avoid the smallest (possibly negative) values of log_Um: 
      where ( log_Um < log_Um(1) )
        log_Um = log_Um(1)
        Zmask = 1
      end where

      ! Take the logarithm again to attain a linear form
      log_Um(2:numD) = log( log_Um(2:numD) )
      ! Dummy values at the surface and below the sea floor
      log_Um(1) = loglog_ref(1)
      log_Um(numD+1:SPND) = loglog_ref(numD+1:SPND)
      
      ! Avoid depths where Um is not decreasing like 'in the vicinity of'
      ! the reference shape Um(1) * exp(-AM * KZ**PM)

      where ( log_Um < loglog_ref - 1. )
        log_Um = loglog_ref - 1.
        Zmask = 1
      end where      
      ! where ( log_Um < 1/e*log_ref )
        
      where ( log_Um > loglog_ref + 1. )
        log_Um = loglog_ref + 1.
        Zmask = 1
      end where
      ! where ( log_Um > e*log_ref )

      if ( index(fit_kind,'deep') > 0 ) then
        amx = PMAXd
        amn = PMINd
      else
        amx = PMAX
        amn = PMIN
      end if
      
      ! Number of depths that do not go into the fitting.
      num_unmasked = numD - 1 - sum(Zmask(2:numD))
      
      ! Linear regression fit
      if (xsmf_verbose > 2 ) write (NDSV,*), '  lin_reg () ... '
      if ( num_unmasked == 0 .OR. ( num_unmasked == 1 .AND. &
            .NOT. ( present(A0) .OR. present(P0) ) ) ) then
        ! Number of depths that go into the fitting. Returned in numD
        numD = num_unmasked
        return
      end if
      
      if ( present(A0) ) then
        ! Best fit of PM, only, given AM=A0
        log_AM=log(A0)
        call lin_reg1( log_Um(2:numD), log_KZ(2:numD), Zmask(2:numD), PM, log_AM)
      else if ( present(P0) ) then
        ! Best fit of log_AM, only, given PM=P0
        call lin_reg0( log_Um(2:numD), log_KZ(2:numD), Zmask(2:numD), P0, log_AM)
        AM = exp(log_AM)
      else
        ! Linear regression best fit PM, log_AM
        call lin_reg( log_Um(2:numD), log_KZ(2:numD), Zmask(2:numD), PM, log_AM, amn, amx)
        AM = exp(log_AM)
      end if

      ! Number of depths that go into the fitting. Returned in numD
      numD = num_unmasked      
      
      end SUBROUTINE fit_expr

      
      SUBROUTINE fit_proj(AM,PM,log_KZ,Umr,Umr_min,fit_kind,limi,P0)
        
        implicit none
        real, intent(inout)         :: AM, PM
        real, intent(in)            :: log_KZ(:), Umr(:), Umr_min
        character(len=5), intent(in):: fit_kind
        integer, intent(out)        :: limi
        real, intent(in), optional  :: P0
        real                        :: tmp
        integer                     :: ii, nneg, numPos

        limi=0
        
        if (present(P0)) then
           numPos=NPD
           call fit_expr(numPos,log_KZ,Umr,AM,PM,Umr_min,fit_kind,P0)
           ! qcr = (NPD - numPos) / NPD
           return
        end if
        
!       Determine the upper range of depths (1..numPos) where Umr has positive
!       values. There must be no more than 0.2 * NPD negative or near zero values
!       counted from the surface
        nneg = int(0.2 * NPD)
        do ii = 1,NPD
          if  (Umr(ii) > Umr_min ) then
             numPos = ii
             cycle
          end if
          nneg = nneg - 1
          if ( nneg == 0 ) exit
        end do
                
        if (numPos < 4 ) then
          if (numPos > 1 ) then
            PM=PMs
            call fit_expr(numPos,log_KZ,Umr,AM,PM,Umr_min,fit_kind,PM)
            ! Fraction of depths that do not go into the fitting:
            ! qcr = (NPD - numPos) / NPD
          end if
        else
!       Determine AM, PM as the best fit of the profile, given the surface
!       drift Up. Although Umr(1) is always positive, any of Umr(2:numPos)
!       might be negative. fit_expr will first exchange negative or near-zero
!       values with the previous (input) profile.
          call fit_expr(numPos,log_KZ,Umr,AM,PM,Umr_min,fit_kind)
          ! qcr = (NPD - numPos) / NPD
          
          ! Limiters to PM
          if ( index(fit_kind,'deep') == 0 &
               .and. ( 2.*PM > PMAX+PMs .or. 2.*PM < 3.*PMIN) ) then
            limi=1
            if (2.*PM > PMAX+PMs) then
              tmp = PMAX-PMs
              PM = PMAX - tmp/(2.*(PM-PMs)/tmp + 1.)
            else
              tmp = PMIN
              PM = PMIN * (1. + 1./(5.-2.*PM/PMIN))              
            end if
            ! Best fit value of AM given the new PM
            call fit_expr(numPos,log_KZ,Umr,AM,PM,Umr_min,fit_kind,PM)
            ! qcr = (NPD - numPos) / NPD
            return
          end if

        end if ! numPos < 4

        if ( index(fit_kind,'deep') > 0 ) return
        
        ! Limiters to AM
        if ( AM < AMIN .or. 2. * AM > AMAX+AMs ) then
          limi=1
          if ( AM < AMIN ) then
            AM = AMIN
          else
            tmp = AMAX-AMs
            AM = AMAX - tmp/(2.*(AM-AMs)/tmp + 1.)
          end if
          ! Best fit value of PM given the new AM
          call fit_expr(numPos,log_KZ,Umr,AM,PM,Umr_min,fit_kind,PM,AM)
          ! qcr = (NPD - numPos) / NPD
        end if
      
      end SUBROUTINE fit_proj

      
      SUBROUTINE Ud_cap(Ud,Up,U0,Us,capped)
        ! Limiter to the magnitude of Ud. We presume that U0 = Up + Ud
        implicit none
        complex, intent(inout)  :: Ud,Up
        complex, intent(in)     :: U0,Us
        logical,optional,intent(out) :: capped
        real                    :: tmp
        complex                 :: Uds, unit
        logical                 :: iscapped

        iscapped = .False.
        Uds = Ud+Us-U0
        !if ( 2.*abs(Ud) > abs(Up) ) then
        if ( abs(Uds) > abs(Us) ) then
          ! We presume U0 = Up + Ud
          unit = Us/abs(Us)
          ! Projection of Ud+Us-U0 in the direction of Us
          tmp = real(Uds) * real(unit) + aimag(Uds) * aimag(unit)
          
          if ( tmp < 0. ) then
            ! Shift Ud in the direction of Us so that Ud and Us become
            ! perpendicular. If |Uds|>|Us| reduce to |Us|
            Uds = Uds - (tmp/abs(Us))*Us
            Ud = Uds*min(1.,abs(Us)/abs(Uds)) - (Us-U0)
            Up = U0 - Ud
            iscapped = .True.
          
          else ! if ( abs(Ud) > abs(Up) ) then
            ! Shift Ud in the direction of -U0 so that it's projection
            ! on U0 becomes 0.5*U0. If |Uds|>|Us| reduce to |Us| 
            Uds = Uds - (tmp/abs(Us) - 0.5)*Us
            Ud = Uds*min(1.,abs(Us)/abs(Uds)) - (Us-U0)
            Up = U0 - Ud
            iscapped = .True.
          end if
        end if
        if ( present(capped) ) capped = iscapped
      end SUBROUTINE Ud_cap

      SUBROUTINE S_span(Ut,Zt)
        ! Check for a deep opposing drift
        implicit none
        complex, intent(out)    :: Ut
        integer, intent(inout)  :: Zt
        real                    :: us0,vs0,d20,d21
        integer                 :: ii
        us0 = U_S(1)
        vs0 = V_S(1)
        d20 = 0.
        do ii = Zt, 2, -1
           d21 = (U_S(ii)-us0)**2 + (V_S(ii)-vs0)**2
           if ( d20 > d21 ) exit
           d20 = d21
         end do
        if ( d20 /= d21 ) Zt = ii + 1
        ! Zt = maxloc( (us0-U_S(2:Zt))**2 + (vs0-V_S(2:Zt))**2 )(1) + 1
        Ut = cmplx( us0-U_S(Zt), vs0-V_S(Zt) )

      end SUBROUTINE S_span

      ! SUBROUTINE S_normal(Zt)
      !   ! Check for the significant angle deviation
      !   implicit none
      !   integer, intent(inout)  :: Zt
      !   real                    :: d20,d21,sin_anglt
      !   complex                 :: P, unit0
      !   integer                 :: ii
! 
      !   ! Projections on the normal direction of the surface drift
      !   unit0=cmplx(-V_S(1),U_S(1))
      !   unit0=unit0/abs(unit0)
      !   d20 = 0.
      !   ! First maximum normal projection
      !   do ii = 2, Zt
      !     d21 = abs( cmplx(U_S(ii),V_S(ii)) * unit0 )
      !     if ( d21 < d20 ) exit
      !     d20 = d21
      !   end do
! 
      !   ! Mean angle to the surface drift
      !   P = cmplx(sum(U_S(ii:Zt)),sum(V_S(ii:Zt)))
      !   sin_anglt = abs((P*unit0)/abs(P))
! 
      !   ! First depth where the mean angle is exceeded
      !   do ii = ii, Zt           
      !     P = cmplx(U_S(ii),V_S(ii))
      !     if ( abs(P*unit0) > sin_anglt*abs(P) ) exit
      !   end do
      !   
      !   Zt = ii
      !   
      ! end SUBROUTINE S_normal          
        
      real FUNCTION calc_msd( Up, Ud, Pp, Pd, NPD)

      ! Mean of squares of modulus deviation between a profile U_S, V_S
      ! and a parametric fit Up * Pp + Ud * Pd
      implicit none
      complex, intent(in)               :: Up, Ud
      real, dimension(:), intent(in)    :: Pp, Pd
      integer, intent(in)               :: NPD
      real, allocatable, save           :: Um(:)
      real                              :: tmp
      integer                           :: ii
      
      if ( .not. allocated(Um) )  allocate( Um(SPND), stat=ii )      
 
      Um(2:NPD) = ( U_S(2:NPD) - real(Up) * Pp(2:NPD)     &
                       - real(Ud) * Pd(2:NPD) )
      tmp = sum( Um(2:NPD) * Um(2:NPD) * DZK_S(2:NPD) )
      Um(2:NPD) = ( V_S(2:NPD) - aimag(Up) * Pp(2:NPD)    &
                       - aimag(Ud) * Pd(2:NPD) )     
      tmp = tmp + sum( Um(2:NPD) * Um(2:NPD) * DZK_S(2:NPD) )
      calc_msd = tmp / ZK_S(NPD)
      
      end FUNCTION calc_msd
        
      SUBROUTINE lin_reg0( yy, xx, mask, a0, bb)
      ! Find the y-intercept (bb) of
      ! yf = a0 * xx + bb, where ||yy - yf|| is a minimum
      IMPLICIT NONE
      real, intent(in)            :: yy(:), xx(:)
      real, intent(in)            :: a0
      real, intent(out)           :: bb
      integer, intent(in)         :: mask(:)
     
      integer :: lenArr

      lenArr =  size(yy) - sum(mask(:))

      bb = sum(yy* (1-mask(:))) / lenArr - a0 * ( sum(xx* (1-mask(:))) / lenArr )

      end SUBROUTINE lin_reg0


      SUBROUTINE lin_reg1( yy, xx, mask, aa, b0)
      ! Find the slope (aa) of
      ! yf = aa * xx + b0, where ||yy - yf|| is a minimum
      IMPLICIT NONE
      real, intent(in)            :: yy(:), xx(:)
      integer, intent(in)         :: mask(:)
      real, intent(out)           :: aa
      real, intent(in)            :: b0
     
      integer :: lenArr

      lenArr = size(yy) - sum(mask(:))

      aa = ( sum(yy * (1-mask(:))) - b0 * lenArr ) / sum(xx * (1-mask(:)))

      end SUBROUTINE lin_reg1


      SUBROUTINE lin_reg( yy, xx, mask, aa, bb, amn, amx )
      ! Find the slope (aa) and y-intercept (bb) of
      ! yf = aa * xx + bb, where ||yy - yf|| is a minimum

      IMPLICIT NONE
      real, intent(in)            :: xx(:)
      real, intent(in)            :: yy(:)
      integer, intent(in)         :: mask(:)
      real, intent(in)            :: amn, amx
      real, intent(out)           :: aa, bb

      real, allocatable, save     :: x(:)
      real                        :: sumx, sumy, xave, yave, ad, ax

      integer                     :: ii, lenArr, nxy

      if ( .not. allocated(x) ) allocate(x(SPND), stat=ii)
      
      nxy = size(yy)
      
      lenArr = nxy - sum(mask(:))

      x(1:nxy) = xx(:) * (1-mask(:))
      
      sumx = sum(x(1:nxy))
      sumy = sum(yy * (1-mask(:)))
      xave = sumx / lenArr
      yave = sumy / lenArr

      ! Find the slope (aa)
      aa = dot_product(x(1:nxy),yy) - sumx*yave
      ad = dot_product(x(1:nxy),xx) - sumx*xave
      ax = amx * ad
      if ( ax > aa .and. amn * ad < aa ) then
        ! aa is between min and max
        aa = aa / ad
      else if ( ax > aa ) then
        aa = amn 
      else
        aa = amx
      end if
      
      ! y-intercept (bb)
      bb = yave - aa*xave

      end SUBROUTINE lin_reg
      

      real FUNCTION IP_calc( KsD, AM, PM )
        !
        ! The integrated, normalized pseudo-momentum of the parameterized
        ! profile has the mathematical form, in the limit of infinitesimal depth
        ! increments,
        !   IP = \int_0^KsD exp(- AM * KZ^PM ) dKZ,
        ! where D is the local depth and Ks=KsD/D is an integral wave number.
        !
        ! This can be reformulated, substituting X = AM * KZ^PM,
        !   IP = q AM^{-q} \int_0^u X^{q-1} exp(-X) dX,
        ! where q = 1/PM and u = AM * KsD^PM
        ! The integral is an incomplete gamma function gammainc(q,u), and thus
        !   IP = AM^{-q} q gammainc(q,u)
        !
        IMPLICIT NONE
        REAL, intent(in)            :: KsD, AM, PM
        real                        :: q, u

        q = 1./PM
        u = AM * KsD**PM

        IP_calc = q / AM**q * gammain(q,u)

      end FUNCTION IP_calc

      real FUNCTION gammain( q, u )

        ! The incomplete gamma function P(q,u) multiplied by Gamma(q) for
        ! u > 0 and q > 0 and q is of order one
        !
        ! Weighted average of four approximate functions depending on the water
        ! depth, or more precisely, depending on the value of xc = u / (q + 1.)
        !
        ! Let xcm < xcd < xca, xcd ~ 1.2+xcm, xca ~ 2.2+xcm:
        !
        ! 1) 0.0 <= xc < xcm : use a series representation, function gi()
        !                        (shallow and transitional water waves)
        ! 2) xcm <= xc < xcd   : Continued fraction factor, function GammaincG()
        !                        (transitional and deep water waves)
        ! 3) xcd <= xc < xca   : Gamma(q) - u**(q-1) * exp(-u)
        !                        (deep water)
        ! 4) xca <= xc         : Gamma(q)
        !                        (most of the deep ocean)
        
        IMPLICIT NONE
        real, intent(in)             :: q, u
        real                         :: numit=0.
        real                         :: xc, xce, xcm, xca, xcd
        real                         :: GmA, gammainG, wP
        ! Number of powers in series, function gi():
        integer                      :: nmP

        xc = u/(q + 1.)

        ! Shallow and transitional water thresholds
        !
        ! A requirement in Numerical Recipes function gcf(a,x) is that
        ! xc > (q-1)/(q+1). Similarly we apply an asymptotic formula
        ! gammainc -> Gamma(q) - u**(q-1) * exp(-u) for xc > xcm =
        ! max(0.4, (q-0.8)/(q+1))+0.4
        ! We apply a power series gi(u) when xc <= xcd = 2.0 + xce, where        
        xce = max ( 0., (q - 0.8)/(q + 1.) - 0.4 )
        ! The two formulas are linearly ramped in the interval xcm < xc <= xcd
        
        ! Shallow water threshold
        xcm = 0.8 + xce

        if ( xc >= xcm ) then
           ! (transitional shallow and deep water waves)
           GmA = exp( GAMMLN(q) )
        end if
        
        ! very deep water threshold
        xca = 3.0 + xce
        
        if ( xc >= xca ) then
           ! asymptotic value for xc -> infinity (very deep water)
           gammain = GmA
           ! For very deep water waves we are finished
           return
        end if
        
        ! Intermediate water
        if ( xc > xcm ) then
           gammain = GmA - u**(q-1) * exp(-u) * (xca-xc)/(xca-xcm)
           ! W = (xca-xc)/(xca-xcm); W->0 as xc -> xca; W -> 1 as xc -> xcm

           ! Intermediate to deep water threshold
           xcd = 2.0 + xce
        
           ! For moderate deep water waves we are finished
           if ( xc > xcd ) return
        end if

        ! xc <= xcd: Shallow and intermediate water
        ! Power series in u, truncated after maximally nmP terms
        nmP = int(10. * (xc + 1.)*(u/100.**(1./q)+2.)) ! Suits u in [0.1;100]
        gammainG = giGR(q, u, nmP)
        
        if ( xc > xcm ) then
           ! (transitional from shallow to deep water waves)
           wP = (xcd-xc)/(xcd-xcm)
           ! wP -> 0 as xc -> xcd; wP -> 1 as xc -> xcm
           gammain = gammainG * wP + gammain * (1. - wP)
        else
           ! xc <= xcm: Shallow water waves, asymptotic as xc -> 0
           gammain = gammainG
        end if

      end FUNCTION gammain
      
      real FUNCTION giGR(q,u,nm)
        !
        ! The incomplete gamma function, power series representation for small
        ! values of u / (q + 1.), to power nm
        !
        ! G&R: Gradshteyn and Ryzhik: Table of integrals, series and products,
        ! Academic Press 1980
        !
        ! Apply G&R 8.354.1, with al=q-1, x=u:
        !           gi(al,x) = sum_{n=0}^nm (-1)**n * frac{x^(al+n)}{n!*(al+n)}
        !                    = x^al * sum_{n=0}^nm (-1)**n * frac{x^n}{n!*(al+n)}
        !
        ! G&R 8.356.1: gi(al + 1, x) = al * gi(al,x) - x**al * exp(-x)

        IMPLICIT NONE
        REAL, intent(in)         :: q, u
        integer, intent(in)      :: nm
        integer                  :: n, ng
        real                     :: Gmq, Gf, al, dn, tn, sm, xt, tni, eps=0.001

        if ( q < 1.001 ) then
           ! 
           ! Gammainc(q,u,Gmq)  = Gmq*(1-GAMMCF(q,u))
           !                    = Gmq - EXP(-u) * u**q * GCF(q,u,Gmq/Gf,ng)
           ! Gmq == Gamma(q).
           ! GCF: Continued fraction method for the complement of the
           ! incomplete gamma function. A routine from Numerical Recipes
           ! ng: Max number of steps in GCF.
           ng = 4. + 10. * (q + 1.) / u
           Gmq = EXP( GAMMLN(q) )
           Gf = EXP(-u + q * LOG(u))
           !giGR = Gammainc(q,u,n)
           giGR = Gmq - Gf * GCF(q,u,Gmq/Gf,ng)
           return
        end if

        al = q - 1.

        dn = al
        tn = 1.
        sm = 0.
        xt = -u
        tni=tn/dn
        do n = 1,nm
           sm = sm + tni
           ! Exit loop when converged
           if ( abs(tni) < eps*sm) exit
           tn = tn * xt/n
           dn = dn + 1.
           tni=tn/dn
        end do

        ! Check convergence
        if (abs(tni) > eps*sm) then
          WRITE (NDSV, *), 'WARNING in w3xsmfmd function giGR: ', \
            'max number of terms (', nm, ') too small for q=', q, 'u=', u

          do n = nm,2*nm
            sm = sm + tni
            ! Exit loop if converged
            if ( abs(tni) < eps*sm) then
              WRITE (NDSV, *), ' - Sufficient number of terms: ', n
              exit
            end if
            tn = tn * xt/n
            dn = dn + 1.
            tni=tn/dn
          end do            
        end if

        giGR = ( al * sm  - exp(-u) ) * u**al

      end FUNCTION giGR

      real FUNCTION GAMMLN(XX)
        ! A routine from Numerical Recipes

        IMPLICIT NONE
        REAL, intent(in)             :: XX
        REAL*8 X,TMP,SER

        REAL*8, dimension(6) :: COF = (/76.18009173, -86.50532033, 24.01409822,&
                                     -1.231739516, .120858003E-2,-.536382E-5 /)
        REAL*8               :: STP = 2.50662827465D0
        REAL*8               :: HALF = 0.5D0, ONE = 1.0D0, FPF = 5.5D0

        INTEGER              :: J

        X = XX - ONE
        TMP = X + FPF
        TMP = ( X + HALF ) * LOG(TMP) - TMP
        SER = ONE
        DO J = 1,6
           X = X + ONE
           SER = SER + COF(J)/X
        END DO
        GAMMLN = TMP + LOG( STP * SER )

      END FUNCTION GAMMLN

      real FUNCTION GammaincG(A,X,NUMIT)
        ! 
        ! GammaincG(A,X,GmA) = GmA*gammainc(A,X) = GmA*(1-GAMMCF(A,X))
        !                    = Gamma(A) - EXP(-X) * X**A * GCF(A,X,NUMIT)
        ! GmA == Gamma(A). The user must supply this value.
        ! GAMMCF: Continued fraction method for the complement of the
        ! incomplete gamma function. GAMMCF is the name adpted in Numerical
        ! Recipes routines

        IMPLICIT NONE
        real, intent(in)            :: A,X
        integer, intent(in)         :: NUMIT
        real                        :: Gf,GmA
        
        GmA = exp( GAMMLN(A) )
        Gf = EXP(-X + A * LOG(X))
        GammaincG = GmA - Gf * GCF(A,X,GmA/Gf,NUMIT) 
        RETURN

      END FUNCTION GammaincG

      real FUNCTION GCF(A,X,Gof,NUMIT)
        ! Continued fraction factor in complementary incomplete gamma function,
        !    GammaincG = Gamma(A) - EXP(-X) * X**A * GCF(A,X,NUMIT)
        ! Based on:
        ! Press, W. H., S. A. Teukolsky, W. T. Vetterling, B. P. Flannery, 1993.
        ! Numerical Recipes in Fortran 77: The Art of Scientific Computing.
        ! Cambridge University Press.
        IMPLICIT NONE
        real, intent(in)            :: A,X,Gof ! Gof = Gamma(A)/(EXP(-X)*X**A)
        integer, intent(in)         :: NUMIT
        real                        :: EPS = 0.00001
        real                        :: G,GOLD,A0,A1,B0,B1
        real                        :: AN,ANA,FAC,ANF
        real                        :: FPMIN
        integer                     :: N

        G=0.
        GOLD=G
        A0=1.
        A1=X
        B0=0.
        B1=1.
        FAC=1.
        FPMIN = X * 1.e-7
        
        DO N=1,NUMIT
           AN=FLOAT(N)
           ANA=AN-A
           A0=(A1+A0*ANA)*FAC
           B0=(B1+B0*ANA)*FAC
           ANF=AN*FAC
           A1=X*A0+ANF*A1
           B1=X*B0+ANF*B1
           IF ( ABS(A1) < FPMIN ) THEN
             write(ndse,*) 'ERROR in W3XSMF FUNCTION GCF: abs(A1) < ', FPMIN
             CALL EXTCDE(1)
           end if
           FAC=1./A1
           GOLD=G
           G=B1*FAC           
           ! Test convergence of Gi = GmA - Gf * G:
           ! Gi > 0 and abs((Gi - GiOLD)/GiOLD) < eps, equivalent to:
           if ( abs(G-GOLD) < EPS*(Gof - GOLD) ) exit
        END DO
        
        GCF=G
        
        IF ( N == NUMIT .and.  xsmf_verbose > 1 ) &
           write (NDSV, *), 'WARNING in w3xsmfmd function GCF: ', &
           'number of iterations (', NUMIT, ') too small for A=', A, 'X=', X


      END FUNCTION GCF

      END MODULE W3XSMFMD
