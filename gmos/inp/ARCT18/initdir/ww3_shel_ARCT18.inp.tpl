$ -------------------------------------------------------------------- $
$ WAVEWATCH III shel input file                                        $
$ -------------------------------------------------------------------- $
$ Define input to be used with F/T/C flag for use or nor or coupling and
$ T/F flag for definition as a homogeneous field.
$
$ Include ice and mud parameters only if IC1/2/3/4 used :
$   F F     Ice parameter 1
$   F F     Ice parameter 2
$   F F     Ice parameter 3
$   F F     Ice parameter 4
$   F F     Ice parameter 5
$   F F     Mud parameter 1
$   F F     Mud parameter 2
$   F F     Mud parameter 3
   F F     Water levels
   F F     Currents
   T F     Winds
   T       Ice concentrations
   F       Assimilation data : Mean parameters
   F       Assimilation data : 1-D spectra
   F       Assimilation data : 2-D spectra
$
$ Time frame of calculations ----------------------------------------- $
$ - Starting time in yyyymmdd hhmmss format.
$ - Ending time in yyyymmdd hhmmss format.
$
$   19680606 000000
$   19680606 060000
$
   __SIMSTART__
   __SIMSTOP__
$
$ Define output data ------------------------------------------------- $
$
$ Define output server mode. This is used only in the parallel version
$ of the model. To keep the input file consistent, it is always needed.
$ IOSTYP = 1 is generally recommended. IOSTYP > 2 may be more efficient
$ for massively parallel computations. Only IOSTYP = 0 requires a true
$ parallel file system like GPFS.
$
$    IOSTYP = 0 : No data server processes, direct access output from
$                 each process (requires true parallel file system).
$             1 : No data server process. All output for each type 
$                 performed by process that performs computations too.
$             2 : Last process is reserved for all output, and does no
$                 computing.
$             3 : Multiple dedicated output processes.
$
   1
$
$ Five output types are available (see below). All output types share
$ a similar format for the first input line:
$ - first time in yyyymmdd hhmmss format, output interval (s), and 
$   last time in yyyymmdd hhmmss format (all integers).
$ Output is disabled by setting the output interval to 0.
$
$ ------------------------------------------------------------------- $
$
$ Type 1 : Fields of mean wave parameters
$          Standard line and line with logical flags to activate output
$          fields as defined in section 2.4 of the manual. The logical
$          flags are not supplied if no output is requested. The logical
$          flags can be placed on multiple consecutive lines. However,
$          the total number and order of the logical flags is fixed.
$                               The raw data file is out_grd.ww3, 
$                               see w3iogo.ftn for additional doc.
$
    __SIMSTART__  3600  __SIMSTOP__
$   19680606 000000   3600  19680608 000000
$----------------------------------------------------------------
$ Output request flags identifying fields.
$ 
$ The table below provides a full definition of field output parameters
$  as well as flags indicating if they are available in different field 
$  output output file types (ASCII, grib, NetCDF). 
$ Further definitions are found in section 2.4 of the manual. 
$
$ Selection of field outputs may be made in two ways: 
$   F/T flags: first flag is set to F, requests made per group (1st line)
$              followed by parameter flags (total of 10 groups).
$   Namelists: first line is set to N, next line contains parameter 
$              symbol as per table below.
$ 
$   Example of F/T flag use is given in this sample ww3_shel.inp, below. 
$    For namelist usage, see the sample ww3_ounf.inp for an example.
$
$ ----------------------------------------
$ Output field parameter definitions table
$ ----------------------------------------
$
$ All parameters listed below are available in output file of the types
$  ASCII and NetCDF. If selected output file types are grads or grib, 
$  some parameters may not be available. The first two columns in the
$  table below identify such cases by flags, cols 1 (GRB) and 2 (GXO)
$  refer to grib (ww3_grib) and grads (gx_outf), respectively.
$
$ Columns 3 and 4 provide group and parameter numbers per group.
$ Columns 5, 6 and 7 provide:
$   5 - code name (internal)
$   6 - output tags (names used is ASCII file extensions, NetCDF 
$       variable names and namelist-based selection (see ww3_ounf.inp)
$   7 - Long parameter name/definition
$
$  G  G
$  R  X Grp  Param Code     Output  Parameter/Group
$  B  O Numb Numbr Name        Tag  Definition 
$  --------------------------------------------------
$        1                          Forcing Fields
$   -------------------------------------------------
$  T  T  1     1   DW         DPT   Water depth.
$  T  T  1     2   C[X,Y]     CUR   Current velocity.
$  T  T  1     3   UA         WND   Wind speed.
$  T  T  1     4   AS         AST   Air-sea temperature difference.
$  T  T  1     5   WLV        WLV   Water levels.
$  T  T  1     6   ICE        ICE   Ice concentration.
$  T  T  1     7   IBG        IBG   Iceberg-induced damping.
$  T  T  1     8   D50        D50   Median sediment grain size.
$  T  T  1     9   IC1        IC1   Ice thickness.
$  T  T  1    10   IC5        IC5   Ice flow diameter.
$   -------------------------------------------------
$        2                          Standard mean wave Parameters
$   -------------------------------------------------
$  T  T  2     1   HS         HS    Wave height.
$  T  T  2     2   WLM        LM    Mean wave length.
$  T  T  2     3   T02        T02   Mean wave period (Tm0,2).
$  T  T  2     4   T0M1       T0M1  Mean wave period (Tm0,-1).
$  T  T  2     5   T01        T01   Mean wave period (Tm0,1).
$  T  T  2     6   FP0        FP    Peak frequency.
$  T  T  2     7   THM        DIR   Mean wave direction.
$  T  T  2     8   THS        SPR   Mean directional spread.
$  T  T  2     9   THP0       DP    Peak direction.
$  T  T  2    10   HIG        HIG   Infragravity height
$  T  T  2    11   STMAXE     MXE   Max surface elev (STE)
$  T  T  2    12   STMAXD     MXES  St Dev of max surface elev (STE)
$  T  T  2    13   HMAXE      MXH   Max wave height (STE)
$  T  T  2    14   HCMAXE     MXHC  Max wave height from crest (STE)
$  T  T  2    15   HMAXD      SDMH  St Dev of MXC (STE)
$  T  T  2    16   HCMAXD     SDMHC St Dev of MXHC (STE)
$  F  T  2    17   WBT        WBT   Dominant wave breaking probability bT
$   -------------------------------------------------
$        3                          Spectral Parameters (first 5)
$   -------------------------------------------------
$  F  F  3     1   EF         EF    Wave frequency spectrum
$  F  F  3     2   TH1M       TH1M  Mean wave direction from a1,b2
$  F  F  3     3   STH1M      STH1M Directional spreading from a1,b2
$  F  F  3     4   TH2M       TH2M  Mean wave direction from a2,b2
$  F  F  3     5   STH2M      STH2M Directional spreading from a2,b2
$  F  F  3     6   WN         WN    Wavenumber array
$   -------------------------------------------------
$        4                          Spectral Partition Parameters 
$   -------------------------------------------------
$  T  T  4     1   PHS        PHS   Partitioned wave heights.
$  T  T  4     2   PTP        PTP   Partitioned peak period.
$  T  T  4     3   PLP        PLP   Partitioned peak wave length.
$  T  T  4     4   PDIR       PDIR  Partitioned mean direction.
$  T  T  4     5   PSI        PSPR  Partitioned mean directional spread.
$  T  T  4     6   PWS        PWS   Partitioned wind sea fraction.
$  T  T  4     7   PTHP0      PDP   Peak wave direction of partition.
$  T  T  4     8   PQP        PQP   Goda peakdedness parameter of partition.
$  T  T  4     9   PPE        PPE   JONSWAP peak enhancement factor of partition.
$  T  T  4    10   PGW        PGW   Gaussian frequency width of partition.
$  T  T  4    11   PSW        PSW   Spectral width of partition.
$  T  T  4    12   PTM1       PTM10 Mean wave period (m-1,0) of partition.
$  T  T  4    13   PT1        PT01  Mean wave period (m0,1) of partition.
$  T  T  4    14   PT2        PT02  Mean wave period (m0,2) of partition.
$  T  T  4    15   PEP        PEP   Peak spectral density of partition.
$  T  T  4    16   PWST       TWS   Total wind sea fraction.
$  T  T  4    17   PNR        PNR   Number of partitions.
$   -------------------------------------------------
$        5                          Atmosphere-waves layer
$   -------------------------------------------------
$  T  T  5     1   UST        UST   Friction velocity.
$  F  T  5     2   CHARN      CHA   Charnock parameter
$  F  T  5     3   CGE        CGE   Energy flux
$  F  T  5     4   PHIAW      FAW   Air-sea energy flux
$  F  T  5     5   TAUWI[X,Y] TAW   Net wave-supported stress
$  F  T  5     6   TAUWN[X,Y] TWA   Negative part of the wave-supported stress
$  F  F  5     7   WHITECAP   WCC   Whitecap coverage
$  F  F  5     8   WHITECAP   WCF   Whitecap thickness
$  F  F  5     9   WHITECAP   WCH   Mean breaking height
$  F  F  5    10   WHITECAP   WCM   Whitecap moment
$  F  F  5    11   FWS        FWS   Wind sea mean period
$   -------------------------------------------------
$        6                          Wave-ocean layer 
$   -------------------------------------------------
$  F  F  6     1   S[XX,YY,XY] SXY  Radiation stresses.
$  F  F  6     2   TAUO[X,Y]  TWO   Wave to ocean momentum flux
$  F  F  6     3   BHD        BHD   Bernoulli head (J term) 
$  F  F  6     4   PHIOC      FOC   Wave to ocean energy flux
$  F  F  6     5   TUS[X,Y]   TUS   Stokes transport
$  F  F  6     6   USS[X,Y]   USS   Surface Stokes drift
$  F  F  6     7   [PR,TP]MS  P2S   Second-order sum pressure 
$  F  F  6     8   US3D       USF   Spectrum of surface Stokes drift
$  F  F  6     9   P2SMS      P2L   Micro seism  source term
$  F  F  6    10   TAUICE     TWI   Wave to sea ice stress
$  F  F  6    11   PHICE      FIC   Wave to sea ice energy flux
$  F  F  6    12   USSP       USP   Partitioned surface Stokes drift
$   -------------------------------------------------
$        7                          Wave-bottom layer 
$   -------------------------------------------------
$  F  F  7     1   ABA        ABR   Near bottom rms amplitides.
$  F  F  7     2   UBA        UBR   Near bottom rms velocities.
$  F  F  7     3   BEDFORMS   BED   Bedforms
$  F  F  7     4   PHIBBL     FBB   Energy flux due to bottom friction 
$  F  F  7     5   TAUBBL     TBB   Momentum flux due to bottom friction
$   -------------------------------------------------
$        8                          Spectrum parameters
$   -------------------------------------------------
$  F  F  8     1   MSS[X,Y]   MSS   Mean square slopes
$  F  F  8     2   MSC[X,Y]   MSC   Spectral level at high frequency tail
$  F  F  8     3   WL02[X,Y]  WL02  East/X North/Y mean wavelength compon
$  F  F  8     4   ALPXT      AXT   Correl sea surface gradients (x,t)
$  F  F  8     5   ALPYT      AYT   Correl sea surface gradients (y,t)
$  F  F  8     6   ALPXY      AXY   Correl sea surface gradients (x,y)
$   -------------------------------------------------
$        9                          Numerical diagnostics  
$   -------------------------------------------------
$  T  T  9     1   DTDYN      DTD   Average time step in integration.
$  T  T  9     2   FCUT       FC    Cut-off frequency.
$  T  T  9     3   CFLXYMAX   CFX   Max. CFL number for spatial advection. 
$  T  T  9     4   CFLTHMAX   CFD   Max. CFL number for theta-advection. 
$  F  F  9     5   CFLKMAX    CFK   Max. CFL number for k-advection. 
$   -------------------------------------------------
$        10                         User defined          
$   -------------------------------------------------
$  F  F  10    1              U1    User defined #1. (requires coding ...)
$  F  F  10    2              U2    User defined #1. (requires coding ...)
$   -------------------------------------------------
$
$     Section 4 consist of a set of fields, index 0 = wind sea, index
$     1:NOSWLL are first NOSWLL swell fields.
$
$ Actual active parameter selection section
$
$ (1) Forcing Fields
  T
$ DPT CUR WND AST WLV ICE IBG D50 IC1 IC5
$  T   T   T   T   T   F   F   F   F   F
  T   F   T   F   F   F   F   F   F   F
$ (2) Standard mean wave Parameters
  T
$ CHA st GMOC 20190417 added 8 parameters HIG MXE ...
$ HS  LM  T02 T0M1 T01 FP DIR SPR DP HIG MXE MXES MXH MXHC SDMH SDMHC WBT
$  T   T   T   T   T   T   T   T   T   F   F   F   F   F    F    F     F
  T   F   T   T   F   T   T   T   T   F   F   F   F   F    F    F     F
$ (3) Frequency-dependent parameters
  F
$ EF TH1M STH1M TH2M STH2M WN
$  T   T   T   F   F   F
$ (4) Spectral Partition Parameters
  F
$ PHS PTP PLP PDIR PSPR PNR PDP PQP PPE PGW PSW PTM10 PT01 PT02 PEP PWS TWS 
$  T   T   T   T    T    T   T   T   T   T   T   T     T    T    T   T   T
$ (5) Atmosphere-waves layer
  F
$ UST CHA CGE FAW TAW TWA WCC WCF WCH WCM FWS
$  T   T   T   T   T   T   T   T   T   T   T
$ (6) Wave-Ocean layer
  T
$ SXY TWO BHD FOC TUS USS P2S USF P2L TWI FIC USP XSP
$  T   T   T   T   T   T   T   F   F   F   F   T   F
  F   F   F   F   T   T   F   F   F   F   F   F   T
$ (7) Wave-bottom layer
  F
$ ABR UBR BED FBB TBB
$  T   T   T   T   T
$ (8) Spectrum parameters
  F
$ MSS MSC WL02 AXT AYT AXY
$  T   T   T   T   T   T
$ (9) Numerical diagnostics
  T
$ DTD FC  CFX CFD CFK
  T   T   T   T   T
$ (10) User defined (NOEXTR flags needed)
  F
$ MFIT  U2 (for NOEXTR==6 + 1)
$  T   T   T   T   T   T   F
$
$----------------------------------------------------------------
$
$ Type 2 : Point output
$          Standard line and a number of lines identifying the 
$          longitude, latitude and name (C*10) of output points.
$          The list is closed by defining a point with the name
$          'STOPSTRING'. No point info read if no point output is
$          requested (i.e., no 'STOPSTRING' needed).
$          Example for spherical grid.
$                               The raw data file is out_pnt.ww3, 
$                               see w3iogo.ftn for additional doc.
$
$   NOTE : Spaces may be included in the name, but this is not
$          advised, because it will break the GrADS utility to 
$          plots spectra and source terms, and will make it more
$          difficult to use point names in data files.
$
    __SIMSTART__  3600  __SIMSTOP__
$   19680606 000000    900  19680608 000000
$
   -16.945   32.622  '1301000   ' $ 1301000
    -8.068   36.398  '6200200   ' $ 6200200
    -6.962   36.491  '62085     ' $ 62085
   -25.721   37.726  '6202402   ' $ 6202402
   -28.539   38.586  '6202404   ' $ 6202404
   -27.963   39.086  '6202400   ' $ 6202400
   -31.167   39.364  '6202403   ' $ 6202403
    -9.640   39.510  '6200192   ' $ 6200192
   -72.618   39.601  '44066     ' $ 44066
   -73.770   39.769  '44091     ' $ 44091
   -73.167   40.250  '44025     ' $ 44025
   -72.048   40.692  '44017     ' $ 44017
   -71.127   40.969  '44097     ' $ 44097
    -9.580   41.150  '6200191   ' $ 6200191
   -70.186   41.443  '44020     ' $ 44020
   -70.329   41.840  '44090     ' $ 44090
    -9.430   42.121  '62084     ' $ 62084
   -70.143   42.206  '44018     ' $ 44018
   -61.998   42.262  '44137     ' $ 44137
   -70.651   42.346  '44013     ' $ 44013
   -64.018   42.505  '44150     ' $ 44150
   -70.570   42.520  '44029     ' $ 44029
   -70.168   42.798  '44098     ' $ 44098
   -70.418   43.183  '44030     ' $ 44030
   -69.128   43.204  '44005     ' $ 44005
    -8.536   43.343  'Langosteir' $ Langosteira-coast-buoy
    -3.132   43.398  'Bilbao-coa' $ Bilbao-coast-buoy
    -1.677   43.403  '62079     ' $ 62079
   -67.879   43.491  '44037     ' $ 44037
    -9.210   43.496  '62083     ' $ 62083
    -1.614   43.530  '62066     ' $ 62066
   -70.144   43.531  '44007     ' $ 44007
    -2.002   43.574  'Donostia-b' $ Donostia-buoy
    -3.083   43.642  '62024     ' $ 62024
   -69.358   43.715  '44032     ' $ 44032
    -6.164   43.753  '62025     ' $ 62025
    -3.793   43.890  '6201030   ' $ 6201030
   -69.000   44.060  '44033     ' $ 44033
   -68.110   44.110  '44034     ' $ 44034
    -7.692   44.136  '62082     ' $ 62082
   -57.103   44.240  '44139     ' $ 44139
   -67.314   44.273  '44027     ' $ 44027
   -63.403   44.502  '44258     ' $ 44258
    -5.000   45.200  '62001     ' $ 62001
   -56.183   46.700  '4400050   ' $ 4400050
    -8.500   47.500  '62163     ' $ 62163
   -12.401   48.701  '62029     ' $ 62029
$ Removed all between 48N and 58N and east of 4W
   -65.710   49.538  '45138     ' $ 45138
   -10.548   51.216  '62092     ' $ 62092
    -6.704   51.690  '62094     ' $ 62094
   -15.881   53.075  'M6        ' $ 62095
   -10.146   54.231  'BelmulB-co' $ BelmulletB-coast-buoy
   -10.278   54.268  'BelmulA-co' $ BelmulletA-coast-buoy
     1.909   58.371  'Sleipner-A' $ Sleipner-A
$    10.932   58.488  'Vaderoarna' $ VaderoarnaWR
     1.500   59.500  '63110     ' $ 63110
     2.227   59.574  'Heimdal   ' $ Heimdal
     2.826   60.491  'Oseberg-A ' $ Oseberg-A
     3.719   60.644  'Troll-A   ' $ Troll-A
     1.708   61.000  '63113     ' $ 63113
     2.269   61.204  'GullfaksC ' $ GullfaksC
     1.149   61.240  'NorthCormo' $ NorthCormorant
   -20.345   63.286  'Surtseyjar' $ Surtseyjardufl
   -22.461   63.813  'Grindaviku' $ Grindavikurdufl
   -22.877   64.052  'Gardskagad' $ Gardskagadufl
   -13.627   65.648  'Kogurdufl ' $ Kogurdufl
   -24.778   65.698  'Blakksnes ' $ Blakksnes
   -18.192   66.292  'Grimseyjar' $ Grimseyjarsund
   -23.367   66.438  'Straumnesd' $ Straumnesdufl
$ Virtual Navigation stations
   -37.0     60.3    'DenmaStr01'
   -31.23    60.31   'DenmaStr02'
   -44.0     59.5    'CapeFare_1'
$ 12 nmi off Greenland Harbours 
   -50.2     59.85   'Paamiut_1 '
   -52.5     63.95   'Nuuk_1    '
   -54.5     66.9    'Sisimiut_1'
   -54.3     68.6    'Aasiaat_1 '
   -37.6     65.25   'Tasiilaq_1'
$
     0.0   0.0  'STOPSTRING'
$
$ Type 3 : Output along track.
$          Flag for formatted input file.
$                         The data files are track_i.ww3 and
$                         track_o.ww3, see w3iotr.ftn for ad. doc.
$
$ CHA at GMOC 20190417: Step '0' means this is ignored
   19680606 000000   0  19680606 013000
$     T
$
$ Type 4 : Restart files (no additional data required).
$                               The data file is restartN.ww3, see
$                               w3iors.ftn for additional doc.
$
    __SIMSTOP__  __HOTDELTA__  __SIMSTOP__
$   19680606 030000  3600  19680607 030000
$
$ Type 5 : Boundary data (no additional data required).
$                               The data file is nestN.ww3, see
$                               w3iobcmd.ftn for additional doc.
$
   __SIMSTART__   1200   __SIMSTOP__
$   19680606 000000   3600  20010102 000000
$
$ Type 6 : Separated wave field data (dummy for now).
$          First, last step IX and IY, flag for formatted file
$
$ CHA at GMOC 20190417: Step '0' means this is ignored
   19680606 000000   0  20010102 000000
$      0 999 1 0 999 1 T
$
$ Type 7 : Coupling. (must be fully commented if not used with switch COU)
$          Namelist type selection is used here.
$          Diagnostic fields to exchange. (see namcouple for more information)
$
$  19680606 000000   3600  20010102 000000
$  N
$
$   - Sent fields by ww3:
$       - Ocean model : T0M1 OCHA OHS DIR BHD TWO UBR FOC TAW TUS USS LM DRY
$       - Atmospheric model : ACHA AHS TP (or FP) FWS
$       - Ice model : IC5 TWI
$
$  CHA
$
$   - Received fields by ww3:
$       - Ocean model : SSH CUR
$       - Atmospheric model : WND
$       - Ice model : ICE IC1 IC5
$
$  WND
$
$ Homogeneous field data --------------------------------------------- $
$ Homogeneous fields can be defined by a list of lines containing an ID
$ string 'LEV' 'CUR' 'WND', date and time information (yyyymmdd
$ hhmmss), value (S.I. units), direction (current and wind, oceanogr.
$ convention degrees)) and air-sea temperature difference (degrees C).
$ 'STP' is mandatory stop string.
$ Also defined here are the speed with which the grid is moved
$ continuously, ID string 'MOV', parameters as for 'CUR'.
$
$   'LEV' 19680606 010000    1.00
$   'CUR' 19680606 073125    2.0    25.
$   'WND' 19680606 000000   20.    145.    2.0
$   'MOV' 19680606 013000    4.0    25.
   'STP'
$
$ -------------------------------------------------------------------- $
$ End of input file                                                    $
$ -------------------------------------------------------------------- $
