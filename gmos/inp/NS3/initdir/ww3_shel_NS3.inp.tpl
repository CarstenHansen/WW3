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
   F       Ice concentrations
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
$ Define output data ------------ ------------------------------------- $
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
    __SIMSTART__  1200  __SIMSTOP__
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
$ CHA at GMOC 20190417 added 8 parameters HIG MXE ...
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
  T
$ MFIT  U2 (for NOEXTR==6 + 1)
  T   T   T   T   T   T   F
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
    __SIMSTART__  1200  __SIMSTOP__
$   19680606 000000    900  19680608 000000
$
    -1.620   49.695  '62059     ' $ 62059
     2.000   51.240  '62170     ' $ 62170
     2.703   51.350  'Kwintebank' $ Kwintebank
     2.439   51.389  'Westhinder' $ Westhinder
     3.370   51.393  'CadzanBoei' $ CadzandBoei
     2.770   51.410  'Akkaert   ' $ Akkaert
     3.415   51.433  'WielinNoor' $ WielingenNoord
     3.242   51.504  'Deurlo    ' $ Deurlo
     3.399   51.620  'DomburRass' $ DomburgerRassen
     3.682   51.620  'Roompotslu' $ Roompotsluis
     3.481   51.644  'Oostersc11' $ Oosterschelde11
     3.694   51.656  'Oostersch4' $ Oosterschelde4
     3.312   51.748  'Schouwenba' $ Schouwenbank, Schouwenbank2
     3.617   51.769  'Brouwersha' $ Brouwershavensegat
     3.670   51.926  'LichteiGoe' $ LichteilandGoeree1, LichteilandGoeree2
     3.000   51.948  'EurogeuDWE' $ EurogeulDWE
     3.276   51.999  'Europlatf2' $ Europlatform2, Europlatform3
     3.737   52.010  'EurogeuE13' $ EurogeulE13
     4.519   52.466  'IJgeulstr1' $ IJgeulstroompaal1
     4.268   52.493  'IJgeul5   ' $ IJgeul5
     4.058   52.550  'IJmuidMuni' $ IJmuidenMunitiestort, IJmuidenMunitiestort2
     4.151   52.926  'Q1        ' $ Q1, Q11
     4.957   52.995  'MeetPBWW11' $ MeetboeiPBW1
     1.700   53.000  '62130     ' $ 62130
     2.800   53.103  '62145     ' $ 62145
     3.220   53.218  'K13a      ' $ K13a, K13a2, K13a3
     3.633   53.267  'K141      ' $ K141
     4.662   53.277  'WaddenEGat' $ WaddenEierlandseGat
     5.576   53.320  'Amelande61' $ Amelander61
     5.102   53.322  'StortemelO' $ StortemelkOost
     4.989   53.323  'StortemelB' $ StortemelkBoei
     5.772   53.374  'Amelande62' $ Amelander62
     5.699   53.394  'Amelander5' $ Amelander51, Amelander52
     1.700   53.400  '62144     ' $ 62144
     5.636   53.425  'Amelander4' $ Amelander42
     5.604   53.452  'Amelander3' $ Amelander31, Amelander32
     6.760   53.470  'Uithuizer3' $ Uithuizerwad3
     5.574   53.478  'Amelander2' $ Amelander21, Amelander22
     5.950   53.500  'AWG       ' $ AWG
     5.482   53.509  'Amelander1' $ Amelander11, Amelander12
     6.165   53.529  'SchierWest' $ SchiermonnikoogWestgat
     6.632   53.570  'MeetbRZGN1' $ MeetboeiRZGN1
     6.167   53.596  'SchierNoor' $ SchiermonnikoogNoord
     0.700   53.600  '62150     ' $ 62150
     4.961   53.614  'L91       ' $ L91
     6.518   53.617  'MeetboWEO1' $ MeetboeiWEO1
     6.368   53.620  'MeetboWEW1' $ MeetboeiWEW1
     1.100   53.700  '62149     ' $ 62149
     2.800   53.800  '62146     ' $ 62146
     2.950   53.817  'J61       ' $ J61
     6.839   53.985  'NO1       ' $ NO1
     0.700   54.000  '62127     ' $ 62127
     1.100   54.000  '62165     ' $ 62165
     8.114   54.017  'Elbe      ' $ ElbeWR
     4.017   54.117  'F161      ' $ F161
     7.891   54.180  'Helgoland ' $ HelgolandWR
     7.818   54.219  'HelgolNort' $ Helgoland-NorthWR
     4.727   54.854  'F3platform' $ F3platform
     8.089   54.919  'Sylt      ' $ Sylt
     3.817   55.417  'A121      ' $ A121, A122
    11.223   58.251  'Brofjorden' $ BrofjordenWR
     1.909   58.371  'Sleipner-A' $ Sleipner-A
    10.932   58.488  'Vaderoarna' $ VaderoarnaWR    
$ Old known sites:
$ KDI:
    8.582    57.1315 'Hanstholm '     
    9.4105   57.577  'Hirtshals '
    9.9622   57.607  'HirtshalsN'    
    7.940    55.810  'Nymindegab'
    8.230    55.350  'FanoeBugt '
    8.060    56.470  'Fjaltring '
$
$ End of stations
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
