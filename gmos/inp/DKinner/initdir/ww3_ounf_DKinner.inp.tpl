$ -------------------------------------------------------------------- $
$ WAVEWATCH III Grid output post-processing                            $
$--------------------------------------------------------------------- $
$ First output time (yyyymmdd hhmmss), increment of output (s), 
$ and number of output times.
$
    __SIMSTART__  1200  __SIMSTOP__
$  19850101 000000  3600.  1000
$
$ Fields requested --------------------------------------------------- $
$
$ Output request flags identifying fields as in ww3_shel.inp. See that
$ file for a full documentation of field output options. Namelist type
$ selection is used here (for alternative F/T flags, see ww3_shel.inp).
$
$ DPT CUR WND AST WLV ICE IBG D50 IC1 IC5 HS LM T02 T0M1 T01 FP DIR SPR
$ DP HIG EF TH1M STH1M TH2M STH2M WN PHS PTP PLP PDIR PSPR PWS TWS PNR
$ UST CHA CGE FAW TAW TWA WCC WCF WCH WCM SXY TWO BHD FOC TUS USS P2S 
$ USF P2L TWI FIC ABR UBR BED FBB TBB MSS MSC DTD FC CFX CFD CFK MFIT U1 
$
 N
 DPT HS FP T02 T0M1 DIR TUS USS DTD MFIT
$
$--------------------------------------------------------------------- $
$ netCDF version [3,4]
$        and variable type 4 [2 = SHORT, 3 = it depends , 4 = REAL]
$ swell partitions [0 1 2 3 4 5]
$ variables in same file [T] or not [F] 
$
 3 2
 0 1 2
 T
$
$ -------------------------------------------------------------------- $
$ File prefix
$ number of characters in date [4(yearly),6(monthly),8(daily),10(hourly)]
$ - or (CHA at FCOO): max length of files [16(month)]
$ IX and IY ranges [regular:IX NX IY NY, unstructured:IP NP 1 1]
$
 ww3.         
 0
 1 143 1 153
$
$ For each field and time a new file is generated with the file name
$ ww3.date_xxx.nc , where date is a conventional time indicator with S3
$ characters, and xxx is a field identifier.
$
$ -------------------------------------------------------------------- $
$ End of input file                                                    $
$ -------------------------------------------------------------------- $
