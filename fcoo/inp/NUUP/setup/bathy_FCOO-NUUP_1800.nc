CDF       
      lon    +   lat    :         CDI       ?Climate Data Interface version 1.9.3 (http://mpimet.mpg.de/cdi)    Conventions       COARDS, CF-1.5     history      0�Wed Nov 29 09:52:32 2023: cdo -setcindexbox,-9999.,43,43,1,58 -setcindexbox,-9999.,1,43,58,58 bathy_FCOO-NUUP_1800.nc bathy_FCOO-NUUP_1800_missedge.nc
Wed Nov 29 09:35:59 2023: ncks -d lon,0,42 -d lat,0,57 bathy_FCOO-NUUP_1800_full.nc bathy_FCOO-NUUP_1800.nc
Wed Nov  1 11:39:14 2023: ncks -O -d lon,-14.05055,-12.61656 -d lat,-.53710,.52636 bathy_NUUP_1800_reg.nc bathy_NUUP_1800.nc
Tue Oct 31 14:14:38 2023: cdo maskregion,selreg_1800.txt -setcindexbox,-9999.,51,51,33,33 -setcindexbox,-9999.,101,112,41,60 -setcindexbox,-9999.,50,50,67,67 bathy_NUUP_1800_std.nc bathy_NUUP_1800_reg.nc
Tue Oct 31 12:58:56 2023: cdo remapcon2,griddes_NUUP_1800_std.txt bathy_NUUP_600_std.nc bathy_NUUP_1800_std.nc
Tue Oct 31 12:50:59 2023: cdo remapcon2,griddes_NUUP_600_std.txt bathy_NUUP_200_std.nc bathy_NUUP_600_std.nc
Thu Oct 26 14:42:54 2023: ncap2 -O -s grid_type=2; lon=(lon+30000.)/200.*0.001790-14.46215; lat=lat/200.*0.001790;  tmp3.nc bathy_NUUP_200_full.nc
Thu Oct 26 14:16:23 2023: ncatted -a long_name,lon,m,c,longitude in rotated pole grid -a long_name,lat,m,c,latitude in rotated pole grid -a units,lon,m,c,degrees -a units,lat,m,c,degrees -a delta_x,lon,d,, -a delta_y,lat,d,, tmp2.nc tmp3.nc
Thu Oct 26 14:08:57 2023: ncatted -a grid_mapping_name,rotated_pole,m,c,rotated_latitude_longitude -a latitude_of_projection_origin,rotated_pole,d,, -a straight_vertical_longitude_from_pole,rotated_pole,d,, -a standard_parallel,rotated_pole,d,, -a false_easting,rotated_pole,d,, -a false_northing,rotated_pole,d,, -a spatial_proj4,rotated_pole,d,, -a grid_north_pole_latitude,rotated_pole,a,d,22.2959 -a grid_north_pole_longitude,rotated_pole,a,d,160.86456 -a long_name,bathymetry,m,c,Final bathymetry at T-points tmp1.nc tmp2.nc
Thu Oct 26 13:53:18 2023: ncrename -O -d x,lon -v x,lon -d y,lat -v y,lat -v mapping,rotated_pole tmp0.nc tmp1.nc
Thu Oct 26 13:52:44 2023: ncks -v bathymetry,z0,grid_type,mapping /mnt/nfs/modeldev02-scratch/cha/Projects/NUUK200_freshwater/topo_nuuk1_200m_ibcao_hadj_gst_combo_sm02_coastmask.nc tmp0.nc
2023-03-20 11:50:59: /home/bjb/Projects/DcooWare/dcoo-bathytools-ssh2023-03-03/pytools/bathy_filt.py -A --curviloops=5 --curvistr=0.0001 --rloops=50 --rstr=0.003 --rcmax=1.9 topo_nuuk1_200m_ibcao_hadj_gst_combo_sm01_coastmask.nc.02.tmp
2023-03-20 11:50:57: /home/bjb/Projects/DcooWare/dcoo-bathytools-ssh2023-03-03/pytools/bathy_filt.py -A --curviloops=5 --curvistr=0.0001 --rloops=50 --rstr=0.003 --rcmax=1.9 topo_nuuk1_200m_ibcao_hadj_gst_combo_sm01_coastmask.nc.02.tmp
2023-03-20 11:50:55: /home/bjb/Projects/DcooWare/dcoo-bathytools-ssh2023-03-03/pytools/bathy_filt.py -A --curviloops=5 --curvistr=0.0001 --rloops=50 --rstr=0.003 --rcmax=1.9 topo_nuuk1_200m_ibcao_hadj_gst_combo_sm01_coastmask.nc.02.tmp
2023-03-20 11:50:53: /home/bjb/Projects/DcooWare/dcoo-bathytools-ssh2023-03-03/pytools/bathy_filt.py -A --curviloops=5 --curvistr=0.0001 --rloops=50 --rstr=0.003 --rcmax=1.9 topo_nuuk1_200m_ibcao_hadj_gst_combo_sm01_coastmask.nc.02.tmp
2023-03-20 11:50:50: /home/bjb/Projects/DcooWare/dcoo-bathytools-ssh2023-03-03/pytools/bathy_filt.py -A --curviloops=5 --curvistr=0.0001 --rloops=50 --rstr=0.003 --rcmax=1.9 topo_nuuk1_200m_ibcao_hadj_gst_combo_sm01_coastmask.nc.02.tmp
2023-03-20 11:49:58: /home/bjb/Projects/DcooWare/dcoo-bathytools-ssh2023-03-03/pytools/bathy_filt.py -A --curviloops=5 --curvistr=0.0001 --rloops=50 --rstr=0.003 --rcmax=1.9 topo_nuuk1_200m_ibcao_hadj_gst_combo_sm01_coastmask.nc.02.tmp
2023-03-20 11:49:56: /home/bjb/Projects/DcooWare/dcoo-bathytools-ssh2023-03-03/pytools/bathy_filt.py -A --curviloops=5 --curvistr=0.0001 --rloops=50 --rstr=0.003 --rcmax=1.9 topo_nuuk1_200m_ibcao_hadj_gst_combo_sm01_coastmask.nc.02.tmp
2023-03-20 11:49:54: /home/bjb/Projects/DcooWare/dcoo-bathytools-ssh2023-03-03/pytools/bathy_filt.py -A --curviloops=5 --curvistr=0.0001 --rloops=50 --rstr=0.003 --rcmax=1.9 topo_nuuk1_200m_ibcao_hadj_gst_combo_sm01_coastmask.nc.02.tmp
2023-03-20 11:49:51: /home/bjb/Projects/DcooWare/dcoo-bathytools-ssh2023-03-03/pytools/bathy_filt.py -A --curviloops=5 --curvistr=0.0001 --rloops=50 --rstr=0.003 --rcmax=1.9 topo_nuuk1_200m_ibcao_hadj_gst_combo_sm01_coastmask.nc.02.tmp
2023-03-20 11:49:49: /home/bjb/Projects/DcooWare/dcoo-bathytools-ssh2023-03-03/pytools/bathy_filt.py -A --curviloops=5 --curvistr=0.0001 --rloops=50 --rstr=0.003 --rcmax=1.9 topo_nuuk1_200m_ibcao_hadj_gst_combo_sm01_coastmask.nc.02.tmp
2023-03-20 11:37:24: /home/bjb/Projects/DcooWare/dcoo-bathytools-ssh2023-03-03/pytools/bathy_filt.py -A --curviloops=0 --rloops=50 --rstr=0.003 --rcmax=1.9 topo_nuuk1_200m_ibcao_hadj_gst_combo_sm01_coastmask.nc.01.tmp
2023-03-20 11:37:13: /home/bjb/Projects/DcooWare/dcoo-bathytools-ssh2023-03-03/pytools/bathy_filt.py -A --curviloops=0 --rloops=50 --rstr=0.003 --rcmax=2.0 topo_nuuk1_200m_ibcao_hadj_gst_combo_sm01_coastmask.nc.01.tmp
2023-03-20 11:36:47: /home/bjb/Projects/DcooWare/dcoo-bathytools-ssh2023-03-03/pytools/bathy_filt.py -A --curviloops=0 --rloops=50 --rstr=0.002 --rcmax=2.0 topo_nuuk1_200m_ibcao_hadj_gst_combo_sm01_coastmask.nc.01.tmp
2023-03-20 11:36:19: /home/bjb/Projects/DcooWare/dcoo-bathytools-ssh2023-03-03/pytools/bathy_filt.py -A --curviloops=0 --rloops=50 --rstr=0.002 --rcmax=2.2 topo_nuuk1_200m_ibcao_hadj_gst_combo_sm01_coastmask.nc.01.tmp
2023-03-20 11:36:00: /home/bjb/Projects/DcooWare/dcoo-bathytools-ssh2023-03-03/pytools/bathy_filt.py -A --curviloops=0 --rloops=30 --rstr=0.002 --rcmax=2.4 topo_nuuk1_200m_ibcao_hadj_gst_combo_sm01_coastmask.nc.01.tmp
2023-03-20 11:35:48: /home/bjb/Projects/DcooWare/dcoo-bathytools-ssh2023-03-03/pytools/bathy_filt.py -A --curviloops=0 --rloops=30 --rstr=0.002 --rcmax=2.6 topo_nuuk1_200m_ibcao_hadj_gst_combo_sm01_coastmask.nc.01.tmp
2023-03-20 11:35:39: /home/bjb/Projects/DcooWare/dcoo-bathytools-ssh2023-03-03/pytools/bathy_filt.py -A --curviloops=0 --rloops=30 --rstr=0.002 --rcmax=2.8 topo_nuuk1_200m_ibcao_hadj_gst_combo_sm01_coastmask.nc.01.tmp
2023-03-20 11:34:56: /home/bjb/Projects/DcooWare/dcoo-bathytools-ssh2023-03-03/pytools/bathy_filt.py -A --curviloops=0 --rloops=90 --rstr=0.001 --rcmax=3.0 topo_nuuk1_200m_ibcao_hadj_gst_combo_sm01_coastmask.nc.01.tmp
2023-03-20 11:34:41: /home/bjb/Projects/DcooWare/dcoo-bathytools-ssh2023-03-03/pytools/bathy_filt.py -A --curviloops=0 --rloops=50 --rstr=0.001 --rcmax=4.0 topo_nuuk1_200m_ibcao_hadj_gst_combo_sm01_coastmask.nc.01.tmp
2023-03-20 11:34:27: /home/bjb/Projects/DcooWare/dcoo-bathytools-ssh2023-03-03/pytools/bathy_filt.py -A --curviloops=0 --rloops=50 --rstr=0.001 --rcmax=5.0 topo_nuuk1_200m_ibcao_hadj_gst_combo_sm01_coastmask.nc.01.tmp
Wed Mar 15 13:57:06 2023: ncap2 -A -s z0[y,x]=float(-10);z0@missing_value=-10.f;z0@units="m";z0@long_name="bottom roughness";z0@min=0.f;z0@max=1.f;where(bathymetry>0){ z0=float(0.01);} topo_nuuk1_200m_ibcao_hadj_gst_combo_sm01_coastmask.nc
Wed Mar 15 13:56:27 2023: ncatted -a long_name,bathymetry,a,c,Final bathymetry at T-points -a lon_name,bathymetry,d,, topo_nuuk1_200m_ibcao_hadj_gst_combo_sm01_coastmask.nc
Wed Mar 15 13:55:39 2023: ncatted -a long_name,bathymetry,a,c,Final bathymetry at T-points -a lon_name,bathymetry,d,, topo_nuuk1_200m_ibcao_hadj_gst_combo_sm01_coastmask.nc
Wed Mar 15 13:54:37 2023: ncatted -a long_name,bathymetry,a,c,Final bathymetry at T-points -a lon_name,bathymetry,d,, topo_nuuk1_200m_ibcao_hadj_gst_combo_sm01_coastmask.nc
Wed Mar 15 13:54:32 2023: ncatted -a long_name,bathymetry,a,c,Final bathymetry at T-points -a lon_name,bathymetry,d,, topo_nuuk1_200m_ibcao_hadj_gst_combo_sm01_coastmask.nc
2023-03-15 09:39:50: /home/bjb/Projects/DcooWare/dcoo-bathytools-ssh2022-11-16/Scripts/remove_lakes.py --seed=10,20 -A topo_nuuk1_200m_ibcao_hadj_gst_combo_sm01_coastmask.nc
Wed Mar 15 09:35:19 2023: ncap2 -A -S nuuk1_200m_addcells.nco topo_nuuk1_200m_ibcao_hadj_gst_combo_sm01_coastmask.nc
Wed Mar 15 09:35:19 2023: ncks -O -3 --no_abc --header_pad=10000 /mnt/nfs/bifrost1-home/bjb/Projects/Bathymetry/Godthaabsfjord2022/topo_nuuk1_200m_ibcao_hadj_gst_combo_sm01.nc topo_nuuk1_200m_ibcao_hadj_gst_combo_sm01_coastmask.nc
2023-03-14 11:38:25: llc2cart.py -A /mnt/nfs/bifrost1-home/bjb/Projects/Bathymetry/Godthaabsfjord2022/topo_nuuk1_200m_ibcao_hadj_gst_combo_sm01.nc
Tue Mar 14 11:38:25 2023: ncap2 -s x@delta_x=200.0d;y@delta_y=200.0d /mnt/nfs/bifrost1-home/bjb/Projects/Bathymetry/Godthaabsfjord2022/topo_nuuk1_200m_ibcao_hadj_gst_combo_sm01.nc
2023-03-14 11:38:10: /home/bjb/Projects/DcooWare/dcoo-bathytools-ssh2022-11-16/pytools/bathy_coarsen.py -O --x0=0 --dx=200 --y0=0 --mask=0.51 /mnt/nfs/bifrost1-home/bjb/Projects/Bathymetry/Godthaabsfjord2022/topo_nuuk1_100m_ibcao_hadj_gst_combo_sm01.nc /mnt/nfs/bifrost1-home/bjb/Projects/Bathymetry/Godthaabsfjord2022/topo_nuuk1_200m_ibcao_hadj_gst_combo_sm01.nc
Thu Dec 15 13:09:17 2022: ncks -O --no-alpha -x -v landmask /mnt/nfs/bifrost1-home/bjb/Projects/Bathymetry/Godthaabsfjord2022/topo_nuuk1_100m_ibcao_hadj_gst_combo_sm01.nc /mnt/nfs/bifrost1-home/bjb/Projects/Bathymetry/Godthaabsfjord2022/topo_nuuk1_100m_ibcao_hadj_gst_combo_sm01.nc-tmpnoll
Thu Dec 15 13:09:10 2022: ncap2 -s where(landmask>0.5) bathymetry=-9999 /mnt/nfs/bifrost1-home/bjb/Projects/Bathymetry/Godthaabsfjord2022/topo_nuuk1_100m_ibcao_hadj_gst_combo_sm01.nc
Thu Dec 15 13:09:03 2022: ncks -A /mnt/nfs/bifrost1-home/bjb/Projects/Bathymetry/Godthaabsfjord2022/landmask_nuuk1_100m_ibcao.nc /mnt/nfs/bifrost1-home/bjb/Projects/Bathymetry/Godthaabsfjord2022/topo_nuuk1_100m_ibcao_hadj_gst_combo_sm01.nc
Thu Dec 15 13:07:20 2022: ncks -O -v landmask /mnt/nfs/bifrost1-home/bjb/Projects/Bathymetry/Godthaabsfjord2022/topo_nuuk1_100m_ibcao.nc /mnt/nfs/bifrost1-home/bjb/Projects/Bathymetry/Godthaabsfjord2022/landmask_nuuk1_100m_ibcao.nc
Thu Dec 15 13:07:14 2022: ncap2 -s landmask[y,x]=float(0) -s where (bathymetry==-9999) landmask=1 -A /mnt/nfs/bifrost1-home/bjb/Projects/Bathymetry/Godthaabsfjord2022/topo_nuuk1_100m_ibcao.nc
2022-12-15 13:06:32: /home/bjb/Projects/DcooWare/dcoo-bathytools-ssh2022-11-16/Scripts/remove_lakes.py --seed=50,50 -A /mnt/nfs/bifrost1-home/bjb/Projects/Bathymetry/Godthaabsfjord2022/topo_nuuk1_100m_ibcao.nc
Thu Dec 15 13:06:07 2022: ncap2 -s landmask[y,x]=float(0) -s where (bathymetry==-9999) landmask=1 -A /mnt/nfs/bifrost1-home/bjb/Projects/Bathymetry/Godthaabsfjord2022/topo_nuuk1_100m_ibcao.nc
Thu Dec 15 12:44:58 2022: ncap2 -s where (tmpmask < 0.5 && bathymetry < 0) bathymetry = -9999; /mnt/nfs/bifrost1-home/bjb/Projects/Bathymetry/Godthaabsfjord2022/topo_nuuk1_100m_ibcao.nc
Thu Dec 15 12:40:57 2022: ncap2 -s tmpmask[y,x]=float(0); -s tmpmask@missing_value=-9999. /mnt/nfs/bifrost1-home/bjb/Projects/Bathymetry/Godthaabsfjord2022/topo_nuuk1_100m_ibcao.nc
2022-12-15 12:40:51: llc2cart.py --noconv --noxgrid -A /mnt/nfs/bifrost1-home/bjb/Projects/Bathymetry/Godthaabsfjord2022/topo_nuuk1_100m_ibcao.nc
2022-12-15 12:38:46: bathy2cart.py -O --dx=100 --tgtproj=nuuk1 /mnt/nfs/bifrost1-home/bjb/Projects/Bathymetry/Godthaabsfjord2022/topo_GODTHAABSFJORD_100m.nc /mnt/nfs/bifrost1-home/bjb/Projects/Bathymetry/Godthaabsfjord2022/topo_nuuk1_100m_ibcao.nc
2022-11-15 15:31:41: bathy_multres.py -O --bat=bathymetry --mult=2 /mnt/nfs/modeldev02-scratch/valid/OBSERVATIONS/ORIGINAL_DATA/INT_IBCAO/topo_GODTHAABSFJORD.nc topo_GODTHAABSFJORD_100m.nc
Tue Nov  8 15:24:30 2022: ncap2 -O -s grid_type=1;mapping=3996 -s grid_type@long_name="Type of horizontal grid" -s grid_type@option_1_="Cartesian" -s grid_type@option_2_="Spherical" -s grid_type@option_3_="Curvilinear" -s grid_type@option_4_="Spherical Curvilinear" -s mapping@grid_mapping_name = "polar_stereographic" -s mapping@epsg = 3996 -s mapping@init = "epsg:3996" -s mapping@latitude_of_projection_origin = 90. -s mapping@straight_vertical_longitude_from_pole = 0. -s mapping@standard_parallel = 75. -s mapping@false_easting = 0. -s mapping@false_northing = 0. -s mapping@spatial_proj4 = "+proj=stere +lat_0=90 +lat_ts=75 +lon_0=0" topo_GODTHAABSFJORD.nc.tmp2 topo_GODTHAABSFJORD.nc
Tue Nov  8 15:23:55 2022: ncatted -A -a _FillValue,bathymetry,o,f,-9999. -a missing_value,bathymetry,o,f,-9999 topo_GODTHAABSFJORD.nc.tmp2
Tue Nov  8 15:23:55 2022: ncks -O -x -v z topo_GODTHAABSFJORD.nc.tmp1 topo_GODTHAABSFJORD.nc.tmp2
Tue Nov  8 15:21:10 2022: ncap2 -O -s bathymetry=-z GODTHAABSFJORD.nc topo_GODTHAABSFJORD.nc.tmp1
Tue Oct 11 10:23:22 2022: ncks -O -d x,-2400000.,-2100000. -d y,-1950000.,-1650000. ibcao_v4_2_200m_t16x0y0.nc GODTHAABSFJORD.nc
grdconvert cut/ibcao-v4.2-200m-t16x0y0.tiff cut/ibcao-v4.2-200m-t16x0y0.nc   title         Produced by grdconvert     description       hTarget grid projection defined by: +proj=stere +lat_0=90 +lat_ts=64 +lon_0=-22 +x_0=1450000 +y_0=2430000   GMT_version       5.4.5 [64-bit]     node_offset             NCO       _netCDF Operators version 4.7.5 (Homepage = http://nco.sf.net, Code = http://github.com/nco/nco)    history_of_appended_files        hThu Dec 15 13:09:03 2022: Appended file /mnt/nfs/bifrost1-home/bjb/Projects/Bathymetry/Godthaabsfjord2022/landmask_nuuk1_100m_ibcao.nc had following "history" attribute:
Thu Dec 15 13:07:20 2022: ncks -O -v landmask /mnt/nfs/bifrost1-home/bjb/Projects/Bathymetry/Godthaabsfjord2022/topo_nuuk1_100m_ibcao.nc /mnt/nfs/bifrost1-home/bjb/Projects/Bathymetry/Godthaabsfjord2022/landmask_nuuk1_100m_ibcao.nc
Thu Dec 15 13:07:14 2022: ncap2 -s landmask[y,x]=float(0) -s where (bathymetry==-9999) landmask=1 -A /mnt/nfs/bifrost1-home/bjb/Projects/Bathymetry/Godthaabsfjord2022/topo_nuuk1_100m_ibcao.nc
2022-12-15 13:06:32: /home/bjb/Projects/DcooWare/dcoo-bathytools-ssh2022-11-16/Scripts/remove_lakes.py --seed=50,50 -A /mnt/nfs/bifrost1-home/bjb/Projects/Bathymetry/Godthaabsfjord2022/topo_nuuk1_100m_ibcao.nc
Thu Dec 15 13:06:07 2022: ncap2 -s landmask[y,x]=float(0) -s where (bathymetry==-9999) landmask=1 -A /mnt/nfs/bifrost1-home/bjb/Projects/Bathymetry/Godthaabsfjord2022/topo_nuuk1_100m_ibcao.nc
Thu Dec 15 12:44:58 2022: ncap2 -s where (tmpmask < 0.5 && bathymetry < 0) bathymetry = -9999; /mnt/nfs/bifrost1-home/bjb/Projects/Bathymetry/Godthaabsfjord2022/topo_nuuk1_100m_ibcao.nc
Thu Dec 15 12:40:57 2022: ncap2 -s tmpmask[y,x]=float(0); -s tmpmask@missing_value=-9999. /mnt/nfs/bifrost1-home/bjb/Projects/Bathymetry/Godthaabsfjord2022/topo_nuuk1_100m_ibcao.nc
2022-12-15 12:40:51: llc2cart.py --noconv --noxgrid -A /mnt/nfs/bifrost1-home/bjb/Projects/Bathymetry/Godthaabsfjord2022/topo_nuuk1_100m_ibcao.nc
2022-12-15 12:38:46: bathy2cart.py -O --dx=100 --tgtproj=nuuk1 /mnt/nfs/bifrost1-home/bjb/Projects/Bathymetry/Godthaabsfjord2022/topo_GODTHAABSFJORD_100m.nc /mnt/nfs/bifrost1-home/bjb/Projects/Bathymetry/Godthaabsfjord2022/topo_nuuk1_100m_ibcao.nc
2022-11-15 15:31:41: bathy_multres.py -O --bat=bathymetry --mult=2 /mnt/nfs/modeldev02-scratch/valid/OBSERVATIONS/ORIGINAL_DATA/INT_IBCAO/topo_GODTHAABSFJORD.nc topo_GODTHAABSFJORD_100m.nc
Tue Nov  8 15:24:30 2022: ncap2 -O -s grid_type=1;mapping=3996 -s grid_type@long_name="Type of horizontal grid" -s grid_type@option_1_="Cartesian" -s grid_type@option_2_="Spherical" -s grid_type@option_3_="Curvilinear" -s grid_type@option_4_="Spherical Curvilinear" -s mapping@grid_mapping_name = "polar_stereographic" -s mapping@epsg = 3996 -s mapping@init = "epsg:3996" -s mapping@latitude_of_projection_origin = 90. -s mapping@straight_vertical_longitude_from_pole = 0. -s mapping@standard_parallel = 75. -s mapping@false_easting = 0. -s mapping@false_northing = 0. -s mapping@spatial_proj4 = "+proj=stere +lat_0=90 +lat_ts=75 +lon_0=0" topo_GODTHAABSFJORD.nc.tmp2 topo_GODTHAABSFJORD.nc
Tue Nov  8 15:23:55 2022: ncatted -A -a _FillValue,bathymetry,o,f,-9999. -a missing_value,bathymetry,o,f,-9999 topo_GODTHAABSFJORD.nc.tmp2
Tue Nov  8 15:23:55 2022: ncks -O -x -v z topo_GODTHAABSFJORD.nc.tmp1 topo_GODTHAABSFJORD.nc.tmp2
Tue Nov  8 15:21:10 2022: ncap2 -O -s bathymetry=-z GODTHAABSFJORD.nc topo_GODTHAABSFJORD.nc.tmp1
Tue Oct 11 10:23:22 2022: ncks -O -d x,-2400000.,-2100000. -d y,-1950000.,-1650000. ibcao_v4_2_200m_t16x0y0.nc GODTHAABSFJORD.nc
grdconvert cut/ibcao-v4.2-200m-t16x0y0.tiff cut/ibcao-v4.2-200m-t16x0y0.nc
   nco_openmp_thread_number            CDO       ?Climate Data Operators version 1.9.3 (http://mpimet.mpg.de/cdo)          lon                 standard_name         	longitude      	long_name         longitude in rotated pole grid     units         degrees    axis      X         �  B�   lat                standard_name         latitude   	long_name         latitude in rotated pole grid      units         degrees    axis      Y         �  CP   
bathymetry                     	long_name         Final bathymetry at T-points   units         m      
_FillValue        �<    missing_value         �<      &�  D8   z0                     	long_name         bottom roughness   units         m      
_FillValue        �      missing_value         �      max       ?�     min                &�  k0�`Υ�`���`J��`��_Ʋ�_���_B��_ ��^���^|��^:��]���]���]t��]2��\���\���\l��\*��[���[���[d��["��Z���Z���Z\��Z��Y��Y��YU�Y�X��X��XM�X�W��W�!�WE$�W'�V�+�V.�V=2�U�5�	xտY�9C��2����b��о�t?��4�����ȵ���u���6e���Ӿ��A��w���8���������|�оls��[�Kud�:�@�*w����	xս��b������ӽ����[������ӻ���</��<��=/��=q�b=���=��@=��=���>��>w�>/��>@u�>P�>at?>q�b>�9C>�x�>��g>���>�7�>�w>���>��@>�5��< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< BR]�BF�KB/XB|6�CB��C�~HC���C��B���B�\PB�o}BĉBw`�B8v�A��&A�8�A��GA>��A@��< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< BS!�BE�BC�B�ķC���C��C�z�COlEB�iB늱B�'B��(B
ڲA�$nAȌ�A�9�ApשA7}�A@��< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< BJ�:B8��Bqz�C\�uC���C���C�P�Ce�$C� CMB��B:ٖB��A���A��A���Ad`A@%N�< A�	�< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< BG-�B;��B«�C�jC�l}C��!Cs��CE�C:ߟCN�(C9��B�2�B8��A�
+A��A��SA�0�< �< A��A�@A���< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< B?0@B6��C��C��C�h�C���C,=�B�VB�}�C30�C_C��B��cBh�?A��)A��A�H��< �< A���< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< B2��B;��CC�C�3wC��YC�шCH�MB�B���C9�C/�Ck%�C)�'B��jBqAǅ�A�R��< �< A�2�< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< B5®B�ReC�i�C�  C�m�C���C|U]B�t�BۢACifBܻ�C6^�C^�<B�6{BD��A�("A�v"A�36�< A�w�< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< B��,C4�@C�C�C��bC��KD�}C��C��C�B俞B��B���CMX*C9AB�B�BJBA�4$A��A�jA���< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< C=�C�8�C��C�� DD��C�{C5�C��B��B{�SBÖ�C��C>�:CE$Bv�pA�>�A��aA�:��< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< C���C���C���C���D��D�C��GC�B�B�KB�.B�<(B��NC�B裠B�K	A�;A�UA����< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< C�w	C�D�C�ɳC�*D�C�V�C@e�B���B�0yBYGRBzIgBq,[B�A�B�ppB�9�C	fMB+�]Bl0�B(U$B�^A��A�ȢA���< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< CЇ�C՟pC�%>C�5fC���C*�kB��B��B��eBr�.B�4	B"�B�{B�`CpNCU�SB�5�B@��B�}�< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< C���C�բDt�C��EC�5yB�-�B��|B��B�kB��&C� B�joBLVB@�sBj�.B��.B��rB"�Bf��< �< �< �< �< �< �< �< �< �< BQ/�< �< �< �< �< �< �< �< �< �< �< �< �< C�	DH^D�C��xC�<�B��B��B�L&B�YPBalB�%{B&��A�t�B��B_�BۍB�=�B��cB�sB
���< �< �< �< �< �< B��CA�C���C'��< �< �< �< �< �< �< �< �< �< �< �< �< D�|D{UD��CܭJC5�DBԮ�B��B�gBFBpBr��BA�02B�B���B�	LBx*B��B���B�w�B~'Bx{B	{�B��B�;C��C�AB��B��C�YC��C��E�< �< �< �< �< Bq &�< �< �< �< �< C���C�,�C䅹C�tC]E�B��SB��B�gLB35�BwOBf�BbA�=�B:�pB���B���B���B��5B��C�hC*�C$�C�.B�i�B��B3��< �< B"�< C
��C�X�CearClaCl�ZC.=�C&2�C"�]C��B�j�A���< �< C݋�C��nC�^D�SC��>C"B��B�H�B/4dB��B�sB�.B boBHeyB\�B]�B=̓ByhB^
�C�C[��C�>�CqdCޱBEZ��< �< �< �< �< �< �< B)4aBތkB��sB*���< �< �< �< �< �< �< C�ҡCĐjC��gC�)�C�9kCH��C3�B�\Bh�BmpB_A�B��A�}B9�>B@T�B'��A��AȺ�B$�C�CMCp�@C}9kCr�B�-y�< B����< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< C~��C�ԜC�>�C�b�C��C(tC%0�B��B�|B���BP�B��A�޿B
��A�+A��AȤA��B�B�,BʔOC��C4�CA�C,%IB��C-B��U�< �< �< �< �< �< �< �< �< �< �< �< �< �< �< Cc1�C��#C��C̍�C���C0��C	��B��?Bc�BH<�B.�aB:� B6�A�I2AȦ�A�/-A��A��A�XA�D�BX�B� �C$`�C,QOB� �C	�PC;CoC��B�b�< �< �< �< �< �< �< �< �< �< �< �< �< �< Ct�C�S=C���C�q�C���C7:xC%w~C-1�B�O�B[SBBCY2BZ�YB.t�Bf-A��A�$�A�z�A��A�!+A�&�A�L�B#�XB�6ICL|-C8}C+VC\�FCT�B���< �< �< �< �< �< �< �< �< �< �< �< �< �< Cq�C���C�>YC͖PC��kC$3@C2��B�BB�M�B��B?��Bc�lBDi�B8$�B#�#B�7Aڙ�A�s�A���A���Aą�A�'B8ܦCG�Ce=�C#p�B�o�B��B���B��< �< �< �< �< �< �< �< �< �< �< �< �< Cv�oC�,�C�x�C�)�C�dUC>��C0�hB�C	"NB�fBX$~B��BZ��BJBIx�B@P"A��Ata�A,�OA)qcA"X��< �< BU��C��B��\A�dB�[�Bx'2�< B�xB�b�< �< �< �< �< �< �< �< �< �< �< CnljC�5`C���C��PCσ<C��IC?�C(��B�uYC"�C
��B��9B���Be�MBd,hBg\B,�AkA@mA @�A �A !A aA'�B���C@jBy��Bo��B� ZB��B�B�6�< �< �< �< �< �< �< �< �< �< �< CL��C�C���C�r1C�)C�ڬC�'\CE��B� �B�Bi��B�pB�a�B���B���B�12BPVaB"��A�LA$�SA &A �A IA $�A��C;!)C-�1B��^B�B��B��< �< �< �< �< �< �< �< �< �< �< �< Ce-CdO1C���C��C�^C��C��FC?3B��B��\B���Bwc�B��vB��B�ŐB��BP��B@��A���A>��A MA <A �A|��B���B��aB�wwB4��< B���< �< �< BFU�Bd��< �< �< �< �< �< �< �< BшcC,�eCl#�Cy��C�ԡC�6�C��Ce8CA�C�C|�BwAB�* C �.B���BuO�BM
B*�rA��A>��A �A �A �A���Bʰ�BUܽA�
�A (j�< �< �< �< �< B�#�< �< �< �< �< �< �< �< �< B���C"�C?Y�CQB�C���C��]C�	C�W'C&-ZC�{C֭B�p_B���C+��B���B��XB`�JA�k�Abv.A!� A f�< A!B6Z�C�B�R��< �< �< �< �< �< B閄B�D��< �< �< �< C��C�YC��W�< �< B�'�B���C2x�CP�"Cp!�C�X�C��:C���C�JCKB���B�5=B�BҠ�Bћ�B�`�B|a�Aȩ�A ��A _A A Ad�`B��XB�i3B9:$�< �< �< �< �< C ZTClk�C;O�C|��C���C��ND�D(&�D-}�D���< �< B0a1B��C�{CE�aCc^�C�[�C�6gCDP<B�NrB���B�dMB��2B���B���B���B�#�BqZ�A��A ��A �A 3�A=�B���B�zBYrB,w�< �< �< �< C?"C��C멌D]D
�=D��C�C��lC���C�X�C2qh�< �< B$�/B��B��C3�]C>4�CU�8C[��C��CjEB��*B��'B4�|Be*�Bn�B��lB�w�BEB1A�cA �aA �A���B��NBҪ^A٢�B45�Bu\B�dhB�bB��CgɑC���Cv!/CRқCY\�Cf+�CnCB�1��< �< �< �< �< �< B+8�B}~�BƵC"x�CGt�C_�~CRB#B٪XB��B�6C�%B�>�BFR�B,��B<�]B?,B��B ýA�4HA�!A�HgC7t�B��BZ�GB�zLBph�B$��BI�2C)�3C�vqC��B�t��< �< �< �< �< �< �< �< �< �< �< BDu�By��B��ZC	�C`-sC��C1�#B�2B�ZqB��Bn2B�*�B���B`)wBHE�B?�A��pAؤ�A棵B1�xC��B݅IB
XA��B�A�0�BF{�CZ�C��C]0�Bz�Y�< �< �< �< �< �< �< �< �< �< �< �< Bb	�B�HB�saB�ECW�C�W�CcAB�\B���B��
B��BB��BNUoBi�:Bw�RBv�B/�pBF��B]Z�B�kwB��GBL7A<B�Aj�gAW�&A�] C5؛CbƗC��B���Av���< �< �< �< A�^A���< �< B�!�CmO��< �< Bb+BR�2BN�B��2C2�	C~��C���CO��B���B�!hB:��BvJ�B��A���BSmBU�B.��B��MB��B˒)B���A�`Ad��Ap���< C�Cu.�C��BzP�B>��< �< �< �< �< B�N0B����< �< C�#BC�P��< �< BM�B-�(B%�Bf%�B��iCD'�C�M�C���C6$�B�l�BR4]B#�AA|>>AI��AM�\A��CB���B���B�UB���B��,Bj�PBV2tBD�WB�TJB��fB�S�B��BRA��>�< �< �< B�WLB�c��< �< C�XtC�y5CU�z�< �< B5�B"�%B&�B]�B��rC�CXF�C�+C��C#��B{�0B�.A��Ai-�A XA �A$I�A�73A�x�A�{B %�B-ZNB$�mBfN�B���B��UBM�A�uzA�� B��BzB)��B�B1�B\���< �< �< C��C�0C�e�< �< B+HB3�)BG�BoИBwB�B���CCQ�7Cs%�CE?�B���B}ߊB#k�BXKA�@�A$��A]f�B	kA�$A;��A�݂B�B��B��Bǅ;B)�;A��BecB��/Bz%k�< �< Ba��B�)��< �< �< �< C��*CĐ�B��w�< �< BE�FBY�Bi�^Be=]B~�B���B���B�YC!qC^�C ��B��DB�:�B�x0B�e�B�r�C�CT~eC:�JC��C$� CN�C^�C��C��C$�B���B�1�CJ	B�b��< �< �< �< �< C6w-C��FCǿ'D_ODw�C�~��< �< B_�`BT��BI4JB8V�BO_^B��wB�%B�	�B�\CC/�kC OYC%�1C:�@C;	�C`VaC8�CP��C��C���C�0TC�;C�f_C`29CyVC��C��C�ԭC�(Bbn=A�>s�< �< �< �< B�C�C�
\D�gCɣ�C�.8C���C��< �< B�R�B]�vB5x�B�=B�,BI�Bǅ�B���CRC�C �B��C	�B� �B��4B6� BH�zB!l�B�B�IB^5�Ap'9A�]A��=B�yB��Cz."C{�MC:�PC�BB��B�/�C=5:C���C��CҺ�C��B�?Z�< �< CIo�< �< B4�B&��B�B�9B��B&��B���B�ѰB�'pB�!CtoCB���B��8B�&�B.�B}Af�+A g�A9?{A����< �< �< �< A���Bn=C:5C���C�J�C�,TC���C�gXCF��B�\C�< BI�A�q�A��A�t�< �< �< B:֑B�BF�B{)B"@�B-|�B�~&B���BԊ�B��JB�@�C?T�B��rB���B��B)�BS�A�fA ��A /0A���< �< @�Vx@��A	�A���B��B�OC���C�ZC�P�C��:C+swB�6��< �< �< �< �< �< �< �< B2�'B&�vB#�Bc�B(lB9}]B�OB�=}B���B��dC	��C'��C ��B��KB���B`�oB��A��A �A#��A��Ao���< @g@��@�Y��< �< �< B��Cf�C�neC��C��C{��C!�B�ypB�K�< �< �< �< �< B/s�B-�B,��B(��B/S�BDHB��B�[@B���B��C&b�C	�B�B��lB���B���B+�A���A%�pA h�A�\Afp��< �< �< @���< �< �< �< �< A�~zBտSC�D&C���C���CSaRCA�UC6/�C/tP�< �< �< B-�>B*�B)�B*LB0�uBO�6B�B�HB�ƑCW�C��C�B�UB�vB�#B��B8~A���ADm�A �BA�<_Ay\��< �< �< �< �< �< �< �< �< �< Av�
Ba�YCa/C��CJ��B���C$6tCW/>CB��< �< B1۬B.e�B-MB%��B%�lBL��B��1B�l�B��C τC�B��B�GqB�7�B�T�B�3�B8AA�JGAS{gA!�BA���A}\5A GK�< �< �< �< �< �< �< �< �< �< �< B��2C��C^��B 5�B+ �BŮB�Q��< �< B@�TB4B�B��B��B|�BKf�B���B�v3Co�C%?Cn�B�&,C��CV/B��dB}�rB0 jA���AB�]A �A���A"�A %%�< �< �< �< �< �< �< �< �< �< �< �< C,�KC�+�CM��C��CRn'B����< �< BK��B2��BY�B=�B�1Bf��B��B�_C 9CK@C
۞B�tEB��LC�}B��B�?6B5?�A��A%�@A �B�Ah�TA �< �< �< �< �< �< �< �< �< �< �< �< B
�<C�V�D�)D�uC��Cq-�< �< BB�sB6[�B$cOB�dBh�B���B��"C�iC2)�CNN�C��B�.C�(C	tB�!,B��)Br{	B�A���A�N�BF�A%�l�< �< �< �< �< �< �< �< �< �< �< �< �< �< B��C`\�C�D�pD���< �< B;r�B6p�B(/Bi�BW�B�HB�H�C'CbQCIUaCI>B��C��C5PC��C��B�w=Bn�,Bnv�B�B�_�A��uAFq��< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< C*��C�z��< �< B@�B2��B.?�B|�B�bB��B�y#C*��C8��C�wB�|�B��B�-B�l�B�$�B���B�[%B�?CB�#8B��,B��B+{�A.��A~O7�< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< B`B<�B+�3B��B'�jBkz�B�}�C)giC>�UC(fC*_Ci�C�B��RB���B���B�PgB�3?B�/B�9;B���B1�1A��aB!U��< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< B��uBY��B6��B,F�B3�B\��B���C!�CJO?C3p�CɣB�(.C��B�EB��B}.cB�AB��B�N(B�]kB��>A���A#"A��	A%���< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< B��B^O B=�9B<��BC=�BW2B�]C0�CZVpC��C�'C5�B���B��B��B��B�HVBs�eB���B���BZ�&A�̑A"��A��JA'մ�< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< B�aKBv�BK��B2�B8�B<�B��C�=CZJ CB�yC`�CuEBΖ�B�:�C�)C
!�Bֿ�BmnuB���Bi\�B24A���A:�~A*f�Ab�A*r��< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   <#�
�   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   <#�
<#�
<#�
�   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   <#�
�   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   <#�
�   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   <#�
�   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �< �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   �   �   �   �   �   <#�
�   �   �   �   �   �   �   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   �   <#�
<#�
<#�
<#�
�   �   �   �   �   �   �   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   <#�
�   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   <#�
�   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   �   �   <#�
<#�
<#�
<#�
�   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   <#�
�   �   �   �   �   �   �   �   �   �   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   �   �   �   �   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   �   �   �   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   �   �   �   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   �   �   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   <#�
<#�
<#�
<#�
<#�
<#�
�   <#�
<#�
�   �   �   �   �   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   �   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   <#�
�   �   �   <#�
<#�
�   �   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   <#�
�   �   �   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   <#�
<#�
<#�
<#�
�   �   �   �   �   �   <#�
<#�
�   �   �   �   <#�
<#�
<#�
�< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   �   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   <#�
<#�
�   �   <#�
<#�
�< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   <#�
<#�
<#�
<#�
<#�
�   �   �   �   �   <#�
<#�
�   �   <#�
<#�
�< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   <#�
<#�
�   �   <#�
<#�
<#�
�< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   <#�
<#�
<#�
�< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   <#�
<#�
�   �   �   �   <#�
<#�
<#�
�< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   <#�
<#�
<#�
<#�
<#�
<#�
�< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
�< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   <#�
�< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   <#�
<#�
<#�
<#�
�   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�< <#�
<#�
<#�
�   �   �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   <#�
�   �   �   �   �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   �   �   �   �   �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   �   �   �   �   �   �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
�< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   �   �   �   �   �   �   �   <#�
<#�
<#�
<#�
<#�
<#�
�< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   �   �   �   �   �   �   �   <#�
<#�
<#�
<#�
<#�
<#�
�< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   �   �   �   �   �   �   �   �   �   <#�
<#�
<#�
<#�
<#�
�< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   <#�
<#�
�< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< 