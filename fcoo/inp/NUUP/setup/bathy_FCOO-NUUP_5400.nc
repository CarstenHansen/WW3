CDF       
      lat    (   lon             CDI       ?Climate Data Interface version 1.9.3 (http://mpimet.mpg.de/cdi)    Conventions       COARDS, CF-1.5     history      /�Wed Nov  1 11:56:09 2023: ncks -O -d lon,-14.50163,-13.24485 -d lat,-.92374,.96133 bathy_NUUP_5400_reg.nc bathy_NUUP_5400.nc
Tue Oct 31 13:17:14 2023: cdo -setcindexbox,-9999.,20,23,2,7 -setcindexbox,-9999.,22,40,14,30 bathy_NUUP_5400_std.nc bathy_NUUP_5400_reg.nc
Tue Oct 31 13:00:03 2023: cdo remapcon2,griddes_NUUP_5400_std.txt bathy_NUUP_1800_std.nc bathy_NUUP_5400_std.nc
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
grdconvert cut/ibcao-v4.2-200m-t16x0y0.tiff cut/ibcao-v4.2-200m-t16x0y0.nc      title         Produced by grdconvert     description       hTarget grid projection defined by: +proj=stere +lat_0=90 +lat_ts=64 +lon_0=-22 +x_0=1450000 +y_0=2430000   GMT_version       5.4.5 [64-bit]     node_offset             NCO       `netCDF Operators version 4.7.5 (Homepage = http://nco.sf.net, Code = http://github.com/nco/nco)    history_of_appended_files        hThu Dec 15 13:09:03 2022: Appended file /mnt/nfs/bifrost1-home/bjb/Projects/Bathymetry/Godthaabsfjord2022/landmask_nuuk1_100m_ibcao.nc had following "history" attribute:
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
   nco_openmp_thread_number            CDO       ?Climate Data Operators version 1.9.3 (http://mpimet.mpg.de/cdo)          
bathymetry                     	long_name         Final bathymetry at T-points   units         m      
_FillValue        �<    missing_value         �<      �  A�   lat                 standard_name         latitude   	long_name         latitude in rotated pole grid      units         degrees    axis      Y         �  R�   lon                standard_name         	longitude      	long_name         longitude in rotated pole grid     units         degrees    axis      X         l  ST   z0                     	long_name         bottom roughness   units         m      
_FillValue        �      missing_value         �      max       ?�     min                �  S��< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< D���D���C�*�B��B��aCpCR�C$�CieoC��0CwOB�>WB��cA� �A�P��< �< �< �< �< �< �< �< �< �< �< �< D�!�D���C�iB�zrB��WB��WB�f�BuE�B���CR�C��YB��>B�}�BT�A���< �< �< �< �< �< �< �< �< �< �< �< D��SD�y�C��B�mB��B��DB�T�BV�;B���C"0C�JB�#B};�BC�BZ�MB8�9B�I�B\��< �< �< �< �< �< �< �< �< Dĺ�D���Cݴ�C��B�B���B� �BVդB��NC!�3C�B��B�(�A��"A�B^�I�< �< �< �< �< �< �< �< �< �< �< D���D��[C�C߿B�+�B��B�٨Ba�B�APC
��Cg�BȐZB��.A��dA���Ap�BQ�J�< �< �< �< �< �< �< �< �< �< D��D���D-�KCI�B�x�B�G8B�XnBhYBP�B�O�C4�HC\ooB�c�Bo�Aݶ�B%��B;h�B6cR�< �< �< �< �< �< �< �< �< D��D��ED��C�,lB��AB�	B���Bq�+BL�hB��C.{�CuFB�B>!A�#�BDIA�x
�< �< �< �< �< �< �< �< �< �< D�yVD��	DA��C��B�rB�B���BvJ�BSG�B;~]C&CUFBמ�B���A�+�AG$��< �< �< �< �< �< �< �< �< �< �< D�JD��oDh��C�kBĺ�B���B���By@BZZ�B_F�C�F�C_q B�q�A�.A���A���< �< �< �< �< �< �< �< �< �< �< D���D�W�C�w�C_B��^B���By��Bg��BN�B��@C�h�C��C+��B��&A��cA���< �< �< �< �< �< �< �< �< �< �< D�lD��8C��WB�n�B�c^C��C'G�C6�C#�C�BC��WCG1_B�!C�B�A����< �< �< �< �< �< �< �< �< �< �< D�IZD�zC֠GCK��C��C�ܿC˴�C�C�phC�~�C���B�4�B���B��9B�+�B#�< �< �< �< �< �< �< �< �< �< �< D�d9D��vC�ȭC��RC��7C��+C���C�5uC�s�D ��C��4B���B@��B.N�B��4B�]/B�f�B�eC.�C�ܿ�< �< �< �< �< �< �< Dǉ{D��MDC�kC���C�cC��{CDC���C��dC��XB�~�BB��B!��B$�B��CQ2�B��'�< �< �< �< �< �< �< �< �< D���D�KD:�C�u�Ca�B�l�BM�]BS�hC=�kC���C�XcC {�BSؤB�AՑ�A�ĘB��C'T�B⡂�< �< �< �< �< �< �< �< D�2D���D���C��2B�Q�B^�zBH`
B%�|B��C�7�C���C##�B�f�B�b3B)q�A=N!A�7}B߿~B�	B���< �< �< �< �< �< �< D���D�50D�@�C�(�B��2B|4�BJ��B1fiB8�C1�C��DCwZ�B�C�B��BWɂAr��A��xB��Bt�< �< �< �< �< �< �< �< D�>D�+�D��0C���B�ǘB��B[4�B'�HB+�Bǟ0Cc��CP��BǏ&B��BB^FA-�*B\��By�LBʝC�M�< �< �< �< �< �< �< D���D��Dq,9Cـ�B���B���Bgx�B+��B(7zB��6CHQ$CX�B�mB[WMB%�Bf!iB_�Bqg>CDN�B��@�< �< �< �< �< �< �< D��D��HD���D&�QC!�7B���Bt��BH��Bb�B;�C:yCf�B��3A��B3B�T�B/�jBܖ�BV��B��< �< �< �< �< �< �< D���D��Dz}aD1h�C.yOB�X�B�JBy��B>b�BL*8B{8tC��C�'B���C8�C!��C0��C2�.B���< �< �< �< �< �< �< �< D�9�D�{eDz��D9j;C@��B�jB�n=BS��BJ��B6�B>UB��C��Bɬ�B8yA�	�An1SB�I�CYاC���< �< �< �< �< �< �< D��D�;D��nD*C%��B���B���Bs��B@��B+��B2s�BȼXCs9B�KB+�ASؐAEi=@�,b�< C%���< �< �< �< �< �< �< D��D�ʁDp�jC�1=CW�Cs�B�o�B��BK
�B1��B.�PBٖXC�>B��B8t�A��&A<H��< �< �< �< �< �< �< �< �< �< D��TDnE?D�gCzd�Cj<jCa֕CӄB�;�B�4>B4�B> �Cg�C��B��VB���B���A����< �< �< �< �< �< �< �< �< �< DgH/DM�	C�ӢCz�DC�\+C���C��aC>rB�=�BVu�B@��ChwC1UB鏵B�s�B�[�AɅe�< �< �< �< �< �< �< �< �< �< DtܱD=��C��C���C�ڈC��bC�T�C� �C|#�B���B;�5CZ�C:�>B�4�Bȸ�Bgp^A���AZ\��< �< �< �< �< �< �< �< �< D��6DBx�C� �C���C��8C�6�C�#)C��,C�|�C��C��Cr[�C���C%b�C
�lB��7B7�vA]JA@�< �< �< �< �< �< �< �< Dw��D+�5C��cC�VC��C�yC���C�<�C�C�qC�
5C�,�C���Cv�wC$E B�
�B�1�BaAP;�A@��< �< �< �< �< �< �< Dk��D9W8C���C���C��=C�]�C��C�[pC�A�C�DC���C��HC���C�eC��B��B�f>BVKA���Ap�wB"�t�< �< �< �< �< �< D<��D3�\C���C��-C�{�C�h�C���C���C�t�C���C��C�H�C�D ��C�ACskeC�B�E@BܣVA�0�A�)`�< �< �< �< �< �< D0iD�vCܮsCF;!C3C�C�AC��C#`EC:��Ck��C�YC��@C�CyC�-D KC��bCur=C9�CB��B;9`�< �< �< �< �< �< D*p�CރoC��SCC�C8B�rB�L�B�W�B]�BW�B��C���C͝:C��C���C�@bC�~C:�C�.B���BO1BB(�B=4�B u\�< �< �< C�dC�-�CSi]C<vvC��B�d�B��CB��;B�j Bv�BRt�BҐ�C�
�C�9�C�w\C@�"Co#NC|B�C{jC*b�B�PfB��B���B��Bh�S�< �< C���CaݹC>�.C&�C
��B�p^B�v?B��=B�c�B�/�B���C�C�cD ]�C�_+C�<|C/FB��C�.B�TB�D�B'z�B��vB�ѝ�< �< �< CkRCS��C0ԠC��C +�B�iB�&�B�fC`�C��C-�3CT)�C���C��D�_C�A&C��JCXv C .@C�]C=5B���B*���< �< �< �< CagUCL�C+��C�?C``C��B�QCiYC�#C 2�C5I�CE�Ct�C��aC�&LC��C�TYC���C'y�C-$CJ�C�)B�u?B}�$B�D$�< �< �< CJg�C'�C�dC
ZAC��B�KC	��CyC!�C/w�C+D�CP�C���C�̅Cٹ�C���C�C�C.��B���C�tCcD�B��:C+.��< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �< �ls��`Q�S���GU��:�@�.��"7���0�	xվ�2���t?�ȵ����Ӿ�8�|�оKud�������[􈻯��=/��=��@>��>@u�>q�b>��g>�w>�5�>��>��>?��?U?w�?+�
?86e?D��?P�?]Tv?i��?v+�hD�g@O�fzY�e�c�d�n�d(x�cb��b���a֗�a��`J��_���^���]���]2��\l��[���Z���Z��YU�X��W��W'�V=2�Uw<�T�F�S�Q�   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �< �< �< �< �   �   �   �   �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �< �< �< �< �   �   �   �   �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�< �< �< �< �   �   �   �   �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �< �< �< �< �   �   �   �   �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �< �< �< �< �   �   �   �   �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�< �< �< �< �   �   �   �   �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�< �< �< �< �< �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �< �< �< �< �< �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �< �< �< �< �< �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�< �< �< �< �< �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �< �< �< �< �< �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�< �< �< �< �< �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�< �< �< �< �< �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�< �< �< �< �< �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �< �< �< �< �< �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�< �< �< �< �< �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
�< �< �< �< �< �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �< �< �< �< �< �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �< �< �< �< �< �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �< �< �< �< �< �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �< �< �< �< �< �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �< �< �< �< �< �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�< �< �< �< �< �< �   <#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
<#�
<#�
<#�
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
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
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
<#�
<#�
<#�
<#�
<#�
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
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
<#�
�   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   �   