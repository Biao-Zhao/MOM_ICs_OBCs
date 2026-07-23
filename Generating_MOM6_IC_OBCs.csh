#!/bin/tcsh
module load matlab/R2024a

set basedir = "/gpfs/f6/bil-coastal-gfdl/scratch/Biao.Zhao/MOM_ICs_OBCs"
set rundir = "/ncrc/home1/Biao.Zhao/grid_prep/MOM_ICs_OBCs/scripts"

#model resolution
set res = 3200
set NK  = 85
# hh_list can be set to (00 06 12 18)
set hh_list = (12)

set year = 2022
set months = (11)#(01 02 03 04 05 06 07 08 09 10 11 12)
set days = (28)

#Specify the number of days of boundary conditions to be generated
set duration = 3

#Reconstruct 3D geostrophic currents based on temperature and salinity, 0 no, 1 yes
set recontr = 1

foreach mon ($months)
  foreach dd ($days)

    set CDATE = ${year}${mon}${dd}
    set EDATE = `date -u -d "${CDATE} +${duration} days" +%Y%m%d`
    set sy  = `echo $CDATE | cut -c1-4`
    set sm  = `echo $CDATE | cut -c5-6`
    set sd  = `echo $CDATE | cut -c7-8`
    set ey = `echo $EDATE | cut -c1-4`
    set em = `echo $EDATE | cut -c5-6`
    set ed = `echo $EDATE | cut -c7-8` 
    set CDATE_FMT = "${sy}-${sm}-${sd}"
    set EDATE_FMT = "${ey}-${em}-${ed}"

    foreach hh ($hh_list)
      echo ${CDATE_FMT}
      echo ${EDATE_FMT}
      set ICS  = ${basedir}/ICs/C${res}/NK${NK}/MOM6_IC_${CDATE}${hh}_C${res}.nc
      set OBCS = ${basedir}/OBCs/C${res}/NK${NK}/${CDATE}${hh}
      
      echo $ICS
      echo $OBCS

      if (! -f $ICS ) then
         set mode = 2
         echo "Generating Initial Conditions: $mode" 
      #${rundir}/prepare_MOM6_inputs.sh ${CDATE_FMT} ${hh} ${EDATE_FMT} ${mode} >>& stdout/making_ICS_OBCs_${CDATE}.log 
      else
         echo " ICS already exists: $ICS, please double-check"
      endif

      if (! -d $OBCS) then
         set mode = 3
         echo "Generating Boundary Conditions: $mode"
      #${rundir}/prepare_MOM6_inputs.sh ${CDATE_FMT} ${hh} ${EDATE_FMT} ${mode} >>& stdout/making_ICS_OBCs_${CDATE}.log
      else
         echo " OBCS already exists: $OBCS, please double-check"
      endif

      if ( $recontr ) then
        set ICS_GEO  = ${basedir}/ICs/C${res}/NK${NK}/MOM6_IC_${CDATE}${hh}_C${res}_geocurrents.nc
        echo "Reconstruct geostrophic currents"
        matlab -nodisplay -nosplash -r "ICS='$ICS'; ICS_GEO='$ICS_GEO'; run('${rundir}/geostrophic_adj/reconstruction_current.m'); exit"  >>& stdout/making_ICS_OBCs_${CDATE}.log
      endif 

    end   # hh
  end     # dd
end       # mon
