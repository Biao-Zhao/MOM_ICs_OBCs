#!/bin/tcsh
module load matlab/R2024a
#model resolution
set res = 3200

# hh_list can be set to (00 06 12 18)
set hh_list = (00)

set year = 2024
set months = (08)#(01 02 03 04 05 06 07 08 09 10 11)
set days = (15)

#Specify the number of days of boundary conditions to be generated
set duration = 10

foreach mon ($months)
  foreach dd ($days)

    set CDATE = ${year}${mon}${dd}

    foreach hh ($hh_list)
      set ICS  = MOM6_IC_${CDATE}${hh}_C${res}

      echo $ICS

      if (! -f $ICS) then 
         matlab -nodisplay -nosplash -r "ICS='$ICS'; run('reconstruction_current.m'); exit"  >&! stdout/MOM6_ICS_adjust_${CDATE}.log 
      else 

        echo "Skip $CDATE $hh :"
        if (-f $ICS) then
          echo " ICS already exists: $ICS, please double-check"
        endif

      endif

    end   # hh
  end     # dd
end       # mon
