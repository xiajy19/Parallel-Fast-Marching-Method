#!/usr/bin/env -S bash -e
# GMT modern mode bash template
# Date:    2020-11-12T21:45:35
# User:    xiajy
# Purpose: Purpose of this script
export GMT_SESSION_NAME=$$	# Set a unique session name
gmt begin figurename
#-R140/148.999990703/-298.5/1.5
bds=`cat boundew.gmt`

gmt figure new6ew jpg

gmt xyz2grd grid2dvew.z -Ggrid2dvew.grd -R138/150.63635055520/-298.5/1.5 -I1252+/169+ -ZLB 
gmt grdimage grid2dvew.grd -R138/150.63635055520/-298.5/1.5 -JX8.0i/2.6i -B -Cvelgrado.cpt


gmt psscale -Cvelgradoms.cpt -Ba100f100 -D10.2/9.1/12.0/0.6h 
gmt end show


