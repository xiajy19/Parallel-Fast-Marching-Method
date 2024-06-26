#!/usr/bin/env -S bash -e
# GMT modern mode bash template
# Date:    2020-11-13T10:28:34
# User:    xiajy
# Purpose: Purpose of this script
export GMT_SESSION_NAME=$$	# Set a unique session name
gmt begin figurename

bds=`cat boundns.gmt`
 
#
#####################################
# Plot velocity field
#####################################
gmt figure new6ns jpg
gmt xyz2grd grid2dvns.z -Ggrid2dvns.grd -R-39.08858000000/-31.5/-298.5/1.5 -I655+/169+ -ZRB 
gmt grdimage grid2dvns.grd -R-39.08858000000/-31.5/-298.5/1.5 -JX8.0i/2.60i -B -Cvelgrado.cpt

gmt psscale -Cvelgradoms.cpt -B -D10.2/9.1/10.0/0.6h 

gmt text -F+f15p+jML -N << EOF
-33.3 98.7 @~\144@~Vp(m/s)
EOF
gmt end show
