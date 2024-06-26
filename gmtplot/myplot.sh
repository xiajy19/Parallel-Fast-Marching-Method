#!/usr/bin/env -S bash -e
# GMT modern mode bash template
# Date:    2020-11-12T20:37:23
# User:    xiajy
# Purpose: Purpose of this script
export GMT_SESSION_NAME=$$	# Set a unique session name
gmt begin figurename


# -R138/150.78/-39.626/-31.5

#bounds = "-R140/149/-38/-32.5"
#proj = "-Jl144:23/-38/-34/-36/2.0"
#psfile = "gmtslice.ps"
#-JM144.5/-35.25/20

set ANOT_FONT Helvetica ANOT_FONT_SIZE 16p

gmt figure new6 jpg

gmt xyz2grd grid2dvd.z -Ggrid2dvd.grd -R138/150.63635055520/-39.08857981470/-31.5 -I835+/655+ -ZLB 

gmt grdimage grid2dvd.grd -R138/150.63635055520/-39.08857981470/-31.5 -JM145/-35/20 -Cvelgrado.cpt  

gmt psscale -Cvelgradoms.cpt -Ba100f100 -D10.5/13.2/8.0/0.4h  

gmt psxy seal.dat -R138/150.63635055520/-39.08857981470/-31.5 -JM20 -: -St0.40 -G200/50/50 -W3  

gmt coast -R138/150.63635055520/-39.08857981470/-31.5 -JM20 -Ia -C -W0.5p -A2 -B -Dh  
	
	
gmt end show
