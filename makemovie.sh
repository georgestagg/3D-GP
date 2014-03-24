#!/bin/bash
# Compile and set up the code ready to run in a separate directory
# usage:  ./makemovie.sh <velocity> <type> <cbrange> <xrange> <yrange>

nx=`grep -F "NX" < params.in | grep -o "[0-9]*"`
ny=`grep -F "NY" < params.in | grep -o "[0-9]*"`
nofiles=`ls $1plot* |tac| grep -o -m 1 "\.[0-9]*" | grep -o "[0-9]*"`
scaledy=`echo "1024*($ny/$nx)" | bc -l | awk '{printf("%d\n",$1 + 0.5)}'`
scaledyvort=`echo "1577*($ny/$nx)*0.927" | bc -l | awk '{printf("%d\n",$1 + 0.5)}'`
dspace=`grep -F "DSPACE" < params.in | grep -o "=.*$" | grep -o "[0-9.-]*" | xargs |awk '{printf("%f\n",$1*(10**$2))}'`

if [ -z "$4" ]; then
	xrange=`echo "$nx*$dspace/2" | bc -l`
else
	xrange=$4
fi

if [ -z "$5" ]; then
	yrange=`echo "$ny*$dspace/2" | bc -l`
else
	yrange=$5
fi

echo "running gnuplot - nx = $nx - ny = $ny - dspace = $dspace - xrange = $xrange - yrange = $yrange"
for count in $(eval echo "{0001..$nofiles}")
do
if [ "$2" = "density" ]; then
echo "
set pm3d map
set title ''
set terminal png enhanced size 1024,$scaledy crop;
set cbrange [0:$3]
set xrange[-$xrange:$xrange]
set yrange[-$yrange:$yrange]
set output '$count.png';
splot '$1plotwf.$count' using 1:2:(sqrt(\$3**2+\$4**2))" | gnuplot;
fi

if [ "$2" = "phase" ]; then
echo "
set pm3d map
set title ''
set terminal png enhanced size 2048,$scaledy  crop;
set xrange[-$xrange:$xrange]
set yrange[-$yrange:$yrange]
set cbrange [-1.57:1.57];
set output '$count.png';
splot '$1plot.$count' using 1:2:4" | gnuplot;
fi
if [ "$2" = "vort" ]; then
echo "
set title '';
unset key;
set terminal png enhanced size 1577,$scaledyvort crop;
set xrange[-$xrange:$xrange]
set yrange[-$yrange:$yrange]
unset colorbox
set output '$count.v.png';
set style line 9 linecolor rgb 'blue' pt 6 ps 2
set style line 7 linecolor rgb 'green' pt 7 ps 2
plot '$1plotvort.$count' u (\$3==0?-1000:\$1):(\$3==0?-1000:\$2) ps 2 pt 7 lc 3,'$1plotvort.$count' u (\$3==1?-1000:\$1):(\$3==1?-1000:\$2) ps 2 pt 9 lc 1" | gnuplot;
fi
if [ "$2" = "densvort" ]; then
echo "
set pm3d map
set title ''
set terminal png enhanced size 2048,$scaledy crop;
set xrange[-$xrange:$xrange]
set yrange[-$yrange:$yrange]
set cbrange [0:$3]
unset colorbox
set output '$count.d.png';
splot '$1plot.$count' using 1:2:3" | gnuplot;
echo "
set title '';
unset key;
set terminal png enhanced size 1577,$scaledyvort crop;
set xrange[-$xrange:$xrange]
set yrange[-$yrange:$yrange]
unset colorbox
set output '$count.v.png';
set ytics nomirror;
unset xtics;
set border 7;
set style line 9 linecolor rgb 'blue' pt 6 ps 2
set style line 7 linecolor rgb 'green' pt 7 ps 2
plot '$1plotvort.$count' u (\$3==0?-1000:\$1):(\$3==0?-1000:\$2) ps 2 pt 7 lc 3,'$1plotvort.$count' u (\$3==1?-1000:\$1):(\$3==1?-1000:\$2) ps 2 pt 9 lc 1" | gnuplot;
convert $count.v.png -transparent white $count.v.png
composite $count.v.png $count.d.png $count.f.png
rm $count.v.png $count.d.png
fi
done

echo "creating movie $2.avi"
mencoder "mf://*.png" -mf fps=30 -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=20971520 -o $2.avi 

rm -Rf *.jpg *.png

