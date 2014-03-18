#!/usr/bin/env bash
echo "running gnuplot
"
for count in {0001..1500};
do
if [ "$2" = "density" ]; then
echo "
set pm3d map
set title ''
set terminal png enhanced size 2048,1024 crop;
set cbrange [0:$3]
set output '$count.png';
splot '$1plot.$count' using 1:2:3" | gnuplot;
fi

if [ "$2" = "phase" ]; then
echo "
set pm3d map
set title ''
set terminal png enhanced size 1000,1000  crop;
set xrange[-20:20]
set yrange[-20:20]
set cbrange [-1.57:1.57];
set output '$count.png';
splot '$1plot.$count' using 1:2:4" | gnuplot;
fi
if [ "$2" = "vort" ]; then
echo "
set title '';
unset key;
set terminal png enhanced size 1600,480 crop;
set xrange[0:420]
set yrange[-60:60]
set output '$count.png';
set object 1 ellipse center 393,0 size 4.2*2.4,4.2*4.8 angle 0. front fillstyle solid 1.0 fc rgb 'black'
set style line 9 linecolor rgb 'blue' pt 6 ps 1
set style line 7 linecolor rgb 'red' pt 7 ps 1
plot '$1plotvort.$count' u (\$3==0?-1000:\$1):(\$3==0?-1000:\$2) ps 1 pt 7 lc 3,'$1plotvort.$count' u (\$3==1?-1000:\$1):(\$3==1?-1000:\$2) ps 1 pt 9 lc 1" | gnuplot;
fi
if [ "$2" = "densvort" ]; then
echo "
set pm3d map
set title ''
set terminal png enhanced size 2048,1024 crop;
set xrange[-200:200]
set yrange[-120:120]
set cbrange [0:$3]
unset colorbox
set output '$count.d.png';
splot '$1plot.$count' using 1:2:3" | gnuplot;
echo "
set title '';
unset key;
set terminal png enhanced size 1577,748 crop;
set xrange[-200:200]
set yrange[-120:120]
unset colorbox
set output '$count.v.png';
set style line 9 linecolor rgb 'blue' pt 6 ps 2
set style line 7 linecolor rgb 'green' pt 7 ps 2
plot '$1plotvort.$count' u (\$3==0?-1000:\$1):(\$3==0?-1000:\$2) ps 2 pt 7 lc 3,'$1plotvort.$count' u (\$3==1?-1000:\$1):(\$3==1?-1000:\$2) ps 2 pt 9 lc 1" | gnuplot;
convert $count.v.png -transparent white $count.v.png
composite $count.v.png $count.d.png $count.f.png
fi
done

echo "creating movie $2.avi
"
if [ "$2" = "density" ]; then
mencoder "mf://*.png" -mf fps=30 -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=20971520 -o density.avi 
fi
if [ "$2" = "phase" ]; then
mencoder "mf://*.png" -mf fps=30 -ovc xvid -ovc lavc  -o phase.avi
fi
if [ "$2" = "vort" ]; then
mencoder "mf://*.png" -mf fps=30 -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=20971520 -o vort.avi 
fi
if [ "$2" = "densvort" ]; then
mencoder "mf://*.f.png" -mf fps=30 -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=20971520 -o densvort.avi 
fi
rm -Rf *.jpg *.png

