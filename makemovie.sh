#!/usr/bin/env bash
echo "running gnuplot
"
for count in {0001..0300};
do
if [ "$2" = "density" ]; then
echo "
set pm3d map
set title ''
set terminal png enhanced size 512,512 crop;
set cbrange [0:$3]
set output '$count.png';
splot '$1plot.$count' using 1:2:3" | gnuplot;
fi

if [ "$2" = "phase" ]; then
echo "
set pm3d map
set title ''
set terminal png enhanced size 512,512  crop;
set cbrange [-1.57:1.57];
set output '$count.png';
splot '$1plot.$count' using 1:2:4" | gnuplot;
fi
if [ "$2" = "vort" ]; then
echo "
set title '';
unset key;
set terminal png enhanced size 1024,1024 crop;
set output '44_$count.png';
set xrange [-60:60];
set yrange [-60:60];
set object 1 ellipse center 0,0 size 1.5,2*1.5 angle 0. front fillstyle solid 1.0 fc rgb 'black'
set style line 9 linecolor rgb 'blue' pt 6 ps 1
set style line 7 linecolor rgb 'red' pt 7 ps 1
plot '$1plotvort.$count' u (\$3==0?-1000:\$1):(\$3==0?-1000:\$2) ps 1 pt 7 lc 3,'$1plotvort.$count' u (\$3==1?-1000:\$1):(\$3==1?-1000:\$2) ps 1 pt 9 lc 1" | gnuplot;
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
rm -Rf *.jpg *.png

