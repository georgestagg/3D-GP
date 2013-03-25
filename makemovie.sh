#!/usr/bin/env bash

echo "running gnuplot

"
for count in {0080..3000};
do
if [ "$2" = "density" ]; then
echo "
set pm3d map
set title ''
set terminal png enhanced size 2048,512 crop;
set cbrange [0:$3]
set output '$count.png';
splot '$1plot.$count' using 1:2:3" | gnuplot;
fi

if [ "$2" = "phase" ]; then
echo "
set pm3d map
set title ''
set terminal png enhanced size 2048,600  crop;
set cbrange [-5:5];
set output '$count.png';
splot '$1plot.$count' using 1:2:4" | gnuplot;
fi
if [ "$2" = "vel" ]; then
echo "
set pm3d map
set title ''
set terminal png enhanced size 2048,600   crop;
set cbrange [0:5];
set output '$count.png';
splot '$1plot.$count' using 1:2:5" | gnuplot;
fi
if [ "$2" = "vort" ]; then
echo "
set title '';
unset key;
set terminal png enhanced size 720,320 crop font arial 20;
set output '44_$count.png';
set xrange [220:400];
set yrange [-40:40];
set object 1 ellipse center 980*0.4,0 size 3.0*2*2,3.0*2*2*3 angle 0. front fillstyle solid 1.0 fc rgb 'black'
set style line 9 linecolor rgb 'blue' pt 6 ps 1
set style line 7 linecolor rgb 'red' pt 7 ps 1
plot '$1plotvort.$count' u (\$3==0?-1000:\$1*0.4):(\$3==0?-1000:\$2*0.4) ps 1 pt 7 lc 3,'$1plotvort.$count' u (\$3==1?-1000:\$1*0.4):(\$3==1?-1000:\$2*0.4) ps 1 pt 9 lc 1" | gnuplot;
fi
if [ "$2" = "velx" ]; then
echo "
set pm3d map
set title ''
set terminal png enhanced size 2048,600  crop;
set cbrange [-2:2];
set output '$count.png';
splot '$1plot.$count' using 1:2:6" | gnuplot;
fi
if [ "$2" = "vely" ]; then
echo "
set pm3d map
set title ''
set terminal png enhanced size 2048,600  crop;
set cbrange [-2:2];
set output '$count.png';
splot '$1plot.$count' using 1:2:7" | gnuplot;
fi
if [ "$2" = "mach" ]; then
echo "
set pm3d map
set title ''
set terminal png enhanced size 2048,600  crop;
set cbrange [0:1.5];
set output '$count.png';
splot '$1plot.$count' using 1:2:(\$5/\$8)" | gnuplot;
fi
done


echo "converting png images to jpg
"
for i in `ls -1| grep png`;
do
convert $i $i.jpg;
done

echo "creating movie $2.avi
"
if [ "$2" = "density" ]; then
mencoder "mf://*.jpg" -mf fps=3 -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=8000 -o density.avi 
fi
if [ "$2" = "phase" ]; then
mencoder "mf://*.jpg" -mf fps=30 -ovc xvid -ovc lavc  -o phase.avi
fi
if [ "$2" = "vort" ]; then
mencoder "mf://*.jpg" -mf fps=30 -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=8000 -o vort.avi 
fi
if [ "$2" = "vel" ]; then
mencoder "mf://*.jpg" -mf fps=30 -ovc xvid -ovc lavc -o vel.avi
fi
if [ "$2" = "velx" ]; then
mencoder "mf://*.jpg" -mf fps=30 -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=800 -o velx.avi
fi
if [ "$2" = "vely" ]; then
mencoder "mf://*.jpg" -mf fps=30 -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=800 -o vely.avi
fi
if [ "$2" = "mach" ]; then
mencoder "mf://*.jpg" -mf fps=30 -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=800 -o mach.avi
fi
rm -Rf *.jpg *.png

