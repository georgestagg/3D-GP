#!/usr/bin/env bash

echo "running gnuplot

"
for count in {0001..3000};
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
set terminal png enhanced size 2048,512 crop;
set cbrange [-5:5];
set output '$count.png';
splot '$1plot.$count' using 1:2:4" | gnuplot;
fi
if [ "$2" = "vel" ]; then
echo "
set pm3d map
set title ''
set terminal png enhanced size 2048,512  crop;
set cbrange [0:5];
set output '$count.png';
splot '$1plot.$count' using 1:2:5" | gnuplot;
fi
if [ "$2" = "vorticity" ]; then
echo "
set title '';
unset key;
set terminal png enhanced size 2000,480 crop;
set output '$count.png';
set xrange [-250:250];
set yrange [-60:60];
plot '$1plotvort.$count' u (\$3==0?-1000:\$1):(\$3==0?-1000:\$2) ps 1 pt 7,'$1plotvort.$count' u (\$3==1?-1000:\$1):(\$3==1?-1000:\$2) ps 1 pt 9" | gnuplot;
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
mencoder "mf://*.jpg" -mf fps=30 -ovc xvid -ovc lavc -lavcopts vcodec=ffv1 -o density.avi 
fi
if [ "$2" = "phase" ]; then
mencoder "mf://*.jpg" -mf fps=30 -ovc xvid -ovc lavc -lavcopts vcodec=ffv1 -o phase.avi
fi
if [ "$2" = "vorticity" ]; then
mencoder "mf://*.jpg" -mf fps=30 -ovc xvid -ovc lavc -lavcopts vcodec=ffv1 -o vorticity.avi
fi
if [ "$2" = "vel" ]; then
mencoder "mf://*.jpg" -mf fps=30 -ovc xvid -ovc lavc -lavcopts vcodec=ffv1 -o vel.avi
fi

rm -Rf *.jpg *.png

