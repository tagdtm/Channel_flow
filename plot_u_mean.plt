set term jpeg

#表示するx, y, zの範囲
set xrange[0.1:150]
set yrange[0:20]

#軸ラベルの設定
set xlabel'y+'
set ylabel'U mean'

set logscale x

#縦横比の設定
set size ratio 20/150

#出力するファイルを設定
set output 'u_mean.jpg'

set key left top

set palette rgb 33,13,10
plot \
  (x<15) ? x:1/0 w line lc 1 title 'theoretical line',\
  (x>10) ? log(x)/0.41+5.5:1/0 w line lc 1 title '',\
  'model data/velocity_mean.dat' using 3:4 w points pt 1 lc 2 title 'tht data', \
  'u_sta.dat' using ($1*300):2 w points pt 5 lc 3 title 'calculation result'
