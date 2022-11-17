set term jpeg

#表示するx, y, zの範囲
set xrange[0:150]
set yrange[0:3]

#軸ラベルの設定
set xlabel'y+'
set ylabel'u-RMS'

#縦横比の設定
set size ratio 3/150

#出力するファイルを設定
set output 'u_rms.jpg'

set key left top

set palette rgb 33,13,10
plot \
  'model data/u_rms.dat' using 3:4 w points pt 7 lc 3 title 'tht data'
  'v_rms.dat' using ($1*300):2 w points pt 6 lc 4 title 'calculation result'
