set term jpeg

#表示するx, y, zの範囲
set xrange[0:1]
set yrange[0:1]

#軸ラベルの設定
set xlabel'x'
set ylabel'y'

#縦横比の設定
set size ratio -1

#出力するファイルを設定
set output 'vector.jpg'

set palette rgb 33,13,10
plot 'Vector\vec_01000.dat' u 1:2:3:4:(sqrt($3*$3+$4*$4)) w vector lc palette ti ""
