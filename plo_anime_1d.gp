#-------------------------------------------------------------------------------
# gnuplotの設定
#-------------------------------------------------------------------------------
reset

set xr[0:300]
set yr[-2:2]
set term gif animate     # 出力をgifアニメに設定
set output "animation_1D.gif" # 出力ファイル名の設定

#-------------------------------------------------------------------------------
# 変数の設定
#-------------------------------------------------------------------------------
n0 = 20    # ループ変数の初期値
n1 = 1000  # ループ変数の最大値
dn = 20    # ループ変数の増加間隔

#-------------------------------------------------------------------------------
# ループの開始
#-------------------------------------------------------------------------------
load "plo_anime_1d.plt" 