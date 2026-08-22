#  差し替え（実体は届いていない）。sink.py の前半が要る名前 CELLS・XYC・O・math・random だけを
#  global_pull.py の前半から供給する。放物線の軌道・ring_of・follow はここには無い。
exec(open('global_pull.py').read().split('for MB in (6, 8):')[0])
print("■ 投げられたものを追う（差し替え・ここから下は無い）")
