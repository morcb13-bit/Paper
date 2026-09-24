# pc5 一式（2026-09-24）
置き場：2026-09-21-plant-sun/code/（b13_chain_units.py・gen10.py と同じフォルダで走る）
中間の pickle は /home/claude/prop/ に書く。走らせる順：
  slit14.py → br2.py（鏡の対の一覧）→ g10b.py（10枚の担体の隣接）→ pc5.py
各検定：if1b（IF1残り2行）ha1（HA1）rg1（RG1）ps1（PS1）add2（AD2）ad3（AD3）cs1（CS1）loop1（LP1）ly1（LY1）pc5（D1〜D3）
ac1.py  AC1 アキュムレータ（pc5.py の後に走らせる）
m1.py   M1 命令列（ADD/NEG/SHL/SHR/JNZ/JNEG/JMP）で掛け算・割り算（ac1.py の後）
sp1.py  SP1 プログラム内蔵（メモリ30語・累算器の機械、m1.py の後）
fl1.py  FL1 フィボナッチとリュカを内蔵プログラムで計算し担体の φ^n と照合（sp1.py の後）
sv1.py  SV1 素数（絶対番地）と幸運数（相対番地）の篩（sp1.py の後、fast.py を読む。幸運数は数分かかる）
vp1.py  VP1 型紙だけで描く一桁の加算器と、桁を描き足す加算器（template_digit.json を書き出す）
nr1-3.py NR ひし形と同じ中心の相似な黄金のひし形（比φ³）／ vp2.py VP2 ひし形の入れ子で次の桁の型紙を置く
alt1.py ALT1 二つの仮想ペンローズ（五角形の層Aと円環の中心の層B）を交互に描いて道を決める
alt3.py/alt3b.py ALT3 三階建て（五芒星の層の噛み合いで ±φ³ の向きを決める）
