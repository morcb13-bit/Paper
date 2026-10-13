# 検定ION：席の配置が同じで核の電荷だけ違う組。予言（事前）：陽イオンは中性より開く、陰イオンは閉じる
# 判定は角の大小（記録の値どうしの比較）だけ。None は記録なし
REC={"NH3":106.769,"H3O+":111.95,"CH3-":110.0,"H2O":104.4776,"H2F+":112.30,"NH2-":None,
     "PH3":93.4,"H3S+":None,"SiH3-":None,"H2S":92.11,"H2Cl+":94.24,"PH2-":None}
NOTE={"H3O+":"予言を作る前に見ていた（当てはめ）","CH3-":"1977年SCFの約110°、弱い記録"}
for neu,cat,ani in [("NH3","H3O+","CH3-"),("H2O","H2F+","NH2-"),("PH3","H3S+","SiH3-"),("H2S","H2Cl+","PH2-")]:
    n=REC[neu]
    for ion,want in ((cat,"開く"),(ani,"閉じる")):
        v=REC[ion]
        if v is None: res="記録なし"
        else:
            ok=(v>n) if want=="開く" else (v<n)
            res=("OK" if ok else "NG")+f"（{v} 対 {n}）"
        print(f"{neu:4} → {ion:6} 予言:{want:3}  {res}  {NOTE.get(ion,'')}")
