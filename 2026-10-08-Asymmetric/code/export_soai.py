# 表示用の座標（観察者側の表示。浮動小数）と、整数の判定値をまとめて書き出す
import json, numpy as np
d = json.load(open("soai_mol.json"))
P = np.array(d["P"], float) / 100.0          # Å
arms = d["arms"]; star = d["star"]
M = P * np.array([-1, 1, 1])                 # 鏡像（S体）
# 長い向きを画面の横に：主軸で回す
c = P.mean(0); u, s, vt = np.linalg.svd(P - c); Rpca = vt
if np.linalg.det(Rpca) < 0: Rpca[2] *= -1
def view(X): return (X @ Rpca.T)
def rigid_from(src, dst):                    # src の3点を dst の3点へ（回転＋平行移動、Kabsch）
    sc, dc = src.mean(0), dst.mean(0); H = (src - sc).T @ (dst - dc)
    U, S, Vt = np.linalg.svd(H); Dg = np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))])
    Rm = Vt.T @ Dg @ U.T; return Rm, dc - Rm @ sc
def axis_angle(Rm):
    ang = np.arccos(np.clip((np.trace(Rm) - 1) / 2, -1, 1))
    ax = np.array([Rm[2,1]-Rm[1,2], Rm[0,2]-Rm[2,0], Rm[1,0]-Rm[0,1]]); n = np.linalg.norm(ax)
    return (ax / n if n > 1e-9 else np.array([0,0,1.0])), ang
idx = lambda ks: [star if k == "C*" else arms[k] for k in ks]
# 場面2：中心・iPr・環 の3点を合わせる（OH と H が入れ替わる）
k2 = idx(["C*", "iPr", "環"])
R2, t2 = rigid_from(M[k2], P[k2])
ax2, an2 = axis_angle(R2)
S2end = M @ R2.T + t2
# 場面3：受け口 OH・iPr・環。板を水平にし、R体の H を上に
k3 = idx(["OH", "iPr", "環"])
R3, t3 = rigid_from(M[k3], P[k3]); S3 = M @ R3.T + t3       # S体を3つの受け口に合わせた姿
tri = P[k3]; n = np.cross(tri[1]-tri[0], tri[2]-tri[0]); n /= np.linalg.norm(n)
if np.dot(P[arms["H"]] - tri[0], n) < 0: n = -n
# n を +y へ回す
y = np.array([0,1.0,0]); v = np.cross(n, y); cth = np.dot(n, y)
K = np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
Ry = np.eye(3) + K + K @ K / (1 + cth)
f3 = lambda X: (X - tri.mean(0)) @ Ry.T
R3d, S3d = f3(P), f3(S3)
print("場面3 H の高さ（板=0）R体 %.2f  S体 %.2f Å" % (R3d[arms["H"]][1], S3d[arms["H"]][1]))
print("場面2 で残る食い違い：OH %.2f Å, H %.2f Å" % (np.linalg.norm(S2end[arms["OH"]]-P[arms["OH"]]),
      np.linalg.norm(S2end[arms["H"]]-P[arms["H"]])))
# 場面2・1 は主軸の向きで見せる
out = {"el": d["el"], "bonds": d["bonds"], "star": star, "arms": arms, "Pint": d["P"],
       "R": view(P).round(3).tolist(), "S": view(M).round(3).tolist(),
       "S2axis": (Rpca @ ax2).round(5).tolist(), "S2ang": round(float(an2), 5),
       "S2t": (Rpca @ t2).round(4).tolist(),
       "R3": R3d.round(3).tolist(), "S3": S3d.round(3).tolist(), "plate": f3(tri).round(3).tolist(), "tri": d["tri"]}
json.dump(out, open("soai_view.json", "w"), ensure_ascii=False, separators=(",", ":"))
