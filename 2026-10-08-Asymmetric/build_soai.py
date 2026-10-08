# 硤合反応の生成物 1-[2-(tert-ブチルエチニル)ピリミジン-5-イル]-2-メチルプロパン-1-オール の三次元模型
# 座標は RDKit の ETKDG＋MMFF（標準的な結合長・角から組んだもの。結晶の実測ではない）
import json, math
from rdkit import Chem
from rdkit.Chem import AllChem
smi = "CC(C)[C@@H](O)c1cnc(nc1)C#CC(C)(C)C"
m = Chem.AddHs(Chem.MolFromSmiles(smi))
AllChem.EmbedMolecule(m, randomSeed=13); AllChem.MMFFOptimizeMolecule(m)
cip = Chem.FindMolChiralCenters(m)[0]
star = cip[0]
conf = m.GetConformer()
# 0.01Å 単位の整数に直す（中心炭素を原点）
c0 = conf.GetAtomPosition(star)
P = [[round((conf.GetAtomPosition(i).x - c0.x) * 100), round((conf.GetAtomPosition(i).y - c0.y) * 100),
      round((conf.GetAtomPosition(i).z - c0.z) * 100)] for i in range(m.GetNumAtoms())]
el = [a.GetSymbol() for a in m.GetAtoms()]
bonds = [[b.GetBeginAtomIdx(), b.GetEndAtomIdx(), int(b.GetBondTypeAsDouble() * 2)] for b in m.GetBonds()]
arms = [n.GetIdx() for n in m.GetAtomWithIdx(star).GetNeighbors()]
def kind(i):
    a = m.GetAtomWithIdx(i)
    if a.GetSymbol() == "H": return "H"
    if a.GetSymbol() == "O": return "OH"
    return "環" if a.GetIsAromatic() else "iPr"
armk = {kind(i): i for i in arms}
def sub(a, b): return [a[k] - b[k] for k in range(3)]
def det3(u, v, w): return u[0]*(v[1]*w[2]-v[2]*w[1]) - u[1]*(v[0]*w[2]-v[2]*w[0]) + u[2]*(v[0]*w[1]-v[1]*w[0])
A, B, C, D = (P[armk[k]] for k in ("H", "OH", "iPr", "環"))
s_R = det3(sub(B, A), sub(C, A), sub(D, A))
mir = [[-x, y, z] for x, y, z in P]
A2, B2, C2, D2 = (mir[armk[k]] for k in ("H", "OH", "iPr", "環"))
s_L = det3(sub(B2, A2), sub(C2, A2), sub(D2, A2))
L2 = {k: sum(c*c for c in P[armk[k]]) for k in armk}
print("CIP", cip, "腕の長さの二乗(0.01Å単位)", L2)
print("s(この分子) =", s_R, " s(鏡像) =", s_L)
# 型紙：受け口 = OH, iPr, 環 の3つ。板の表側に H が来るか = det[C-B, D-B, A-B] の符号
side = lambda Q: det3(sub(Q[armk["iPr"]], Q[armk["OH"]]), sub(Q[armk["環"]], Q[armk["OH"]]), sub(Q[armk["H"]], Q[armk["OH"]]))
tri = lambda Q: [sum(c*c for c in sub(Q[armk[x]], Q[armk[y]])) for x, y in (("OH","iPr"),("iPr","環"),("環","OH"))]
print("受け口の三角形（辺の二乗）", tri(P), tri(mir), " H の側", side(P), side(mir))
json.dump({"el": el, "P": P, "bonds": bonds, "star": star, "arms": armk, "cip": cip[1],
           "sR": s_R, "sL": s_L, "sideR": side(P), "sideL": side(mir), "tri": tri(P)},
          open("soai_mol.json", "w"), ensure_ascii=False, separators=(",", ":"))
print(Chem.MolToSmiles(Chem.RemoveHs(m)), len(el), "原子")
