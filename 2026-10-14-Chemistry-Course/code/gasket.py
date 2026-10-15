# シェルピンスキーのギャスケット（平面、6 段）。三角形の番地は整数、画面に置く高さの換算だけ浮動小数（観察者側の表示）
import math
from PIL import Image, ImageDraw
K=6; N=2**K                      # 一辺を 64 に割る
def tri(x,y,s,k,out):            # 左下 (x,y)、一辺 s（整数）
    if k==0: out.append((x,y,s)); return
    h=s//2
    tri(x,y,h,k-1,out); tri(x+h,y,h,k-1,out); tri(x+h//2*1,y+h,h,k-1,out) if h%2==0 else None
out=[]
def rec(x2,y,s,k):              # x2 は横位置の 2 倍（半分ずらしを整数で持つ）
    if k==0: out.append((x2,y,s)); return
    h=s//2
    rec(x2,y,h,k-1); rec(x2+2*h,y,h,k-1); rec(x2+h,y+h,h,k-1)
rec(0,0,N,K)
W=1600; pad=60; unit=(W-2*pad)/N; hh=math.sqrt(3)/2   # 観察者側の表示
H=int(N*unit*hh+2*pad)
def P(x2,y): return (pad+x2*unit/2, H-pad-y*unit*hh)
img=Image.new("RGB",(W,H),"white"); d=ImageDraw.Draw(img)
svg=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}"><rect width="100%" height="100%" fill="#ffffff"/>']
for x2,y,s in out:
    pts=[P(x2,y),P(x2+2*s,y),P(x2+s,y+s)]
    d.polygon(pts,fill=(47,111,143))
    svg.append('<polygon fill="#2f6f8f" points="'+" ".join(f"{a:.2f},{b:.2f}" for a,b in pts)+'"/>')
svg.append("</svg>")
img.save("/home/claude/v287/img/sierpinski_gasket.png")
open("/home/claude/v287/img/sierpinski_gasket.svg","w").write("\n".join(svg))
print(len(out),"triangles")
