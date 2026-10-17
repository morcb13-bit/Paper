from common import page
model=r'''
const N8=[];for(const a of[-6,6])for(const b of[-6,6])for(const c of[-6,6])N8.push([a,b,c]);
const isAdir=d=>d.filter(x=>x>0).length%2===1;
const NG=3, C0=[6,6,6], G=grow(C0,()=>N8,NG);
const BIGA=tet([0,0,0],12,true), BIGB=tet([0,0,0],12,false);
const FACE=[[6,6,0],[6,6,12],[6,0,6],[6,12,6],[0,6,6],[12,6,6]];
// 偶数世代の原子を中心とする辺 12 の立方体（角 8 つがそろったもの）
const CUBES=G.AT.filter(a=>a.g%2===0).map(a=>{const cs=N8.map(d=>G.GEN.get(key(add(a.p,d))));return cs.every(x=>x!==undefined)?{o:a.p.map(x=>x-6),g:Math.max(...cs)}:null;}).filter(x=>x);
const CC=[...Array(NG+1).keys()].map(k=>CUBES.filter(c=>c.g<=k).length);
'''
draw=r'''
const NST=NG+2, DUR=[...Array(NST).keys()].map(i=>i===0?8:i===NST-1?9:5);
function draw(){
  g.clearRect(0,0,cv.width,cv.height);
  const R=cv.height/50, items=[], tops=[], gN=stage-1;
  const ts=stage===0?cv.height/24:cv.height/(24+15*Math.max(1,gN));
  for(let i=0;i<3;i++)CEN[i]+=(C0[i]-CEN[i])*0.05; SCALE+=(ts-SCALE)*0.06;
  if(stage===0){
    cubeLines(items,[0,0,0],12,0.6,1.4);
    const kA=ease((tS-0.5)/1.3),kB=ease((tS-2)/1.3),kC=ease((tS-3.5)/1.5);
    if(kA>0)tetDraw(items,BIGA,css("--ta"),true,kA,2,0.12);
    if(kB>0)tetDraw(items,BIGB,css("--tb"),false,kB,2,0.12);
    N8.forEach(d=>{if(kC>0)lineP(items,C0,add(C0,d.map(x=>x*kC)),css(isAdir(d)?"--ta":"--tb"),R*0.3,false,-0.3);});
    BIGA.concat(BIGB).forEach(p=>atomP(tops,p,css("--at"),R*0.55));
    if(kC>0)atomP(tops,C0,css("--at"),R*0.7*kC);
    FACE.forEach(p=>{const q=proj(p);tops.push({z:q[2]-1,k:"atom",P:[q],c:rgba(css("--muted"),0.4),r:R*0.3*q[3]});});
  } else {
    const kN=ease((tS-0.6)/1.6);
    CUBES.forEach(c=>{if(c.g>gN)return;const a=c.g===gN?kN:1;cubeLines(items,c.o,12,(c.g===gN?0.9:0.45)*a,c.g===gN?1.6:1);});
    G.BD.forEach(b=>{if(b.g>gN)return;const a=b.g===gN?ease((tS-0.1)/1.2):1;const d=[b.b[0]-b.a[0],b.b[1]-b.a[1],b.b[2]-b.a[2]];
      lineP(items,b.a,b.b.map((x,i)=>b.a[i]+(x-b.a[i])*a),rgba(css(isAdir(d)!==(b.g%2===0)?"--ta":"--tb"),b.g===gN?0.95:0.5),b.g===gN?R*0.24:R*0.14,false,-0.1);});
    G.AT.forEach(at=>{if(at.g>gN)return;const a=at.g===gN?ease((tS-1)/0.8):1;if(a>0)atomP(tops,at.p,css("--at"),R*(at.g===gN?0.55:0.42)*a);});
  }
  render(items,tops);
  texts(R,stage===0?"単位格子：角 8 と真ん中の穴に原子":`世代 ${gN}：原子 ${G.SH[gN]} 個`,
    stage===0?"角 8＝大きな A と B の頂点。真ん中から A 4 本・B 4 本。小さな点は空いた面の中心":gN===0?"真ん中の原子 1 個から、小立方体の対角線の向き 8 本で育てる":"角 8 がそろった立方体（辺 12）が、すき間なく並ぶ",
    stage===0?"結合の長さ² 108（ダイヤモンドの 27 の 4 倍）":`原子 ${G.SH.slice(0,gN+1).join("・")}`, stage===0?"原子 2（格子あたり）":`辺 12 の立方体 累計 ${CC[gN]}`,
    [["--at","原子","atom"],["--ta","A の向きの結合（大きな A）","line"],["--tb","B の向きの結合（大きな B）","dash"],["--cube","辺 12 の立方体","line"]],
    ["格子",...[...Array(NG+1).keys()].map(i=>"世代 "+i)]);
}
'''
ps=["教科書の体心立方（Na・K・Fe）は、立方体の角と真ん中に原子がある。辺 12 の格子で言えば、角 8 個は大きな A の正四面体の頂点 4 個と鏡像 B の頂点 4 個で、真ん中 (6,6,6) はシェルピンスキーの穴（正八面体）の中心にあたる。面心立方で原子のあった面の中心（シェルピンスキーの共有点）は、ここでは空いている。",
"真ん中の原子から角へ出る 8 本は、小立方体の対角線で、A の向き 4 本と B の向き 4 本に分かれる。長さ² は 108 で、ダイヤモンドの結合 27 の 4 倍（長さで 2 倍）。ダイヤモンドの炭素は A か B のどちらか 4 本だけを持ったが、体心立方の原子は両方の 8 本を持つ。",
"真ん中の原子 1 個から育てると、原子は世代ごとに 1・8・26・56 個（その先は 98・152）。偶数の世代の原子を真ん中とする辺 12 の立方体が、すき間なく並んでいく。単位格子あたりの原子は 2。角と真ん中に別のイオンを置くと、教科書の塩化セシウム（配位数 8）になる。"]
open('/home/claude/v288/bond29_bcc.html','w',encoding='utf-8').write(page("体心立方の金属","立方体の角 8 個（大きな A と鏡像 B の頂点）と真ん中の穴に原子を置き、真ん中から A と B の向き 8 本で体心立方を育てる図。ドラッグで回せる","体心立方の金属 ― 穴の中心から A と B の 8 本",ps,model,draw,3))
