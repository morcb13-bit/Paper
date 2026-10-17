from common import page
model=r'''
const FCC=p=>p.every(x=>x%6===0)&&((((p[0]+p[1]+p[2])/6)%2)+2)%2===0;
const N12=[];for(const a of[-6,0,6])for(const b of[-6,0,6])for(const c of[-6,0,6]){const d=[a,b,c];if(dot(d,d)===72)N12.push(d);}
const NG=3, G=grow([0,0,0],()=>N12,NG);
const TETS=[];for(let i=-4;i<=3;i++)for(let j=-4;j<=3;j++)for(let k=-4;k<=3;k++){const n=[i,j,k],o=[6*i,6*j,6*k];
  const T=BITS.map(b=>add(o,[6*b[0],6*b[1],6*b[2]])).filter(FCC);const gs=T.map(v=>G.GEN.get(key(v)));
  if(gs.every(x=>x!==undefined))TETS.push({T,A:((i+j+k)%2+2)%2===1,g:Math.max(...gs)});}
const TC=[...Array(NG+1).keys()].map(k=>TETS.filter(t=>t.g<=k).length);
const CELLT=BITS.map(n=>({T:BITS.map(b=>[6*n[0]+6*b[0],6*n[1]+6*b[1],6*n[2]+6*b[2]]).filter(FCC),A:odd(n)}));
const CELLA=[];for(const a of[0,6,12])for(const b of[0,6,12])for(const c of[0,6,12])if(FCC([a,b,c]))CELLA.push([a,b,c]);
'''
draw=r'''
const NST=NG+2, DUR=[...Array(NST).keys()].map(i=>i===0?8:i===NST-1?9:5);
function draw(){
  g.clearRect(0,0,cv.width,cv.height);
  const R=cv.height/50, items=[], tops=[], gN=stage-1;
  const tc=stage===0?[6,6,6]:[0,0,0], ts=stage===0?cv.height/24:cv.height/(20+11*Math.max(1,gN));
  for(let i=0;i<3;i++)CEN[i]+=(tc[i]-CEN[i])*0.05; SCALE+=(ts-SCALE)*0.06;
  if(stage===0){
    cubeLines(items,[0,0,0],12,0.55,1.3);BITS.forEach(n=>cubeLines(items,[6*n[0],6*n[1],6*n[2]],6,0.2,1));
    CELLT.forEach((t,i)=>{const k=ease((tS-0.5-(t.A?0:2))/1.4);if(k>0)tetDraw(items,t.T,css(t.A?"--ta":"--tb"),t.A,k,2,0.15,0.3+0.7*k);});
    CELLA.forEach(p=>atomP(tops,p,css("--at"),R*0.55));
    const k=ease((tS-4.5)/1);if(k>0){[[6,6,6],[6,0,0],[0,6,0],[0,0,6],[6,12,12],[12,6,12],[12,12,6],[12,6,0],[6,12,0],[0,6,12],[0,12,6],[12,0,6],[6,0,12]].forEach(p=>{const q=proj(p);tops.push({z:q[2]-1,k:"atom",P:[q],c:rgba(css("--muted"),0.35*k),r:R*0.32*q[3]});});}
  } else {
    const kN=ease((tS-0.6)/1.6);
    TETS.forEach(t=>{if(t.g>gN)return;const a=t.g===gN?kN:1;tetDraw(items,t.T,css(t.A?"--ta":"--tb"),t.A,(t.g===gN?1:0.55)*a,t.g===gN?1.8:1,0.12,0.3+0.7*a);});
    G.BD.forEach(b=>{if(b.g>gN)return;const a=b.g===gN?ease((tS-0.1)/1.2):1;lineP(items,b.a,b.b.map((x,i)=>b.a[i]+(x-b.a[i])*a),rgba(css("--share"),b.g===gN?0.8:0.3),b.g===gN?R*0.18:R*0.1,false,-0.1);});
    G.AT.forEach(at=>{if(at.g>gN)return;const a=at.g===gN?ease((tS-1)/0.8):1;if(a>0)atomP(tops,at.p,css("--at"),R*(at.g===gN?0.5:0.4)*a);});
  }
  render(items,tops);
  texts(R,stage===0?"単位格子：シェルピンスキーの頂点だけに原子":`世代 ${gN}：原子 ${G.SH[gN]} 個`,
    stage===0?"小立方体 8 個に正四面体が 1 個ずつ（奇に A、偶に B）。小さな点は空いた正八面体の中心":gN===0?"原子 1 個から、正四面体の辺の向き 12 本で育てる":"新しい正四面体は、隣の正四面体と辺を 1 本分け合う",
    stage===0?"体積 1728＝正四面体 8×72＋正八面体 4×288":`原子 ${G.SH.slice(0,gN+1).join("・")}`, stage===0?"原子 4（格子あたり）":`正四面体 累計 ${TC[gN]}`,
    [["--at","原子","atom"],["--ta","A の正四面体","line"],["--tb","鏡像 B の正四面体","dash"],["--share","金属結合（正四面体の辺）","thick"]],
    ["格子",...[...Array(NG+1).keys()].map(i=>"世代 "+i)]);
}
'''
ps=["教科書の面心立方（Cu・Al・Ag・Au）は、立方体の角と各面の中心に原子がある。辺 12 の格子で言えば、番地の和が偶数の格子点すべてで、シェルピンスキーの 2 段目の頂点 10 個と、残りの角 4 個がそのまま原子の位置になる。",
"格子を 2×2×2 に割ると、どの小立方体にも原子 4 個が作る正四面体が 1 個ずつ入る。番地の和が奇数の小立方体では A、偶数では鏡像 B。A と B のシェルピンスキーの 2 段目を重ねた元素の枠そのもので、隣どうしの正四面体は辺を 1 本分け合う。真ん中と辺の中点は正八面体の穴で、何も入らない。格子の体積 1728 は、正四面体 8 個（72 ずつ）と正八面体 4 個（288 ずつ）でちょうど埋まる。",
"原子 1 個から、正四面体の辺の向き 12 本（長さ² 72）で育てると、原子は世代ごとに 1・12・42・92 個。原子 1 個は正四面体 8 個と正八面体 6 個の頂点になる。塩化ナトリウムの Na だけを残し、正八面体の穴の Cl を抜いた形でもある。"]
open('/home/claude/v288/bond28_fcc.html','w',encoding='utf-8').write(page("面心立方の金属","教科書の面心立方の格子をシェルピンスキーの A と B の 2 段目として描き、原子 1 個から正四面体の辺の向きで育てる図。ドラッグで回せる","面心立方の金属 ― A と B のシェルピンスキーを重ねた枠",ps,model,draw,2))
