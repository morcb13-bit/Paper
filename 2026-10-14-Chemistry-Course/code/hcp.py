from common import page
model=r'''
const E1=[6,-6,0], E2=[0,6,-6];
const inDisk=(p,lim)=>{const s=p[0]+p[1]+p[2];return 3*dot(p,p)-s*s<=lim;};   // 面内の距離² の 3 倍で比べる（整数）
function layer(base,lim){const out=[];for(let i=-8;i<=8;i++)for(let j=-8;j<=8;j++){const p=add(base,[i*E1[0]+j*E2[0],i*E1[1]+j*E2[1],i*E1[2]+j*E2[2]]);
  const c=p.map((x,t)=>x-base[t]);if(inDisk(c,lim))out.push(p);}return out;}
const LIM=3*72*4;
const L0=layer([0,0,0],LIM), L1=layer([6,6,0],LIM), L2=layer([8,8,8],LIM), LC=layer([12,12,0],LIM);
const d2=(a,b)=>{const d=[a[0]-b[0],a[1]-b[1],a[2]-b[2]];return dot(d,d);};
function tets(lo,hi){const out=[];
  for(const p of hi){const nb=lo.filter(q=>d2(p,q)===72);if(nb.length===3)out.push({T:[...nb,p],up:true,face:nb.map(key).sort().join("|")});}
  for(const p of lo){const nb=hi.filter(q=>d2(p,q)===72);if(nb.length===3)out.push({T:[...nb,p],up:false,face:nb.map(key).sort().join("|")});}
  return out;}
const T01=tets(L0,L1), T12=tets(L1,L2), T1C=tets(L1,LC);
const PAIR=new Set(T01.filter(t=>!t.up).map(t=>t.face).filter(f=>T12.some(u=>u.up&&u.face===f)));
const PAIRC=T01.filter(t=>!t.up).map(t=>t.face).filter(f=>T1C.some(u=>u.up&&u.face===f)).length;
// 層の中の三角形（L0）
const TRI0=[];for(const a of L0)for(const b of L0)for(const c of L0){if(key(a)<key(b)&&key(b)<key(c)&&d2(a,b)===72&&d2(b,c)===72&&d2(a,c)===72)TRI0.push([a,b,c]);}
const mid=L1.reduce((m,p)=>dot(p,p)<dot(m,m)?p:m,L1[0]);
const NB=[L1.filter(q=>d2(mid,q)===72).length,L0.filter(q=>d2(mid,q)===72).length,L2.filter(q=>d2(mid,q)===72).length];
'''
draw=r'''
UPAX=true; pitch=0.5;
const NST=4, DUR=[6,7,8,9];
function draw(){
  g.clearRect(0,0,cv.width,cv.height);
  const R=cv.height/50, items=[], tops=[];
  const tc=stage===0?[0,0,0]:stage===1?[3,3,3]:[5,5,5]; for(let i=0;i<3;i++)CEN[i]+=(tc[i]-CEN[i])*0.05; SCALE+=(cv.height/58-SCALE)*0.06;
  // 層 0
  TRI0.forEach(t=>{const V=t.map(proj);items.push({z:(V[0][2]+V[1][2]+V[2][2])/3+0.4,k:"face",P:V,c:rgba(css("--cube"),0.12)});
    for(let i=0;i<3;i++)lineP(items,t[i],t[(i+1)%3],rgba(css("--cube"),0.7),1.2,false,0.2);});
  L0.forEach(p=>atomP(tops,p,css("--at"),R*0.5));
  if(stage>=1){const k=stage===1?ease((tS-0.3)/1.6):1;
    T01.forEach(t=>tetDraw(items,t.T,css(t.up?"--ta":"--tb"),t.up,(stage===1?ease((tS-2)/1.5):0.45)*(PAIR.has(t.face)&&stage===2&&!t.up?2:1),1.4,0.1));
    L1.forEach(p=>{const q=add(p,[Math.round(0),0,0]);const lift=(1-k)*8;atomP(tops,add(p,[lift,lift,lift]),css("--e"),R*0.5*k);});}
  if(stage===2){const k=ease((tS-0.3)/1.6);
    T12.forEach(t=>{const pr=t.up&&PAIR.has(t.face);tetDraw(items,t.T,css(pr?"--share":(t.up?"--ta":"--tb")),t.up,ease((tS-2)/1.5)*(pr?1:0.5),pr?2.2:1.2,pr?0.2:0.08);});
    T01.forEach(t=>{if(!t.up&&PAIR.has(t.face))tetDraw(items,t.T,css("--share"),false,ease((tS-3)/1.5),2.2,0.2);});
    L2.forEach(p=>{const lift=(1-k)*8;atomP(tops,add(p,[lift,lift,lift]),css("--at"),R*0.5*k);});}
  if(stage===3){const k=ease((tS-0.3)/1.6);
    T1C.forEach(t=>tetDraw(items,t.T,css(t.up?"--ta":"--tb"),t.up,0.5*ease((tS-1.5)/1.5),1.2,0.08));
    LC.forEach(p=>{const lift=(1-k)*8;atomP(tops,add(p,[lift,lift,lift]),rgba(css("--at"),0.55),R*0.45*k);});}
  render(items,tops);
  const T0=["1 層目：三角形に並べた層（A）","2 層目：くぼみに入れる（B）","3 層目：A の真上に戻す（六方最密）","比べる：3 層目を C に置くと面心立方"][stage];
  const T1=["黒鉛と同じ三角形。黒鉛は三角形の中心に原子、ここは頂点だけ","上の層は (6,6,0)・(6,0,6)・(0,6,6) の向きに乗る。長さ² 72、立方体の枠の向き","上の層は (2,2,8)・(2,8,2)・(8,2,2) の向きに乗る。長さ² 72、枠の外の向き",`(6,6,0) の向きで乗り続けると 3 層目は C。B の層をはさんで面を分け合う組は ${PAIRC}`][stage];
  texts(R,T0,T1,stage>=2?`B 層の原子の隣：同じ層 ${NB[0]}・下 ${NB[1]}・上 ${NB[2]}＝12`:"",stage===2?`B の三角形をはさんで面を分け合う正四面体の組 ${PAIR.size}（紫）`:"",
    [["--at","A 層の原子","atom"],["--e","B 層の原子","atom"],["--ta","上向きの正四面体","line"],["--tb","下向きの正四面体","dash"],["--share","面を分け合う 2 個組","thick"]],
    ["A","B","A（六方）","C（面心）"]);
}
'''
ps=["教科書の六方最密（Mg・Zn・Ti）は、三角形に並べた層を、くぼみに次の層を入れながら A・B・A・B と 2 枚ごとに戻るように積む。3 枚ごとに A・B・C と戻すと面心立方になる。三角形の層は、黒鉛の層と同じ向き（面 x＋y＋z＝一定）にとった。",
"B の層は (6,6,0)・(6,0,6)・(0,6,6) の向きで A の層に乗る。長さ² は 72 で、立方体の枠の面の対角線そのもの。3 層目を同じ向きで乗せると C の層になり、面心立方の A と B の正四面体が続く。六方最密では 3 層目を (2,2,8)・(2,8,2)・(8,2,2) の向きで乗せる。長さ² は同じ 72 だが、立方体の枠の向きから外れ、3 層目は A の真上 (8,8,8) に来る。",
"どちらでも B の層の原子の隣は、同じ層 6・下 3・上 3 の 12 個。違いは正四面体のつながり方で、六方最密では B の層の三角形をはさんで、下向きと上向きの正四面体が面を分け合う 2 個組ができる。面心立方ではこの組は 0 個で、正四面体どうしは辺までしか分け合わない。"]
open('/home/claude/v288/bond30_hcp.html','w',encoding='utf-8').write(page("六方最密の金属","三角形の層を A・B・A と積んで六方最密を作り、3 層目を C に置く面心立方と比べる図。六方最密では B の層をはさんで正四面体が面を分け合う。ドラッグで回せる","六方最密の金属 ― 3 層目で枠の向きから外れる",ps,model,draw,2))
