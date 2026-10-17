s=open('/home/claude/v288/in/bond15_h2o.html',encoding='utf-8').read()
HEAD0=s[:s.index('<div class="wrap">')].replace('aspect-ratio:4/3','aspect-ratio:1/1').replace('--share:#7a3fd1;','--share:#7a3fd1; --at:#3d4744;').replace('--share:#b48cff; color-scheme:dark}}','--share:#b48cff; --at:#c9d2ce; color-scheme:dark}}').replace('--share:#b48cff; color-scheme:dark}','--share:#b48cff; --at:#c9d2ce; color-scheme:dark}')
def head(title): return HEAD0.replace('<title>シェルピンスキーの水</title>',f'<title>{title}</title>')
CTRL='''  <div class="ctrl">
    <button id="play" aria-pressed="true">再生中</button>
    <button id="step">次の段</button>
    <button id="spin" aria-pressed="true">回転</button>
    <label class="spd" for="spd">速さ <input id="spd" type="range" min="1" max="100" value="45"></label>
  </div>'''
HELP=r'''
const BITS=[];for(const a of[0,1])for(const b of[0,1])for(const c of[0,1])BITS.push([a,b,c]);
const odd=n=>(n[0]+n[1]+n[2])%2===1;
const tet=(o,e,wantOdd)=>BITS.filter(n=>odd(n)===wantOdd).map(n=>[o[0]+n[0]*e,o[1]+n[1]*e,o[2]+n[2]*e]);
const key=v=>v.join(",");
const add=(a,b)=>[a[0]+b[0],a[1]+b[1],a[2]+b[2]];
const dot=(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
function grow(seed,nbr,K){const GEN=new Map([[key(seed),0]]),AT=[{p:seed,g:0}],BD=[];let front=[seed];
  for(let k=1;k<=K;k++){const nxt=[];for(const p of front)for(const d of nbr(p)){const q=add(p,d),kq=key(q);
    if(!GEN.has(kq)){GEN.set(kq,k);AT.push({p:q,g:k});nxt.push(q);} if(GEN.get(kq)===k)BD.push({a:p,b:q,g:k});}front=nxt;}
  return{GEN,AT,BD,SH:[...Array(K+1).keys()].map(k=>AT.filter(a=>a.g===k).length)};}
'''
VIEW=r'''
const cv=document.getElementById("cv"), g=cv.getContext("2d");
const css=n=>getComputedStyle(document.documentElement).getPropertyValue(n).trim();
let yaw=0.6, pitch=0.38, spin=true, play=true, speed=0.45, stage=0, tS=0, last=performance.now(), SCALE=1, CEN=[0,0,0], UPAX=false;
const ease=t=>t<0?0:t>1?1:t*t*(3-2*t);
const RGBA={};function rgba(c,a){let h=RGBA[c];if(!h){const d=document.createElement("canvas").getContext("2d");d.fillStyle=c;h=RGBA[c]=d.fillStyle;}
  if(h[0]==="#")return`rgba(${parseInt(h.slice(1,3),16)},${parseInt(h.slice(3,5),16)},${parseInt(h.slice(5,7),16)},${a})`;return h;}
const S2=Math.SQRT2,S6=Math.sqrt(6),S3=Math.sqrt(3);
function proj(q){
  let p=[q[0]-CEN[0],q[1]-CEN[1],q[2]-CEN[2]];
  if(UPAX)p=[(p[0]-p[1])/S2,(p[0]+p[1]+p[2])/S3,(p[0]+p[1]-2*p[2])/S6];
  const cy=Math.cos(yaw), sy=Math.sin(yaw), cp=Math.cos(pitch), sp=Math.sin(pitch);
  const x1=cy*p[0]+sy*p[2], z1=-sy*p[0]+cy*p[2];
  const y2=cp*p[1]-sp*z1, z2=sp*p[1]+cp*z1;
  const f=1/(1+z2/180);
  return [cv.width/2+x1*SCALE*f, cv.height*0.55-y2*SCALE*f, z2, f];
}
const TE=[[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]], TF=[[0,1,2],[0,1,3],[0,2,3],[1,2,3]];
const CE=[];BITS.forEach((a,i)=>BITS.forEach((b,j)=>{if(j<=i)return;let d=0;for(let t=0;t<3;t++)if(a[t]!==b[t])d++;if(d===1)CE.push([i,j]);}));
function cubeLines(items,o,e,alpha,w){const V=BITS.map(n=>proj([o[0]+n[0]*e,o[1]+n[1]*e,o[2]+n[2]*e]));
  CE.forEach(([i,j])=>items.push({z:(V[i][2]+V[j][2])/2+0.05,k:"line",P:[V[i],V[j]],c:rgba(css("--cube"),alpha),w}));}
function tetDraw(items,T,col,solid,alpha,w,fillA,grow=1){
  const c=[0,1,2].map(i=>(T[0][i]+T[1][i]+T[2][i]+T[3][i])/4);
  const V=T.map(v=>proj(v.map((x,i)=>c[i]+(x-c[i])*grow)));
  TF.forEach(f=>items.push({z:(V[f[0]][2]+V[f[1]][2]+V[f[2]][2])/3,k:"face",P:f.map(i=>V[i]),c:rgba(col,fillA*alpha)}));
  TE.forEach(([i,j])=>items.push({z:(V[i][2]+V[j][2])/2-0.02,k:"line",P:[V[i],V[j]],c:rgba(col,0.9*alpha),w,dash:!solid}));}
function lineP(items,a,b,c,w,dash,dz=0){const p=proj(a),q=proj(b);items.push({z:(p[2]+q[2])/2-0.05+dz,k:"line",P:[p,q],c,w,dash});}
function atomP(tops,p,c,r){const q=proj(p);tops.push({z:q[2]-1,k:"atom",P:[q],c,r:r*q[3]});}
function render(items,tops){items.sort((p,q)=>q.z-p.z); tops.sort((p,q)=>q.z-p.z);
  for(const it of items.concat(tops)){const P=it.P;
    if(it.k==="line"){g.beginPath();g.moveTo(P[0][0],P[0][1]);g.lineTo(P[1][0],P[1][1]);g.strokeStyle=it.c;g.lineWidth=it.w;g.lineCap="round";if(it.dash)g.setLineDash([7,6]);g.stroke();g.setLineDash([]);}
    else if(it.k==="face"){g.beginPath();g.moveTo(P[0][0],P[0][1]);P.slice(1).forEach(p=>g.lineTo(p[0],p[1]));g.closePath();g.fillStyle=it.c;g.fill();}
    else if(it.k==="atom"){const q=P[0];g.beginPath();g.arc(q[0],q[1],Math.max(0.5,it.r),0,7);g.fillStyle=it.c;g.fill();g.lineWidth=1.2;g.strokeStyle=css("--panel");g.stroke();}}}
function texts(R,T0,T1,L2,L3,leg,names){
  g.textAlign="left";g.textBaseline="top";
  g.fillStyle=css("--fg");g.font=`700 ${R*1.6}px ${css("--f-body")}`;g.fillText(T0,R,R*0.9);
  g.fillStyle=css("--muted");g.font=`${R*1.1}px ${css("--f-body")}`;g.fillText(T1,R,R*2.9);
  g.fillStyle=css("--fg");g.font=`${R*1.1}px ${css("--f-mono")}`;if(L2)g.fillText(L2,R,R*4.5);if(L3)g.fillText(L3,R,R*5.9);
  let y=cv.height-R*1.45*leg.length+R*0.2;g.font=`${R*1.1}px ${css("--f-body")}`;g.textBaseline="middle";
  for(const[c,t,k]of leg){g.strokeStyle=css(c);g.fillStyle=css(c);
    if(k==="atom"){g.beginPath();g.arc(R*2,y,R*0.5,0,7);g.fill();}
    else{g.lineWidth=k==="thick"?R*0.26:2.5;if(k==="dash")g.setLineDash([6,5]);g.beginPath();g.moveTo(R,y);g.lineTo(R*3,y);g.stroke();g.setLineDash([]);}
    g.fillStyle=css("--fg");g.fillText(t,R*3.4,y);y+=R*1.45;}
  g.textAlign="right";g.textBaseline="top";
  names.forEach((n,i)=>{g.fillStyle=i===stage?css("--share"):css("--muted");g.fillText((i===stage?"● ":"○ ")+n,cv.width-R,R*(8+i*1.45));});}
function setStage(n){stage=(n+NST)%NST;tS=0;}
const pb=document.getElementById("play");
const stopUI=()=>{play=false;pb.setAttribute("aria-pressed",false);pb.textContent="止めている";};
pb.onclick=()=>{play=!play;pb.setAttribute("aria-pressed",play);pb.textContent=play?"再生中":"止めている";};
document.getElementById("step").onclick=()=>{stopUI();setStage(stage+1);tS=DUR[stage];};
document.getElementById("spin").onclick=e=>{spin=!spin;e.target.setAttribute("aria-pressed",spin);};
document.getElementById("spd").oninput=e=>{speed=e.target.value/100;};
let drag=null;
cv.addEventListener("pointerdown",e=>{drag=[e.clientX,e.clientY];cv.setPointerCapture(e.pointerId);});
cv.addEventListener("pointermove",e=>{if(!drag)return;yaw+=(e.clientX-drag[0])*0.01;pitch=Math.max(-1.2,Math.min(1.5,pitch+(e.clientY-drag[1])*0.01));drag=[e.clientX,e.clientY];});
cv.addEventListener("pointerup",()=>drag=null);
function tick(t){
  const dt0=Math.min(0.05,(t-last)/1000); last=t;
  const k=0.3+speed*1.7, dt=dt0*k;
  if(spin) yaw+=dt0*0.25*k;
  tS=Math.min(tS+dt, play?1e9:DUR[stage]);
  if(play&&tS>DUR[stage]) setStage(stage+1);
  draw(); requestAnimationFrame(tick);
}
'''
def page(title,aria,h1,ps,model,draw,reduced):
    body='<div class="wrap">\n  <canvas id="cv" width="1000" height="1000" aria-label="'+aria+'"></canvas>\n'+CTRL+'\n  <h1>'+h1+'</h1>\n'+''.join('  <p>'+p+'</p>\n' for p in ps)+'</div>\n'
    js='<script>\n// ---- 整数の模型 ----'+HELP+model+'\n// ---- ここから先は観察者側の表示（回転・透視・補間だけ浮動小数） ----'+VIEW+draw+'\nif(matchMedia("(prefers-reduced-motion: reduce)").matches){spin=false;document.getElementById("spin").setAttribute("aria-pressed",false);stopUI();setStage('+str(reduced)+');tS=DUR[stage];}\nrequestAnimationFrame(tick);\n</script>\n'
    return head(title)+body+js
