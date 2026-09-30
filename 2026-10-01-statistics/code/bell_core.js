const PHI=(1+Math.sqrt(5))/2;
function penrose(n){let tris=[];for(let i=0;i<10;i++){let b=[Math.cos((2*i-1)*Math.PI/10),Math.sin((2*i-1)*Math.PI/10)];let c=[Math.cos((2*i+1)*Math.PI/10),Math.sin((2*i+1)*Math.PI/10)];if(i%2===0){const t=b;b=c;c=t;}tris.push([0,[0,0],b,c]);}
for(let k=0;k<n;k++){const out=[];for(const [col,A,B,C] of tris){if(col===0){const P=[A[0]+(B[0]-A[0])/PHI,A[1]+(B[1]-A[1])/PHI];out.push([0,C,P,B],[1,P,C,A]);}else{const Q=[B[0]+(A[0]-B[0])/PHI,B[1]+(A[1]-B[1])/PHI];const R=[B[0]+(C[0]-B[0])/PHI,B[1]+(C[1]-B[1])/PHI];out.push([1,R,C,A],[1,Q,R,B],[0,R,Q,A]);}}tris=out;}return tris;}
function carrier(){const tris=penrose(8);const e0=tris[0];const s=1/Math.hypot(e0[2][0]-e0[1][0],e0[2][1]-e0[1][1]);
const WH=Math.min(34,Math.floor(s*0.72)),HH=Math.round(WH*0.6);const idx=new Map(),X=[],Y=[],edges=new Set();
const inR=p=>Math.abs(p[0]*s)<=WH&&Math.abs(p[1]*s)<=HH;const vid=p=>{const k=Math.round(p[0]*s*1e3)+','+Math.round(p[1]*s*1e3);let i=idx.get(k);if(i===undefined){i=X.length;idx.set(k,i);X.push(p[0]*s);Y.push(p[1]*s);}return i;};
for(const [,A,B,C] of tris)for(const Q of [B,C]){if(!inR(A)||!inR(Q))continue;const a=vid(A),b=vid(Q);edges.add(a<b?a+'_'+b:b+'_'+a);}
const n=X.length,adj=Array.from({length:n},()=>[]),E=[];for(const e of edges){const [a,b]=e.split('_').map(Number);adj[a].push(b);adj[b].push(a);E.push([a,b]);}return {n,X,Y,adj,E,WH,HH};}
// 1刻：自分の数 ＋ 隣の数 を足す（その場に留まる道と、隣へ移る道）
function step(W,adj,blocked){const n=W.length,V=new Array(n);for(let v=0;v<n;v++){if(blocked[v]){V[v]=0n;continue;}let t=W[v];for(const w of adj[v])t+=W[w];V[v]=t;}return V;}
if(typeof module!=='undefined')module.exports={carrier,step};
