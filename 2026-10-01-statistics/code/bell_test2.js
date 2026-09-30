const {carrier,step}=require('./bell_core.js');
const c=carrier();
for(const [sx,sy] of [[0,0],[5.3,2.1]]){
let src=0,bd=1e9;for(let v=0;v<c.n;v++){const d=Math.hypot(c.X[v]-sx,c.Y[v]-sy);if(d<bd){bd=d;src=v;}}
const blocked=new Uint8Array(c.n);let W=new Array(c.n).fill(0n);W[src]=1n;
for(let T=1;T<=60;T++)W=step(W,c.adj,blocked);
let mx=0n;for(const w of W)if(w>mx)mx=w;
let S=0n,Sr2=0n;for(let v=0;v<c.n;v++){const r2=BigInt(Math.round(((c.X[v]-c.X[src])**2+(c.Y[v]-c.Y[src])**2)*1000));S+=W[v];Sr2+=W[v]*r2;}
const s2=Number(Sr2*1000n/S)/1e6/2; // 一軸あたり
const rows={};for(let v=0;v<c.n;v++){const r=Math.hypot(c.X[v]-c.X[src],c.Y[v]-c.Y[src]);const b=Math.floor(r);const q=Number(W[v]*100000n/mx)/100000/Math.exp(-r*r/(2*s2));(rows[b]=rows[b]||[]).push(q);}
console.log('src',sx,sy,'deg',c.adj[src].length,'sigma',Math.sqrt(s2).toFixed(2));
for(let b=0;b<=14;b++){const a=rows[b]||[];if(!a.length)continue;console.log(b,a.length,Math.min(...a).toFixed(2),Math.max(...a).toFixed(2));}
// x 帯（幅1）の和を、帯の番地数で割らずに・割って
}
{let src=0,bd=1e9;for(let v=0;v<c.n;v++){const d=Math.hypot(c.X[v]-5.3,c.Y[v]-2.1);if(d<bd){bd=d;src=v;}}
let W=new Array(c.n).fill(0n);W[src]=1n;for(let T=1;T<=60;T++)W=step(W,c.adj,new Uint8Array(c.n));
let mx=0n;for(const w of W)if(w>mx)mx=w;let S=0n,Sr2=0n;for(let v=0;v<c.n;v++){const r2=BigInt(Math.round(((c.X[v]-c.X[src])**2+(c.Y[v]-c.Y[src])**2)*1000));S+=W[v];Sr2+=W[v]*r2;}
const s2=Number(Sr2*1000n/S)/1e6/2;const g={};
for(let v=0;v<c.n;v++){const r=Math.hypot(c.X[v]-c.X[src],c.Y[v]-c.Y[src]);if(r>14)continue;const q=Number(W[v]*100000n/mx)/100000/Math.exp(-r*r/(2*s2));const d=c.adj[v].length;(g[d]=g[d]||[]).push(q);}
for(const d in g){const a=g[d].sort((x,y)=>x-y);console.log('次数',d,'個数',a.length,'min',a[0].toFixed(2),'中央',a[a.length>>1].toFixed(2),'max',a[a.length-1].toFixed(2));}}
