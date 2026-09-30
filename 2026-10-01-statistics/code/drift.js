const {carrier,step}=require('./bell_core.js');
const c=carrier(); const T=30;
console.log('床',c.WH,c.HH);
const out=[];
for(const d of [0,3,6,9,12,15]){
 for(let a=0;a<10;a++){
  const ang=a*Math.PI/10; const tx=d*Math.cos(ang), ty=d*Math.sin(ang)*0.9;
  if(d===0&&a>0)break;
  let src=0,bd=1e9;for(let v=0;v<c.n;v++){const q=Math.hypot(c.X[v]-tx,c.Y[v]-ty);if(q<bd){bd=q;src=v;}}
  const xs=c.X[src],ys=c.Y[src];
  if(Math.abs(xs)+18>c.WH||Math.abs(ys)+14>c.HH) continue;   // 端から 4σ 以上離す
  let W=new Array(c.n).fill(0n);W[src]=1n;for(let t=0;t<T;t++)W=step(W,c.adj,new Uint8Array(c.n));
  let S=0n,Sx=0n,Sy=0n,edge=0n;
  for(let v=0;v<c.n;v++){const w=W[v];if(!w)continue;S+=w;Sx+=w*BigInt(Math.round((c.X[v]-xs)*1000));Sy+=w*BigInt(Math.round((c.Y[v]-ys)*1000));
    if(Math.abs(c.X[v])>c.WH-1||Math.abs(c.Y[v])>c.HH-1)edge+=w;}
  const mx=Number(Sx*1000n/S)/1e6, my=Number(Sy*1000n/S)/1e6;
  const r=Math.hypot(xs,ys); const ux=r?xs/r:0, uy=r?ys/r:0;
  const rad=mx*ux+my*uy, tan=-mx*uy+my*ux;   // rad<0 なら中心へ寄る
  out.push([d,a,r.toFixed(1),rad.toFixed(3),tan.toFixed(3),(edge*1000000n/S).toString(),c.adj[src].length]);
 }
}
console.log('設定距離 向き 実距離 半径方向のずれ(負=中心へ) 接線方向のずれ 端の百万分率 源の隣数');
for(const o of out)console.log(o.join('  '));
