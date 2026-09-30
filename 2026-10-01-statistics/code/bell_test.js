const {carrier,step}=require('./bell_core.js');
const c=carrier();let src=0,bd=1e9;for(let v=0;v<c.n;v++){const d=Math.hypot(c.X[v],c.Y[v]);if(d<bd){bd=d;src=v;}}
const blocked=new Uint8Array(c.n);let W=new Array(c.n).fill(0n);W[src]=1n;
for(let T=1;T<=90;T++){W=step(W,c.adj,blocked);
 if(T%30===0){let S=0n,Sx2=0n;for(let v=0;v<c.n;v++){const x=BigInt(Math.round((c.X[v]-c.X[src])*1000));S+=W[v];Sx2+=W[v]*x*x;}
  const s2=Number(Sx2*1000n/S)/1000/1e6, sg=Math.sqrt(s2);
  let inside=0n;for(let v=0;v<c.n;v++)if(Math.abs(c.X[v]-c.X[src])<=sg)inside+=W[v];
  // x の帯ごとの和
  const bins={};for(let v=0;v<c.n;v++){const b=Math.round((c.X[v]-c.X[src])/2);bins[b]=(bins[b]||0n)+W[v];}
  let mx=0n;for(const k in bins)if(bins[k]>mx)mx=bins[k];
  const row=[];for(let b=-8;b<=8;b++)row.push(Number(((bins[b]||0n)*9n)/mx));
  let edge=0n;for(let v=0;v<c.n;v++)if(Math.abs(c.X[v])>c.WH-1||Math.abs(c.Y[v])>c.HH-1)edge+=W[v];
  console.log('T',T,'sigma',sg.toFixed(2),'1σ内 百万分率',(inside*1000000n/S).toString(),'端の割合',(edge*1000000n/S).toString(),'形',row.join(''));}}
