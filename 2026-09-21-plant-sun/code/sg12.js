// 検定SG12  10枚の担体で昼と夜。1日20コマ、太陽は18°ずつ一周するが、光るのは昼の半周（180°→90°→0°、11コマ）だけ。
// 夜（9コマ）は太陽の桁が0、円環の +1 だけが続く。規則は SG11 と同じ。
// 判定：東→西／西→東が鏡映で一致／帳尻／SG11（一周とも光る）との違いを数える
const G=require('./geo10.json'),mk=require('./core9.js');const [cx0,cy0]=G.center;
const RAD=Math.max(...G.P.map(p=>Math.hypot(p[0]-cx0,p[1]-cy0)))+22;
const pos=d=>[cx0+RAD*Math.cos(d*Math.PI/180),cy0+RAD*Math.sin(d*Math.PI/180)];
const at={};const key=(x,y)=>Math.round(x*100)+','+Math.round(y*100);G.R.forEach((r,i)=>at[key(r[1],r[2])]=i);
const mir=i=>at[key(2*cx0-G.R[i][1],G.R[i][2])];
const EW=t=>{const f=t%20;return f<=10?pos(180-18*f):null;}, WE=t=>{const f=t%20;return f<=10?pos(18*f):null;};
function run(sun,T){const s=mk(G,null);const hist=[];let lb=-1;for(let t=0;t<T;t++){s.setSun(sun(t));s.step();if(s.born.length)lb=t;const E=[];s.E.forEach((v,i)=>v&&E.push(i));hist.push(E);}return {s,hist,lb};}
const T=20*40;const a=run(EW,T),b=run(WE,T);
const rows=r=>{const d=Array(14).fill(0);r.s.E.forEach((v,i)=>{if(v)d[G.R[i][0]]++;});return d.join(' ');};
for(const [k,r] of [['東→西',a],['西→東',b]]){const per=Array(10).fill(0);r.s.E.forEach((v,i)=>{if(v)per[G.R[i][5]]++;});
 console.log(k,`生きている${per.reduce((x,y)=>x+y)} 生まれた${r.s.births} 消えた${r.s.deaths} 二葉${r.s.splitP} 帳尻${r.s.ledger} 最後に生まれたコマ${r.lb} 枚ごと${per.join('/')}\n   段ごとの生きている環 ${rows(r)}`);}
let first=-1;for(let t=0;t<T;t++){if(a.hist[t].map(mir).sort((x,y)=>x-y).join()!==b.hist[t].join()){first=t;break;}}
console.log('鏡映が崩れるコマ',first);
for(const d of [1,5,10,20,40]){const E=a.hist[d*20-1];console.log(`  ${d}日目の終わり 生きている${E.length}`);}
