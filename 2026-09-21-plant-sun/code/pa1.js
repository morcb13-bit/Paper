// 検定PA1  目的地だけの地図（壁なし）。規則は core10 のまま、太陽＝目的地（昼11コマだけ、夜9コマは無し）
// 事前登録
//   基準  目的地の祖先（円錐）＝ 種から目的地へ届く全経路の和。この担体では一歩ごとに段が一つ下がるので
//         歩数ではどの経路も同じ長さ → 歩数の最短路の和集合＝円錐。長さ（短 a 歩＋長 b 歩＝a+bφ）の最短は整数対で比べる
//   OK    目的地が生まれる、かつ 40日目に目的地の翼で生きている環がすべて円錐の中
//   NG    目的地に届かない／円錐の外に生き残る／全部埋まる
//   負の対照  太陽なし（全部埋まるはず）
const G=require('./geo10.json'),mk=require('./core10.js');const R=G.R,N=R.length;
const kids=R.map(()=>[]);R.forEach((r,i)=>r[3].forEach(p=>kids[p].push(i)));
const w0=[...Array(N).keys()].filter(i=>R[i][5]===0);
const bottom=w0.filter(i=>R[i][0]===13).sort((a,b)=>R[a][1]-R[b][1]);
function cone(g){const s=new Set([g]),q=[g];while(q.length){const i=q.pop();for(const p of R[i][3])if(!s.has(p)){s.add(p);q.push(p);}}return s;}
function run(goal,T,useSun){const s=mk(G,null);const n=[];let reach=-1;
 for(let t=0;t<T;t++){s.setSun(useSun&&t%20<=10?[R[goal][1],R[goal][2]]:null);s.step();
  if(reach<0&&s.E[goal])reach=t; if(t%20===19)n.push(w0.filter(i=>s.E[i]).length);}
 return {s,n,reach};}
for(const k of [0,6,13]){const g=bottom[k],C=cone(g),r=run(g,800,true);
 const alive=w0.filter(i=>r.s.E[i]),out=alive.filter(i=>!C.has(i));
 console.log(`目的地 翼0・段13の左から${k} 円錐${C.size}環  届いたコマ${r.reach} 生${alive.length} 円錐の外${out.length} 生${r.s.births}消${r.s.deaths} 帳尻${r.s.ledger}`);
 console.log('   翼0の生（1日ごと）',r.n.slice(0,20).join(' '));
 const other=R.map((x,i)=>i).filter(i=>R[i][5]!==0&&r.s.E[i]).length; console.log('   他の翼の生',other);}
const z=run(bottom[6],800,false);console.log('負の対照 太陽なし 翼0の生',w0.filter(i=>z.s.E[i]).length,'/105');
