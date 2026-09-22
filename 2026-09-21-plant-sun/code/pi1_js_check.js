// 動く図の計算（logic_walls.js）が PI1 と同じ結果になるかの確認。ジグザグも含む
const G=require('./wing0.json'),mk=require('./logic_walls.js');
const RG=(a,b)=>{const o=[];for(let i=a;i<=b;i++)o.push(i);return o;};
const PRE={one:RG(30,35),right:RG(28,33),zig:[...RG(10,12),...RG(30,35),...RG(55,63)],zig2:[...RG(12,14),...RG(28,33),...RG(57,65)]};
for(const k in PRE)for(const ig of [false,true]){const m=mk(G,PRE[k],ig);const T=Math.max(100,3*Math.max(...m.dist));for(let t=0;t<T;t++)m.step();
 console.log(k,ig?'壁を無視':'道に沿う',G.seed.map(s=>{const r=m.climb(s);return (r.ok?'着く':'止まる段'+r.stopRow)+' '+r.steps+'/'+m.dist[s]}).join(' | '));}
