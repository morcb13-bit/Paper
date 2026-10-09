const h=require('fs').readFileSync('chiral_three.html','utf8');
const s=h.split('<script>')[1];
const core=s.split('// ── 盤 ──')[0].replace('let extra = 1','var extra = 1');
eval(core+`
const P=[[1,0,0],[1,1,0],[1,1,1]];
for (const e of [1,-1,0]) { extra=e; const out=[];
 for (const [c,p,st] of P){ let x=seed(), snap=null, n=0;
  for(;n<30;n++){ if(c)x=copyStep(x); if(p)x=pairStep(x)[0]; if(st)x=stirStep(x)[0];
   const [a,b]=count(x); if(st?(a+b===60||a+b===0):(snap&&snap.join()===x.join()))break; snap=x.slice(); }
  out.push(count(x).join(',')+' 歩'+(n+1)); }
 console.log(e,out.join(' | ')); }`);
