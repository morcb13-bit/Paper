function makeModel(G,walls,ignore){
  const N=G.rings.length, wall=new Set(walls);
  const open=i=>!wall.has(i), spreadOK=i=>ignore||open(i);
  let D=new Array(N).fill(0n), p4=1n, t=0;
  const dist=new Array(N).fill(-1); dist[G.food]=0; const q=[G.food];
  while(q.length){const i=q.shift();for(const j of G.rings[i].n) if(open(j)&&dist[j]<0){dist[j]=dist[i]+1;q.push(j);}}
  function step(){p4*=4n;const E=new Array(N).fill(0n);
    for(let i=0;i<N;i++) if(spreadOK(i)){let s=0n;for(const j of G.rings[i].n) if(spreadOK(j)) s+=D[j];E[i]=s;}
    E[G.food]+=p4; D=E; t++;}
  function climb(start){ // 同点は二つとも
    let cur=new Set([start]); const layers=[[start]], edges=[];
    for(let k=0;k<200;k++){
      if(cur.size===1&&cur.has(G.food)) return {ok:true,steps:k,layers,edges};
      const nxt=new Set();
      for(const c of cur){ if(c===G.food){nxt.add(c);continue;}
        const nb=G.rings[c].n.filter(open); let m=-1n; for(const j of nb) if(D[j]>m) m=D[j];
        if(m<=D[c]) return {ok:false,steps:k,layers,edges,stopRow:G.rings[c].row};
        for(const j of nb) if(D[j]===m){nxt.add(j);edges.push([c,j,k]);} }
      cur=nxt; layers.push([...cur]); }
    return {ok:false,steps:200,layers,edges};}
  return {get D(){return D},get t(){return t},step,climb,dist,wall,open};
}
if(typeof module!=='undefined') module.exports=makeModel;
