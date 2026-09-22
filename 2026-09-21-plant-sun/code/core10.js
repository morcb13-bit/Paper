// 対称にした成長オートマトン（検定SG3）
//  円環 i（生きている）：表面の桁の和を溜め、+3 で −5 して +1、−3 で +5 して −1 を、太陽に最も向いた子へ渡す
//  子の場所 j：受けた ±1 を平衡5進の一桁として溜める
//     +3 → −5：まだ無く親がすべて生きていれば生まれる（そうでなければ外へ）
//     −3 → +5：生きていれば消える（そうでなければ外へ）
//  種（1・2段目）は誰の子でもないので受けない
function makeSim6(GEO,S0){
  let S=S0; const R=GEO.R,N=R.length,P=GEO.P;
  const kids=R.map(()=>[]); R.forEach((r,i)=>r[3].forEach(p=>kids[p].push(i)));
  const ang=(ax,ay,bx,by)=>{let d=Math.abs(Math.atan2(ay,ax)-Math.atan2(by,bx))%(2*Math.PI);return Math.min(d,2*Math.PI-d)*180/Math.PI;};
  const dig=t0=>{const t=Math.round(t0*1e6)/1e6;return t<36?2:t<72?1:t<108?0:t<144?-1:-2;};
  const st={E:new Uint8Array(N),D:new Uint8Array(N),acc:new Int32Array(N),sl:new Int32Array(N),t:0,inp:0,
    fp:0,fn:0,bp:0,bn:0,births:0,deaths:0,outp:0,outn:0,splitP:0,splitN:0,tieFrames:new Set(),splitRings:[],ledger:true,fired:new Uint8Array(N),born:[],died:[]};
  for(let i=0;i<N;i++) if(R[i][0]<=1) st.E[i]=1;
  st.setSun=x=>{S=x;};
  st.step=function(){
    const cnt=new Int32Array(P.length);
    for(let i=0;i<N;i++) if(st.E[i]) for(const k of R[i][4]) cnt[k]++;
    st.fired.fill(0); st.born=[]; st.died=[]; st.splitRings=[];
    const tgt=(i,sg)=>{const cx=R[i][1],cy=R[i][2];const C=(sg>0?kids[i].filter(j=>!st.E[j]&&R[j][3].every(p=>st.E[p])):R[i][3].filter(p=>st.E[p]));
      let bd=1e9,B=[];for(const j of C){let d=S?Math.round(ang(R[j][1]-cx,R[j][2]-cy,S[0]-cx,S[1]-cy)*1e6)/1e6:0;if(sg<0)d=180-d;
        if(d<bd){bd=d;B=[j];}else if(d===bd)B.push(j);}
      if(B.length>1){st.tieFrames.add(st.t);if(sg>0)st.splitRings.push(i);if(sg>0)st.splitP+=B.length-1;else st.splitN+=B.length-1;}
      return B;};
    for(let i=0;i<N;i++){ if(!st.E[i]) continue;
      const cx=R[i][1],cy=R[i][2]; let s=0;
      const root=GEO.rootY!==undefined&&cy<GEO.rootY; s=0;
      if(!root){ s=1; if(S) for(const k of R[i][4]) if(cnt[k]===1) s+=dig(ang(P[k][0]-cx,P[k][1]-cy,S[0]-cx,S[1]-cy)); }
      else if(GEO.rootMode==='A'){ s=1; }
      else if(GEO.rootMode==='B'){ s=S?0:1; }
      else if(GEO.rootMode==='C'){ s=1; if(S) for(const k of R[i][4]) if(cnt[k]===1) s-=dig(ang(P[k][0]-cx,P[k][1]-cy,S[0]-cx,S[1]-cy)); }
      st.inp+=s; st.acc[i]+=s;
      while(st.acc[i]>=3){st.acc[i]-=5;st.fp++;st.fired[i]=1;const B=tgt(i,1);if(B.length)for(const j of B)st.sl[j]++;else st.outp++;}
      while(st.acc[i]<=-3){st.acc[i]+=5;st.fn++;st.fired[i]=2;const B=tgt(i,-1);if(B.length)for(const j of B)st.sl[j]--;else st.outn++;}
    }
    for(let j=0;j<N;j++){
      while(st.sl[j]>=3){st.sl[j]-=5;st.bp++; if(!st.E[j]&&R[j][3].every(p=>st.E[p])){st.E[j]=1;st.births++;st.born.push(j);} }
      while(st.sl[j]<=-3){st.sl[j]+=5;st.bn++; if(st.E[j]){st.E[j]=0;st.D[j]=1;st.deaths++;st.died.push(j);} }
    }
    let a=0,b=0; for(let i=0;i<N;i++){a+=st.acc[i];b+=st.sl[i];}
    st.ledger=st.ledger&&(st.inp===a+5*(st.fp-st.fn))&&((st.fp-st.fn-st.outp+st.outn+st.splitP-st.splitN)===b+5*(st.bp-st.bn)); st.t++;
  };
  return st;
}
if(typeof module!=='undefined') module.exports=makeSim6;
