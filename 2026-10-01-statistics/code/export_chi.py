import json,runpy,sys
sys.argv=['x','30','20000']
g=runpy.run_path('chi_penrose.py')
X,Y,adj,src,land,shell,RMAX,deg,shuf,n=[g[k] for k in ['X','Y','adj','src','land','shell','RMAX','deg','shuf','n']]
inner=[v for v in range(n) if shell[v]<=RMAX]
pos={v:i for i,v in enumerate(inner)}
V=[[round(X[v]-X[src],3),round(Y[v]-Y[src],3),deg[v],shell[v],land[v],shuf[v]] for v in inner]
E=[[pos[a],pos[b]] for a in inner for b in adj[a] if b in pos and pos[a]<pos[b]]
json.dump({'V':V,'E':E},open('chi_data.json','w'),separators=(',',':'))
print(len(V),len(E),sum(v[4] for v in V))
