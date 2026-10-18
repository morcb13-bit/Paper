import numpy as np, json, base64
phi=(1+5**0.5)/2
Z=np.load('mbox43.npz'); x=Z['xpar']; xf=x[...,0]+x[...,1]*phi   # 表示用の座標
N=69; lo=-43; h=86/N
slab=np.abs(xf[:,2])<6
ix=np.clip(((xf[:,0]-lo)/h).astype(int),0,N-1); iy=np.clip(((xf[:,1]-lo)/h).astype(int),0,N-1)
cnt=np.zeros((N,N)); np.add.at(cnt,(ix[slab],iy[slab]),1)
def field(body,mode):
    S=np.load('vfield_%s_%s.npz'%(mode,body))['snap']            # (20,V,3,2)
    S=S[...,0]+S[...,1]*phi                                     # 表示用
    return S
def smooth(a):
    p=np.pad(a,1,mode='edge'); return sum(p[i:i+N,j:j+N] for i in range(3) for j in range(3))/9
out={}; gmax=0; tmp={}
for body in ['ball','wing+','wing-']:
    F=field(body,'fan')-field(body,'off')                       # 風あり − 風なし
    frames=[]
    for k in range(2,20):
        m=F[k-2:k+1].sum(0)/300.0                               # 300歩の平均（頂点あたり・1歩あたり）
        gx=np.zeros((N,N)); gy=np.zeros((N,N))
        np.add.at(gx,(ix[slab],iy[slab]),m[slab,0]); np.add.at(gy,(ix[slab],iy[slab]),m[slab,1])
        gx=smooth(gx/np.maximum(cnt,1)); gy=smooth(gy/np.maximum(cnt,1))
        frames.append((gx,gy))
    tmp[body]=frames
allv=np.abs(np.array([[f[0],f[1]] for fr in tmp.values() for f in fr[-8:]]))
gmax=float(np.percentile(allv,99))
for body,frames in tmp.items():
    arr=np.array([[f[0],f[1]] for f in frames])                 # (18,2,N,N)
    q=np.clip(np.round(arr/gmax*127),-127,127).astype(np.int8)
    out[body]=base64.b64encode(q.tobytes()).decode()
meta=dict(N=N,lo=lo,h=h,frames=18,scale=gmax,steps=[(k+1)*100 for k in range(2,20)])
json.dump(dict(meta=meta,data=out),open('fig2_data.json','w'))
print('gmax',gmax,'bytes',len(json.dumps(out)))
