import json
def load(f): return json.load(open(f))
def test(a,b,key):
    d=[x-y for x,y in zip(a[key],b[key])]; S=sum(d); Q=sum(v*v for v in d)
    return S, S*S>9*Q
rows=[('球','w1_fan.json','w1_off.json'),('翼＋','w1_fan_wing+.json','w1_off_wing+.json'),('翼−','w1_fan_wing-.json','w1_off_wing-.json')]
for nm,f,o in rows:
    a,b=load(f),load(o); out=[]
    for key,lab in (('D','抗力ΔD'),('Ly','上下ΔLy'),('Lz','横ΔLz')):
        S,big=test(a,b,key); out.append(f"{lab}={S:+d}{'（雑音越え）' if big else '（雑音の内）'}")
    print(nm,'帳簿',a['ok'] and b['ok'],' '.join(out))
