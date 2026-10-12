import sys
from PIL import Image, ImageDraw, ImageFont
R='/home/claude/mol/images/fib-shells/raw/'; O='/home/claude/mol/images/fib-shells/'
FONT=ImageFont.truetype('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc',26,index=0)
def grid(items,cols,out,cell=400,title=None):
    rows=(len(items)+cols-1)//cols; cap=44; top=56 if title else 0
    W=cols*cell+(cols+1)*12; H=top+rows*(cell+cap)+(rows+1)*12
    im=Image.new('RGB',(W,H),'white'); d=ImageDraw.Draw(im)
    if title: d.text((16,12),title,fill=(28,35,32),font=FONT)
    for k,(f,lab) in enumerate(items):
        r,c=divmod(k,cols); x=12+c*(cell+12); y=top+12+r*(cell+cap+12)
        p=Image.open(R+f).convert('RGB').resize((cell,cell)); im.paste(p,(x,y))
        d.rectangle([x,y,x+cell-1,y+cell-1],outline=(205,211,206))
        d.text((x+6,y+cell+6),lab,fill=(28,35,32),font=FONT)
    im.save(O+out); print(out,im.size)
grid([("td_1.png","H：辺 1"),("td_1m.png","He：辺 1 と鏡像"),("td_2.png","L：辺 2"),("td_3.png","M：辺 3")],4,"fig1_edge_n.png")
grid([("fs_1.png","1  H  辺 1（正立）"),("fs_2.png","2  He  辺 1（鏡像）"),("fs_3.png","3  L  辺 2（正立）"),("fs_4.png","4  M  辺 3（鏡像）"),
      ("fs_5.png","5  N  辺 5（正立）"),("fs_6.png","6  O  辺 8（鏡像）"),("fs_7.png","7  P  辺 13（正立）"),("fs_8.png","8  Q  辺 21（鏡像）")],4,"fig2_alternate.png")
grid([("bv_bohr.png","ボーア模型：円軌道（4 : 9）"),("bv_face.png","B13：席は正四面体の網の点"),("bv_sphere.png","B13：席を二つの球へ写す")],3,"fig3_bohr_b13.png")
