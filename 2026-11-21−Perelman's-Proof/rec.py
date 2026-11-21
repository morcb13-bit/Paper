# Record a GIF from an HTML figure via window.renderFrame / window.FRAME_DUR (canvas pixels)
import sys, base64, io, asyncio
from PIL import Image
from playwright.async_api import async_playwright
src, out = sys.argv[1], sys.argv[2]; colors=int(sys.argv[3]) if len(sys.argv)>3 else 48; mode=sys.argv[4] if len(sys.argv)>4 else "canvas"
async def main():
    async with async_playwright() as p:
        b=await p.chromium.launch(executable_path="/usr/bin/google-chrome",args=["--no-sandbox"])
        pg=await b.new_page(viewport={"width":(860 if mode=="card" else 960),"height":900},device_scale_factor=1)
        await pg.goto("file://"+src+"?snap=0"); await pg.wait_for_function("window.FRAME_DUR!==undefined",timeout=120000)
        durs=await pg.evaluate("window.FRAME_DUR"); fr=[]
        for k in range(len(durs)):
            if mode=="card":
                await pg.evaluate(f"window.renderFrame({k})"); png=await pg.locator(".card").screenshot()
                fr.append(Image.open(io.BytesIO(png)).convert("RGB"))
            else:
                u=await pg.evaluate(f"(window.renderFrame({k}),document.getElementById('cv').toDataURL('image/png'))")
                fr.append(Image.open(io.BytesIO(base64.b64decode(u.split(',')[1]))).convert("RGB"))
        info=await pg.evaluate("window.SNAPS||null"); await b.close()
        pal=[f.convert("P",palette=Image.ADAPTIVE,colors=colors) for f in fr]
        pal[0].save(out,save_all=True,append_images=pal[1:],duration=durs,loop=0,optimize=True)
        for k in (0,len(fr)//3,len(fr)-17 if len(fr)>20 else 0,len(fr)-1): fr[k].save(out.replace(".gif",f"_chk{k}.png"))
        if info: print("first",info[0],"\nlast",info[-1])
        print(len(fr),"frames")
asyncio.run(main())
