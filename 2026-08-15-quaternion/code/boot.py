#  セッションの立ち上げ ── 一行で済ませる
#      python3 boot.py            … リポジトリを取ってきて道具を一箇所に集め、自己検査まで走らせる
#  やること：clone/pull → *.py を作業ディレクトリへ集める → b13_chain_units の検定7本
#            → 親が欠けているスクリプトを名指しする
import os, subprocess, sys, re, glob

W = os.path.expanduser("~/w"); os.makedirs(W, exist_ok=True); os.chdir(W)
if not os.path.isdir("Paper/.git"):
    subprocess.run(["git", "clone", "--depth", "1",
                    "https://github.com/morcb13-bit/Paper"], check=True)
else:
    subprocess.run(["git", "-C", "Paper", "pull", "--depth", "1", "-f"], check=False)

n = 0
for p in glob.glob("Paper/**/*.py", recursive=True):
    subprocess.run(["cp", p, W]); n += 1
print("道具 %d 本を集めた" % n)

r = subprocess.run([sys.executable, "b13_chain_units.py", "sieve.svg"],
                   capture_output=True, text=True)
tail = [l for l in r.stdout.splitlines() if l.startswith("NG")]
print("自己検査：%s" % (tail[0] if tail else "走らなかった → " + r.stderr.strip()[:200]))

need = set()
for p in glob.glob("*.py"):
    for m in re.finditer(r"open\('([A-Za-z0-9_]+\.py)'\)", open(p, encoding='utf-8', errors='replace').read()):
        if not os.path.exists(m.group(1)):
            need.add((p, m.group(1)))
if need:
    print("親が欠けている：")
    for a, b in sorted(need):
        print("   %s ← %s" % (a, b))
else:
    print("欠けている親はない")
