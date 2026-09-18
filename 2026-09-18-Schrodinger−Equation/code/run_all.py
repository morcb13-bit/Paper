#!/usr/bin/env python3
"""run_all.py --- 検定をすべて走らせ、OK/NG の合計を出す

  python3 run_all.py

標準ライブラリだけで走る。外部依存はない。
どれか一つでも NG が出たら終了コードが 1 になる。
"""

import subprocess
import sys
import re
from pathlib import Path

# NG が結果そのものである検定。件数と理由を明示しておく
EXPECTED_NG = {
    "crossing3.py": (1, "検定F3 は「(A, Q) では閉じない」ことを示す検定。"
                        "閉じ方は crossing4.py の (S, O) が与える"),
}

ORDER = [
    ("b13_check.py", "平衡13進の桁列と Z[φ] の整数対"),
    ("pentagram_step.py", "正十二面体の分類と一歩の作用素"),
    ("crossing2.py", "c(k) の二次の形と三領域"),
    ("crossing3.py", "A' = L^2 A + c(k) Q と固有空間の分解"),
    ("crossing4.py", "(S, O) の閉じ方と二階漸化式"),
    ("nonregular.py", "次数が不揃いなグラフでの同じ式"),
    ("kernel.py", "階数39と核81"),
    ("conic.py", "対に乗る二次形式と曲線の種類"),
]

here = Path(__file__).resolve().parent
total_ok = total_all = 0
failed = []

for name, desc in ORDER:
    path = here / name
    if not path.exists():
        print(f"見つからない: {name}")
        failed.append(name)
        continue
    r = subprocess.run([sys.executable, str(path)], capture_output=True,
                       text=True, cwd=str(here))
    tail = [ln for ln in r.stdout.strip().splitlines() if re.match(r"^\d+/\d+ OK$", ln)]
    if r.returncode != 0 or not tail:
        print(f"NG  {name:<20} 走らなかった")
        failed.append(name)
        continue
    ok, all_ = (int(x) for x in tail[-1].split()[0].split("/"))
    total_ok += ok
    total_all += all_
    exp, why = EXPECTED_NG.get(name, (0, ""))
    good = (ok == all_ - exp)
    mark = "OK " if good else "NG "
    print(f"{mark} {name:<20} {ok}/{all_}   {desc}")
    if exp:
        print(f"{'':4}{'':20} 予定された NG {exp}件 ── {why}")
    if not good:
        failed.append(name)

print()
exp_total = sum(e for e, _ in EXPECTED_NG.values())
print(f"合計 {total_ok}/{total_all}（うち予定された NG {exp_total}件）")
sys.exit(1 if failed else 0)
