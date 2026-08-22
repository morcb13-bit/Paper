#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""エルデシュ問題から、B13 の道具と重なる未解決問題を絞り込む。

二つの公開リポジトリを引く（どちらも github.com 配下・許可済み）。

  teorth/erdosproblems              … 状態・タグ・OEIS。本文は入っていない
  google-deepmind/formal-conjectures … 問題文（Lean。informal もコメントに入る）

絞りの条件は三つ。順番に効かせる。

  1. status が open 系（open / falsifiable / verifiable / decidable）
  2. 文が Lean で形式化済み  ← #728 で揉めた「対象の確定」が済んでいる印
  3. 本文に B13 側の語が入る（フィボナッチ・φ・完全性・準周期・タイリング…）

条件2 を外さないこと。外すと、何を解いたことになるのかを自分で決める仕事が
戻ってくる。それは #728 の教訓のちょうど裏返しになる。

使い方:
    python3 shortlist.py            # 手元に clone 済みなら再利用
    python3 shortlist.py --fetch    # clone し直す
"""

import os
import subprocess
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(BASE, "erdosproblems")
FC = os.path.join(BASE, "formal-conjectures")

REPOS = [
    (DB, "https://github.com/teorth/erdosproblems.git"),
    (FC, "https://github.com/google-deepmind/formal-conjectures.git"),
]

OPEN_STATES = {"open", "falsifiable", "verifiable", "decidable"}

# B13 側の語。増やしてよいが、増やしたら当たりが増えるのは当然なので、
# 増やした語をここに記録して、どの語で当たったかを必ず出力に残すこと。
KEYWORDS = [
    "fibonacci", "lucas number", "golden", "zeckendorf", "beatty",
    "complete", "quasi", "penrose", "pisano", "continued fraction",
    "irrational", "pentagon", "tiling", "aperiodic", "recurrence",
]


def fetch():
    for path, url in REPOS:
        if os.path.isdir(path):
            print(f"  既にある: {path}")
            continue
        print(f"  clone: {url}")
        subprocess.run(["git", "clone", "--depth", "1", url, path], check=True)


def load_status():
    try:
        import yaml
    except ImportError:
        sys.exit("PyYAML が要る:  pip install pyyaml --break-system-packages")
    p = os.path.join(DB, "data", "problems.yaml")
    return {q["number"]: q for q in yaml.safe_load(open(p))}


def scan():
    db = load_status()
    d = os.path.join(FC, "FormalConjectures", "ErdosProblems")
    rows = []
    for fn in os.listdir(d):
        if not fn.endswith(".lean"):
            continue
        num = fn[:-5]
        text = open(os.path.join(d, fn), encoding="utf-8").read().lower()
        hit = [w for w in KEYWORDS if w in text]
        if not hit:
            continue
        q = db.get(num, {})
        st = q.get("status", {}).get("state", "?")
        if st not in OPEN_STATES:
            continue
        rows.append((int(num), st, q.get("prize", "no"),
                     ",".join(q.get("tags", [])), ",".join(hit)))
    return sorted(rows)


def main():
    if "--fetch" in sys.argv:
        fetch()
    if not (os.path.isdir(DB) and os.path.isdir(FC)):
        sys.exit("リポジトリが無い。--fetch を付けて走らせること。")
    rows = scan()
    print(f"未解決 × 文が形式化済み × B13 の語  … {len(rows)} 件")
    print()
    print(f"{'#':>6} {'状態':<12} {'賞金':<8} {'当たった語'}")
    for n, st, pz, tags, hit in rows:
        print(f"{n:>6} {st:<12} {pz:<8} {hit}")
        if tags:
            print(f"{'':>6} {'':<12} {'':<8} tags: {tags}")


if __name__ == "__main__":
    main()
