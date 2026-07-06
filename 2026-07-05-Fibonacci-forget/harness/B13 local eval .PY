#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
b13_local_eval.py — フィボナッチ忘却のハーネス検証 (BpVM規律)

写像M (判定対象):
  「フィボナッチ間隔の保持は、同予算の他の忘却則より
    研究想起の性能/コスト比で優れる」

事前登録 (FROZEN — 結果を見る前に固定。変更したら再凍結として明記):
  - 主対抗サロゲート: DOUBLE (0,1,2,4,8,16,... 同じ対数成長)
  - FORWARD: FIBの総合ヒット率 − DOUBLEの総合ヒット率 >= +MARGIN
             かつ FIB > RANDOM平均 (健全性)
  - ABSTAIN: |FIB − DOUBLE| < MARGIN (φ固有の優位なし、写像Mの棄却)
             または FIB <= RANDOM平均 (構造が偶然以下)
  - BLIND:   いずれかのビンの質問数 < MIN_PER_BIN
  - 判定は検索段階のみ: 標的版のチャンクが選択コンテキストに入ったか
  - ビン: 標的の新しい順インデックスで recent(0-3) / mid(4-20) / old(21+)

使い方:
  python3 b13_local_eval.py --make-qa            QA候補の雛形をDBから生成
  python3 b13_local_eval.py --freeze qa.json     QAセットをハッシュ凍結
  python3 b13_local_eval.py --run qa.json        評価実行 (凍結必須)
  python3 b13_local_eval.py --run qa.json --budget 6000 --top-k 8
"""
import argparse
import hashlib
import json
import os
import random
import sys

import psycopg2

# b13_local.py と同じディレクトリに置くこと
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import b13_local as core  # noqa: E402

# ============================ FROZEN (事前登録) ============================
MARGIN = 0.10          # FIB−DOUBLE がこれ以上で FORWARD
MIN_PER_BIN = 5        # 各ビンの最少質問数 (未満なら BLIND)
RANDOM_SEEDS = 5       # RANDOM サロゲートのシード数
BINS = [('recent', 0, 3), ('mid', 4, 20), ('old', 21, 10**9)]
# ==========================================================================


# ---------------------------------------------------------------- 保持則
def idx_fib(n):
    keep, a, b = {0}, 1, 2
    while a < n:
        keep.add(a)
        a, b = b, a + b
    return sorted(keep)


def idx_double(n):
    keep, a = {0}, 1
    while a < n:
        keep.add(a)
        a *= 2
    return sorted(keep)


def idx_uniform(n, k):
    if k >= n:
        return list(range(n))
    step = (n - 1) / (k - 1)
    return sorted({round(i * step) for i in range(k)})


def strategies(n):
    fib = idx_fib(n)
    k = len(fib)
    st = {
        'FULL': list(range(n)),
        'FIB': fib,
        'DOUBLE': idx_double(n),
        'UNIFORM': idx_uniform(n, k),
        'LATEST': list(range(min(k, n))),
    }
    for s in range(1, RANDOM_SEEDS + 1):
        rng = random.Random(s)
        st[f'RANDOM{s}'] = sorted(rng.sample(range(n), min(k, n)))
    return st


# ---------------------------------------------------------------- 検索(版付き)
def retrieve_meta(rows_subset, question, top_k, budget):
    """core と同じスコアリングで、選択チャンクの版番号を返す。"""
    strong, bigrams = core.query_features(question)
    scored = []
    for hid, ver, content, created, fn in rows_subset:
        if not content:
            continue
        for ch in core.split_chunks(content):
            s = core.score_chunk(ch, strong, bigrams)
            if s > 0:
                scored.append((s, str(ver), ch))
    scored.sort(key=lambda x: -x[0])
    picked, used = [], 0
    for s, ver, ch in scored:
        if len(picked) >= top_k:
            break
        if used + len(ch) > budget:
            continue
        picked.append((s, ver, ch))
        used += len(ch)
    return picked, used


# ---------------------------------------------------------------- QA
def sha256_file(path):
    with open(path, 'rb') as f:
        return hashlib.sha256(f.read()).hexdigest()


def cmd_make_qa(rows, out='qa_candidates.json'):
    """各引継書の「確定」を含む行からQA候補の雛形を作る(人手で質問文を書き直す)。"""
    items = []
    for i, (hid, ver, content, created, fn) in enumerate(rows):
        if not content:
            continue
        for line in content.split('\n'):
            if '確定' in line and 10 < len(line.strip()) < 200:
                items.append({
                    'q': 'TODO: この事実を自分の言葉で問いに直す → ' + line.strip()[:100],
                    'target_ver': str(core.normalize_existing_version(ver) or ver),
                    'source_line': line.strip()[:150],
                })
                break  # 1版につき1候補
    with open(out, 'w', encoding='utf-8') as f:
        json.dump(items, f, ensure_ascii=False, indent=1)
    print(f"QA候補 {len(items)}件 → {out}")
    print("重要: q を必ず自分の言葉に書き直すこと(原文のままだと検索が自明になり判定が甘くなる)。")
    print("完成したら --freeze で凍結してから --run。")


def cmd_freeze(qa_path):
    h = sha256_file(qa_path)
    with open(qa_path + '.sha256', 'w') as f:
        f.write(h + '\n')
    with open(qa_path, encoding='utf-8') as f:
        n = len(json.load(f))
    print(f"凍結完了: {qa_path} ({n}問)")
    print(f"SHA-256: {h}")
    print("以後このファイルを変更すると --run は拒否されます。")


# ---------------------------------------------------------------- 評価
def bin_of(idx):
    for name, lo, hi in BINS:
        if lo <= idx <= hi:
            return name
    return 'old'


def cmd_run(rows, qa_path, top_k, budget):
    sha_path = qa_path + '.sha256'
    if not os.path.exists(sha_path):
        sys.exit("エラー: 凍結されていません。先に --freeze してください。")
    frozen = open(sha_path).read().strip()
    now = sha256_file(qa_path)
    if frozen != now:
        sys.exit("エラー: QAファイルが凍結後に変更されています。実行拒否。\n"
                 f"  凍結時: {frozen}\n  現在:   {now}")
    with open(qa_path, encoding='utf-8') as f:
        qa = json.load(f)
    qa = [x for x in qa if not str(x.get('q', '')).startswith('TODO')]
    if not qa:
        sys.exit("エラー: 有効な質問がありません (TODO のままの項目は除外されます)。")

    # 版番号 → 新しい順インデックス
    ver_to_idx = {}
    for i, (hid, ver, content, created, fn) in enumerate(rows):
        nv = core.normalize_existing_version(ver)
        if nv:
            ver_to_idx[nv] = i

    # 質問をビンに割当
    for item in qa:
        tv = str(item['target_ver']).lstrip('vV')
        if tv not in ver_to_idx:
            sys.exit(f"エラー: 標的版 v{tv} がDBにありません: {item['q'][:40]}")
        item['idx'] = ver_to_idx[tv]
        item['bin'] = bin_of(item['idx'])
        item['tv'] = tv

    bin_counts = {b: sum(1 for x in qa if x['bin'] == b) for b, _, _ in BINS}
    print(f"QA: {len(qa)}問  ビン内訳: {bin_counts}")
    print(f"凍結ハッシュ: {frozen[:16]}...")
    print(f"事前登録: MARGIN={MARGIN}, MIN_PER_BIN={MIN_PER_BIN}, "
          f"RANDOM_SEEDS={RANDOM_SEEDS}, budget={budget}, top_k={top_k}\n")

    n = len(rows)
    st = strategies(n)
    results = {}   # name -> {'hits': {bin: [0/1,...]}, 'chars': [..], 'kept': k}
    for name, idxs in st.items():
        subset = [rows[i] for i in idxs]
        hits = {b: [] for b, _, _ in BINS}
        chars = []
        for item in qa:
            picked, used = retrieve_meta(subset, item['q'], top_k, budget)
            hit = any(v == item['tv'] for _, v, _ in picked)
            hits[item['bin']].append(1 if hit else 0)
            chars.append(used)
        results[name] = {'hits': hits, 'chars': chars, 'kept': len(idxs)}

    # RANDOM を平均に集約
    rnd_names = [k for k in results if k.startswith('RANDOM')]

    def overall(name):
        allh = [h for b, _, _ in BINS for h in results[name]['hits'][b]]
        return sum(allh) / len(allh)

    def binrate(name, b):
        hs = results[name]['hits'][b]
        return sum(hs) / len(hs) if hs else float('nan')

    # ---- 表 ----
    order = ['FULL', 'FIB', 'DOUBLE', 'UNIFORM', 'LATEST']
    print(f"{'strategy':>9} | {'kept':>4} | "
          + ' | '.join(f'{b:>7}' for b, _, _ in BINS)
          + f" | {'total':>6} | {'平均文字':>7}")
    print('-' * 78)
    for name in order:
        r = results[name]
        cells = ' | '.join(f'{binrate(name, b)*100:6.1f}%' for b, _, _ in BINS)
        print(f"{name:>9} | {r['kept']:>4} | {cells} | "
              f"{overall(name)*100:5.1f}% | {sum(r['chars'])//len(r['chars']):>7}")
    if rnd_names:
        kept = results[rnd_names[0]]['kept']
        cells = []
        for b, _, _ in BINS:
            v = sum(binrate(x, b) for x in rnd_names) / len(rnd_names)
            cells.append(f'{v*100:6.1f}%')
        rnd_total = sum(overall(x) for x in rnd_names) / len(rnd_names)
        rnd_chars = sum(sum(results[x]['chars']) for x in rnd_names) // (
            len(rnd_names) * len(qa))
        print(f"{'RANDOM*':>9} | {kept:>4} | " + ' | '.join(cells)
              + f" | {rnd_total*100:5.1f}% | {rnd_chars:>7}")
        print(f"(RANDOM* は {len(rnd_names)} シードの平均)")

    # ---- 判定 (ヒット数の整数比較 — 浮動小数点の境界誤差を排除) ----
    print('\n' + '=' * 78)
    small_bins = [b for b, c in bin_counts.items() if c < MIN_PER_BIN]

    def hitcount(name):
        return sum(h for b, _, _ in BINS for h in results[name]['hits'][b])

    nq = len(qa)
    fib_h, dbl_h = hitcount('FIB'), hitcount('DOUBLE')
    fib_t, dbl_t = fib_h / nq, dbl_h / nq
    rnd_h = (sum(hitcount(x) for x in rnd_names) / len(rnd_names)
             if rnd_names else 0.0)
    rnd_t = rnd_h / nq
    diff_h = fib_h - dbl_h                # 整数
    need_h = MARGIN * nq                  # FORWARD に必要なヒット差
    print(f"FIB={fib_h}/{nq} ({fib_t*100:.1f}%)  DOUBLE={dbl_h}/{nq} "
          f"({dbl_t*100:.1f}%)  差={diff_h}問 ({diff_h/nq*100:+.1f}pt)  "
          f"必要差={need_h:.1f}問  RANDOM平均={rnd_h:.1f}/{nq}")
    boundary = abs(diff_h - need_h) < 1   # 1問差で反転する境界
    if small_bins:
        verdict = 'BLIND'
        reason = (f"ビン {small_bins} の質問数が {MIN_PER_BIN} 未満。"
                  "判定材料不足につき結論保留。")
    elif fib_h <= rnd_h + 1e-9:
        verdict = 'ABSTAIN'
        reason = "FIBが同予算RANDOM以下。保持構造が偶然を上回らない。"
    elif diff_h + 1e-9 >= need_h:
        verdict = 'FORWARD'
        reason = (f"FIB−DOUBLE = {diff_h}問 ≥ 必要差{need_h:.1f}問。"
                  "事前登録マージンを満たし、φ固有の優位を支持。")
        if boundary:
            reason += " ただし境界ちょうど: 1問の反転で判定が変わる最弱のFORWARD。"
    else:
        verdict = 'ABSTAIN'
        reason = (f"FIB−DOUBLE = {diff_h}問 < 必要差{need_h:.1f}問。"
                  "写像M(φ固有の優位)は棄却。対数成長という核は本判定の対象外で無傷。")
        if boundary:
            reason += " ただし境界近傍: 1問の反転で判定が変わりうる。"
    print(f"\n判定: {verdict}")
    print(f"理由: {reason}")
    print("注記: ABSTAINは写像Mの棄却であり、フィボナッチ忘却システム自体の"
          "棄却ではない(公理は聖域、写像は判定台)。")


# ---------------------------------------------------------------- main
def main():
    p = argparse.ArgumentParser(
        description='フィボナッチ忘却のハーネス検証',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--make-qa', action='store_true', help='QA候補の雛形を生成')
    p.add_argument('--freeze', metavar='QA.json', help='QAセットをハッシュ凍結')
    p.add_argument('--run', metavar='QA.json', help='評価実行(凍結必須)')
    p.add_argument('--out', default='qa_candidates.json')
    p.add_argument('--budget', type=int, default=6000)
    p.add_argument('--top-k', type=int, default=8)
    args = p.parse_args()

    if args.freeze:
        cmd_freeze(args.freeze)
        return

    conn = psycopg2.connect(**core.DB)
    cur = conn.cursor()
    try:
        m = core.detect_columns(cur)
        rows = core.fetch_all(cur, m)
        if args.make_qa:
            cmd_make_qa(rows, args.out)
        elif args.run:
            cmd_run(rows, args.run, args.top_k, args.budget)
        else:
            p.print_help()
    finally:
        conn.close()


if __name__ == '__main__':
    main()
