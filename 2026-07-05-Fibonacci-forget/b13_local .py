#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
b13_local.py — 引継書DB × ローカルLLM (Ollama) 統合ツール

モード:
  --register FILE        新しい引継書を登録 (FILE='-' で標準入力から)
                         版番号はファイル名/本文から自動抽出、--ver で明示指定可
                         同じ版が既にあると中断 (--force で強行)
  --register-dir DIR     フォルダ内の .md/.txt をまとめて登録
                         登録済みの版・版番号不明のファイルは自動スキップ
                         --date mtime (既定): created_at にファイル更新日時を使用
                         --date now: 現在時刻+版順の秒差で登録
                         --dry-run: 登録せずに計画だけ表示
  --version              引継書一覧を表示
  --search REGEX         正規表現で全文検索(従来機能)
  --preset phi-sqrt2     φ↔√2 対立軸のプリセット検索(近接ペア強調)
  --ask "質問"           DBから関連箇所を検索し、ローカルLLMに回答させる(簡易RAG)
  --next                 最新引継書から「次にやるべきこと」を答えさせる
  --summary [N]          最新N件(既定1件)を要約させる
  --chat                 対話モード(毎ターンDBを再検索して回答)
  --raw                  --ask/--chat でLLMを呼ばず検索結果のみ表示
                         (Claude や Chappy に貼り付ける用)

フィボナッチ忘却:
  --fib                  想起対象を新しい方から 0,1,2,3,5,8,13,21,... 番目の
                         引継書に限定する(密→疎の忘却則)
  --decay                --fib 使用時、古い階層(tier)ほど関連度を 1/φ ずつ減衰
  --version --fib        どの引継書が KEEP / FORGET かを一覧表示

共通オプション:
  --model NAME           Ollamaモデル (既定: lucas2024/llama-3-elyza-jp-8b:q5_k_m)
  --url URL              OllamaのURL (既定: http://localhost:11434)
  --budget N             LLMに渡す参照資料の最大文字数 (既定: 6000)
  --top-k N              使用する最大チャンク数 (既定: 8)
  --num-ctx N            Ollamaのコンテキスト長 (既定: 8192)
  --context N            検索モードでの前後表示行数 (既定: 2)
  --show-context         --ask 時にLLMへ渡した参照資料も表示

環境変数で上書き可: B13_DB_HOST, B13_DB_NAME, B13_DB_USER, B13_DB_PASSWORD,
                    B13_OLLAMA_URL, B13_MODEL
"""
import argparse
import json
import os
import re
import sys

import psycopg2
import requests

# ---------------------------------------------------------------- 設定
DB = dict(
    host=os.environ.get('B13_DB_HOST', 'localhost'),
    dbname=os.environ.get('B13_DB_NAME', 'b13db'),
    user=os.environ.get('B13_DB_USER', 'moroc'),
    password=os.environ.get('B13_DB_PASSWORD', 'b13research2026'),
)
DEFAULT_URL = os.environ.get('B13_OLLAMA_URL', 'http://localhost:11434')
DEFAULT_MODEL = os.environ.get('B13_MODEL', 'lucas2024/llama-3-elyza-jp-8b:q5_k_m')

# φ↔√2 プリセット
PHI = [r'φ', r'黄金比', r'黄金角', r'フィボナッチ', r'Fibonacci',
       r'Cassini', r'ビネ', r'Binet', r'(?<![A-Za-z])phi(?![A-Za-z])']
SQRT2 = [r'√2', r'ルート2', r'sqrt\(?2\)?', r'sqrt2', r'ペル',
         r'(?<![A-Za-z])Pell(?![A-Za-z])', r'魔法数', r'マジック\s*ナンバー',
         r'magic\s*number', r'閉殻', r'1\s*\+\s*√2']
OPP = [r'対立', r'相容れない', r'相反', r'対抗',
       r'(?:vs|VS|対)\s*(?:√2|φ|黄金|ペル|フィボナッチ)',
       r'(?:√2|φ|黄金|ペル|フィボナッチ)\s*(?:vs|VS|対)']

# ---------------------------------------------------------------- DBアクセス
CONTENT_CANDS = ['content', 'body', 'text', 'document', 'doc']
VERSION_CANDS = ['version', 'ver', 'v', 'rev']
ID_CANDS = ['id', 'handover_id', 'pk']
DATE_CANDS = ['created_at', 'created', 'timestamp', 'ts', 'date']
FILE_CANDS = ['filename', 'file', 'name', 'title']


def detect_columns(cur):
    """information_schema から handover テーブルの実カラム名を検出する。"""
    cur.execute(
        "SELECT column_name, data_type FROM information_schema.columns "
        "WHERE table_name = 'handover'"
    )
    types = dict(cur.fetchall())
    cols = list(types)
    if not cols:
        sys.exit("エラー: handover テーブルが見つかりません。"
                 " psql で \\dt を確認してください。")

    def pick(cands):
        for c in cands:
            if c in cols:
                return c
        return None

    m = {
        'content': pick(CONTENT_CANDS),
        'version': pick(VERSION_CANDS),
        'id': pick(ID_CANDS),
        'created': pick(DATE_CANDS),
        'filename': pick(FILE_CANDS),
        'types': types,
    }
    if not m['content']:
        sys.exit(f"エラー: 本文カラムを特定できません。実カラム: {cols}\n"
                 f"CONTENT_CANDS にカラム名を追加してください。")
    return m


def fetch_all(cur, m):
    """全引継書を新しい順に取得。無いカラムは NULL で埋める。"""
    sel = [
        m['id'] or 'NULL',
        m['version'] or 'NULL',
        m['content'],
        m['created'] or 'NULL',
        m['filename'] or 'NULL',
    ]
    order = f"ORDER BY {m['created']} DESC" if m['created'] else ''
    cur.execute(f"SELECT {', '.join(sel)} FROM handover {order}")
    return cur.fetchall()


# ------------------------------------------------------- フィボナッチ忘却
GOLDEN = 0.6180339887  # 1/φ

def fib_keep_indices(n):
    """新しい方から数えて保持するインデックス列: 0,1,2,3,5,8,13,21,..."""
    keep = [0]
    a, b = 1, 2
    while a < n:
        keep.append(a)
        a, b = b, a + b
    return keep


def fibonacci_select(rows):
    """フィボナッチ忘却則を適用し (tier, row) のリストを返す。
    rows は新しい順であること。tier=0 が最新、番号が増えるほど古い階層。"""
    idx = fib_keep_indices(len(rows))
    return [(t, rows[i]) for t, i in enumerate(idx)]


def with_tiers(rows):
    """忘却なし: 全引継書を tier=0 として扱う。"""
    return [(0, r) for r in rows]


def label(ver, created, fn):
    parts = []
    if ver is not None:
        v = str(ver)
        parts.append(v if v.lower().startswith('v') else f'v{v}')
    if created is not None:
        parts.append(str(created)[:16])
    if fn:
        parts.append(str(fn))
    return ' / '.join(parts) if parts else '(不明)'


# ---------------------------------------------------------------- 登録
VER_RX = re.compile(r'[vV](\d{1,4})')


def guess_version(path, text):
    """ファイル名 → 本文冒頭の順で v数字 を探す。"""
    for src in (os.path.basename(path or ''), text[:500]):
        mt = VER_RX.search(src or '')
        if mt:
            return mt.group(1)
    return None


def insert_handover(cur, m, text, ver, fname, ts=None):
    """実カラム構成に合わせて1件INSERTする。ts=None なら NOW()。"""
    cols, vals = [m['content']], [text]
    if m['version'] and ver is not None:
        vtype = m['types'].get(m['version'], '')
        vals.append(int(ver) if 'int' in vtype else str(ver))
        cols.append(m['version'])
    if m['filename'] and fname:
        cols.append(m['filename'])
        vals.append(fname)
    ph = ['%s'] * len(vals)
    if m['created']:
        cols.append(m['created'])
        if ts is not None:
            ph.append('%s')
            vals.append(ts)
        else:
            ph.append('NOW()')
    sql = f"INSERT INTO handover ({', '.join(cols)}) VALUES ({', '.join(ph)})"
    cur.execute(sql, vals)


def cmd_register(cur, conn, m, path, ver=None, force=False):
    # 本文の読み込み
    if path == '-':
        print("標準入力から読み込みます (終了は Ctrl-D):")
        text = sys.stdin.read()
        path = None
    else:
        if not os.path.exists(path):
            sys.exit(f"エラー: ファイルが見つかりません: {path}")
        with open(path, encoding='utf-8') as f:
            text = f.read()
    if not text.strip():
        sys.exit("エラー: 本文が空です。登録を中止しました。")

    # 版番号
    ver = ver or guess_version(path, text)
    if ver is None:
        sys.exit("エラー: 版番号を特定できません。--ver 177 のように指定するか、"
                 "ファイル名か本文冒頭に v177 の形式で含めてください。")
    ver = str(ver).lstrip('vV')

    # 重複チェック
    if m['version'] and not force:
        cur.execute(
            f"SELECT COUNT(*) FROM handover WHERE {m['version']}::text = %s",
            (ver,))
        if cur.fetchone()[0] > 0:
            sys.exit(f"中断: v{ver} は既に登録されています。"
                     " 上書きではなく追加登録するなら --force を付けてください。")

    # INSERT
    insert_handover(cur, m, text, ver,
                    os.path.basename(path) if path else None)
    conn.commit()

    print(f"登録完了: v{ver} ({len(text)}文字"
          + (f", {os.path.basename(path)}" if path else ", 標準入力") + ")")
    print("この引継書がフィボナッチ忘却の tier0 (最新) になります。")
    rows = fetch_all(cur, m)
    print("\n[最新3件]")
    cmd_versions(rows[:3])


def cmd_register_dir(cur, conn, m, dirpath, date_mode='mtime', dry_run=False):
    """フォルダ内の .md/.txt をまとめて登録。登録済み版はスキップ(冪等)。"""
    import datetime
    import glob as _glob
    if not os.path.isdir(dirpath):
        sys.exit(f"エラー: フォルダが見つかりません: {dirpath}")
    files = sorted(_glob.glob(os.path.join(dirpath, '*.md'))
                   + _glob.glob(os.path.join(dirpath, '*.txt')))
    if not files:
        sys.exit(f"エラー: {dirpath} に .md / .txt ファイルがありません。")

    # 登録済みの版番号を取得
    existing = set()
    if m['version']:
        cur.execute(f"SELECT {m['version']}::text FROM handover "
                    f"WHERE {m['version']} IS NOT NULL")
        for (v,) in cur.fetchall():
            v = str(v).lstrip('vV')
            if v.isdigit():
                existing.add(str(int(v)))

    plan, skipped, seen = [], [], set()
    for f in files:
        try:
            with open(f, encoding='utf-8') as fh:
                text = fh.read()
        except UnicodeDecodeError:
            skipped.append((f, '文字コードがUTF-8でない'))
            continue
        if not text.strip():
            skipped.append((f, '本文が空'))
            continue
        ver = guess_version(f, text)
        if ver is None:
            skipped.append((f, '版番号を特定できない'))
            continue
        ver = str(int(ver))
        if ver in existing:
            skipped.append((f, f'v{ver} は登録済み'))
            continue
        if ver in seen:
            skipped.append((f, f'v{ver} がフォルダ内で重複'))
            continue
        seen.add(ver)
        ts = (datetime.datetime.fromtimestamp(os.path.getmtime(f))
              if date_mode == 'mtime' else None)
        plan.append((int(ver), f, text, ts))

    plan.sort()  # 版番号の昇順で登録(古い→新しい)
    if date_mode == 'now':
        base = datetime.datetime.now()
        plan = [(v, f, t, base + datetime.timedelta(seconds=i))
                for i, (v, f, t, _) in enumerate(plan)]

    tag = '[DRY-RUN] ' if dry_run else ''
    print(f"{tag}対象フォルダ: {dirpath}  (ファイル {len(files)}件)")
    print(f"{tag}登録予定: {len(plan)}件 / スキップ: {len(skipped)}件\n")
    if plan:
        print(f"{'ver':>6} | {'created_at に使う日時':>19} | {'len':>7} | file")
        print('-' * 80)
        for v, f, text, ts in plan:
            print(f"{'v'+str(v):>6} | {str(ts)[:19]:>19} | "
                  f"{len(text):>7} | {os.path.basename(f)}")
    if skipped:
        print("\n[スキップ]")
        for f, why in skipped:
            print(f"  {os.path.basename(f)}: {why}")
    if dry_run:
        print("\nDRY-RUN のため登録していません。問題なければ --dry-run を外して実行してください。")
        return
    if not plan:
        print("\n登録するものがありません。")
        return

    for v, f, text, ts in plan:
        insert_handover(cur, m, text, str(v), os.path.basename(f), ts)
    conn.commit()
    print(f"\n登録完了: {len(plan)}件")
    if date_mode == 'mtime':
        print("created_at にはファイル更新日時(mtime)を使いました。"
              "フィボナッチ忘却の順序はこの日時基準になります。")
    rows = fetch_all(cur, m)
    print("\n[--version --fib での現在の保持状況]")
    cmd_versions(rows, fib=True)



# ---------------------------------------------------------------- 検索(従来)
def hits(text, rxs):
    out = []
    for i, line in enumerate(text.split('\n'), 1):
        for r in rxs:
            mt = r.search(line)
            if mt:
                out.append((i, line, mt.group(0)))
                break
    return out


def cmd_versions(rows, fib=False):
    keep = {}
    if fib:
        keep = {i: t for t, i in enumerate(fib_keep_indices(len(rows)))}
        print(f"フィボナッチ忘却: {len(rows)}件中 {len(keep)}件を想起対象に保持")
    head = f"{'id':>4} | {'ver':>8} | {'created_at':>20} | {'len':>7} | filename"
    if fib:
        head += ' | fib'
    print(head)
    print('-' * 100)
    for i, (hid, ver, content, created, fn) in enumerate(rows):
        line = (f"{str(hid):>4} | {str(ver):>8} | {str(created):>20} | "
                f"{len(content or ''):>7} | {fn}")
        if fib:
            line += (f" | KEEP(tier{keep[i]})" if i in keep else ' | forget')
        print(line)


def cmd_search(rows, pattern, ctx):
    rx = re.compile(pattern, re.I)
    total = 0
    for hid, ver, content, created, fn in rows:
        if not content:
            continue
        lines = content.split('\n')
        hs = [(i, l) for i, l in enumerate(lines, 1) if rx.search(l)]
        if not hs:
            continue
        print(f"\n■ {label(ver, created, fn)} (id={hid})  ヒット: {len(hs)}件")
        for ln, _ in hs[:20]:
            lo, hi = max(0, ln - 1 - ctx), min(len(lines), ln + ctx)
            print(f"  --- L{ln} ---")
            for i in range(lo, hi):
                pre = '>>> ' if i + 1 == ln else '    '
                print(f"  {pre}L{i+1}: {lines[i].rstrip()[:130]}")
        total += len(hs)
    print(f"\n総ヒット: {total}")


def cmd_preset_phi_sqrt2(rows, ctx):
    phi_re = [re.compile(p, re.I) for p in PHI]
    sq_re = [re.compile(p, re.I) for p in SQRT2]
    op_re = [re.compile(p, re.I) for p in OPP]
    total = 0
    for hid, ver, content, created, fn in rows:
        if not content:
            continue
        ph, sq, op = hits(content, phi_re), hits(content, sq_re), hits(content, op_re)
        if not ph and not sq:
            continue
        pairs = []
        for pl, pt, pm in ph:
            for sl, st, sm in sq:
                if abs(pl - sl) <= 5:
                    pairs.append((pl, pm, sl, sm))
        if not op and not pairs:
            continue
        print(f"\n■ {label(ver, created, fn)} (id={hid})")
        print(f"  φ系: {len(ph)}, √2系: {len(sq)}, 対立語: {len(op)}, 近接ペア: {len(pairs)}")
        lines = content.split('\n')
        for ln, t, mtext in op[:10]:
            lo, hi = max(0, ln - 1 - ctx), min(len(lines), ln + ctx)
            print(f"  --- 対立語 L{ln} ('{mtext}') ---")
            for i in range(lo, hi):
                pre = '>>> ' if i + 1 == ln else '    '
                print(f"  {pre}L{i+1}: {lines[i].rstrip()[:130]}")
        for pl, pm, sl, sm in pairs[:5]:
            lo, hi = max(0, min(pl, sl) - 1 - ctx), min(len(lines), max(pl, sl) + ctx)
            print(f"  --- 近接 L{pl}('{pm}') × L{sl}('{sm}') ---")
            for i in range(lo, hi):
                mk = 'φ> ' if i + 1 == pl else '√> ' if i + 1 == sl else '   '
                print(f"  {mk}L{i+1}: {lines[i].rstrip()[:130]}")
        total += len(ph) + len(sq) + len(op)
    print(f"\n総ヒット: {total}")


# ---------------------------------------------------------------- RAG部分
HEAD_RX = re.compile(r'^\s*(#{1,6}\s|■|【|---|\*\*|={3,})')
PUNCT_RX = re.compile(r'[\s、。,．.!?！?()()「」『』:;:;・\[\]{}]')
STRONG_RX = re.compile(
    r'[A-Za-z][A-Za-z0-9_]{1,}'      # ASCII語
    r'|[ァ-ヶー]{2,}'                 # カタカナ語
    r'|[一-龥]{2,}'                   # 漢字語
    r'|φ|√2|√5|π'                    # 数学記号
)
STOP_TERMS = {'について', 'こと', 'もの', 'ため', 'よう', 'それ', 'これ',
              'どこ', 'なに', '何', 'ですか', 'ください', '教えて', '引継',
              '引継書', 'とは'}


def split_chunks(content, min_len=200, max_len=900):
    """引継書を見出し・空行を手がかりに数百字のチャンクへ分割する。"""
    chunks, buf = [], []

    def flush():
        t = '\n'.join(buf).strip()
        if t:
            chunks.append(t)
        buf.clear()

    size = 0
    for line in content.split('\n'):
        boundary = HEAD_RX.match(line) or (not line.strip())
        if boundary and size >= min_len:
            flush()
            size = 0
        buf.append(line)
        size += len(line) + 1
        if size >= max_len:
            flush()
            size = 0
    flush()
    return chunks


def query_features(question):
    strong = [t for t in STRONG_RX.findall(question) if t not in STOP_TERMS]
    plain = PUNCT_RX.sub('', question)
    bigrams = {plain[i:i + 2] for i in range(len(plain) - 1)}
    return strong, bigrams


def score_chunk(chunk, strong, bigrams):
    s = 0
    low = chunk.lower()
    for t in strong:
        n = low.count(t.lower())
        s += 3 * min(n, 5)
    s += sum(1 for b in bigrams if b in chunk)
    return s


def retrieve(tiered, question, top_k, budget, decay=False):
    """(tier, row) 群をチャンク化し、質問との関連度上位を予算内で選ぶ。
    decay=True なら tier ごとに 1/φ を掛けて古い記憶を減衰させる。"""
    strong, bigrams = query_features(question)
    scored = []
    for tier, (hid, ver, content, created, fn) in tiered:
        if not content:
            continue
        w = (GOLDEN ** tier) if decay else 1.0
        lab = label(ver, created, fn)
        if tier > 0:
            lab += f' | tier{tier}'
        for ch in split_chunks(content):
            s = score_chunk(ch, strong, bigrams) * w
            if s > 0:
                scored.append((s, lab, ch))
    scored.sort(key=lambda x: -x[0])
    picked, used = [], 0
    for s, lab, ch in scored:
        if len(picked) >= top_k:
            break
        if used + len(ch) > budget:
            continue
        picked.append((s, lab, ch))
        used += len(ch)
    return picked, strong


def build_context(picked):
    blocks = []
    for s, lab, ch in picked:
        blocks.append(f"【出典: {lab} | 関連度 {s:.1f}】\n{ch}")
    return '\n\n----\n\n'.join(blocks)


# ---------------------------------------------------------------- Ollama
def ollama_generate(url, model, prompt, num_ctx):
    try:
        resp = requests.post(
            f'{url}/api/generate',
            json={
                'model': model,
                'prompt': prompt,
                'stream': False,
                'options': {'num_ctx': num_ctx, 'temperature': 0.2},
            },
            timeout=600,
        )
    except requests.ConnectionError:
        sys.exit(f"エラー: Ollama に接続できません ({url})。"
                 " `ollama serve` が動いているか確認してください。")
    if resp.status_code == 404:
        try:
            tags = requests.get(f'{url}/api/tags', timeout=10).json()
            names = [m['name'] for m in tags.get('models', [])]
            sys.exit(f"エラー: モデル '{model}' が見つかりません。\n"
                     f"利用可能: {', '.join(names) or '(なし)'}")
        except Exception:
            sys.exit(f"エラー: モデル '{model}' が見つかりません (404)。")
    resp.raise_for_status()
    return resp.json().get('response', '').strip()


ASK_TEMPLATE = """あなたはB13理論の研究引継書データベースを参照して回答するアシスタントです。
以下の【参照資料】は引継書DBから質問に関連する箇所を抜粋したものです。
資料に書かれている内容だけに基づいて、日本語で簡潔に回答してください。
資料にない事柄は推測せず「資料内に見つかりませんでした」と答えてください。
回答の最後に、根拠にした出典(例: v109)を列挙してください。

【参照資料】
{context}

【質問】
{question}

【回答】"""


def cmd_ask(tiered, question, args, history=None):
    picked, strong = retrieve(tiered, question, args.top_k, args.budget,
                              decay=args.decay)
    if not picked:
        print("関連する箇所がDB内に見つかりませんでした。"
              f"(抽出キーワード: {strong})")
        return
    context = build_context(picked)
    if args.raw or args.show_context:
        print("=" * 70)
        print(f"抽出キーワード: {strong}")
        print(f"使用チャンク: {len(picked)}件 / {sum(len(c) for _, _, c in picked)}文字")
        print("=" * 70)
        print(context)
        print("=" * 70)
        if args.raw:
            return
    q = question
    if history:
        hist = '\n'.join(f"Q: {h[0]}\nA: {h[1][:200]}" for h in history[-2:])
        q = f"(直前のやり取り)\n{hist}\n\n(今回の質問)\n{question}"
    prompt = ASK_TEMPLATE.format(context=context, question=q)
    answer = ollama_generate(args.url, args.model, prompt, args.num_ctx)
    print(answer)
    return answer


def middle_truncate(text, budget):
    if len(text) <= budget:
        return text
    head = budget // 3
    tail = budget - head
    return text[:head] + '\n...(中略)...\n' + text[-tail:]


def cmd_next(rows, args):
    for hid, ver, content, created, fn in rows:
        if content:
            body = middle_truncate(content, args.budget)
            prompt = (f"以下はB13研究の最新引継書({label(ver, created, fn)})です。\n"
                      f"「次にやるべきこと」を箇条書きで3項目以内、日本語で答えてください。\n\n{body}")
            print(f"[最新引継書: {label(ver, created, fn)}]\n")
            print(ollama_generate(args.url, args.model, prompt, args.num_ctx))
            return
    print("引継書が見つかりません。")


def cmd_summary(tiered, n, args):
    targets = [r for _, r in tiered if r[2]][:n]
    if not targets:
        print("引継書が見つかりません。")
        return
    per = max(800, args.budget // max(len(targets), 1))
    blocks = [f"【{label(v, c, f)}】\n{middle_truncate(t, per)}"
              for _, v, t, c, f in targets]
    prompt = ("以下のB13研究引継書を、確定事項・未解決点・次の作業の3観点で"
              "日本語で要約してください。\n\n" + '\n\n----\n\n'.join(blocks))
    print(ollama_generate(args.url, args.model, prompt, args.num_ctx))


def cmd_chat(tiered, args):
    print(f"対話モード (model: {args.model})")
    if args.fib:
        print(f"  フィボナッチ忘却: ON ({len(tiered)}件を想起対象)"
              + (' + 黄金比減衰' if args.decay else ''))
    print("  :q で終了 / :raw で検索結果のみ表示を切替 / 空行は無視")
    history = []
    while True:
        try:
            q = input('\nB13> ').strip()
        except (EOFError, KeyboardInterrupt):
            print()
            break
        if not q:
            continue
        if q in (':q', ':quit', ':exit'):
            break
        if q == ':raw':
            args.raw = not args.raw
            print(f"raw表示: {'ON' if args.raw else 'OFF'}")
            continue
        a = cmd_ask(tiered, q, args, history=history)
        if a:
            history.append((q, a))


# ---------------------------------------------------------------- main
def main():
    p = argparse.ArgumentParser(
        description='引継書DB × ローカルLLM 統合ツール',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--register', metavar='FILE',
                   help="新しい引継書を登録 ('-' で標準入力)")
    p.add_argument('--register-dir', metavar='DIR',
                   help='フォルダ内の .md/.txt をまとめて登録')
    p.add_argument('--date', choices=['mtime', 'now'], default='mtime',
                   help='一括登録時の created_at (既定: mtime)')
    p.add_argument('--dry-run', action='store_true',
                   help='一括登録の計画だけ表示して登録しない')
    p.add_argument('--ver', metavar='N', help='登録時の版番号を明示指定')
    p.add_argument('--force', action='store_true',
                   help='同じ版番号があっても追加登録')
    p.add_argument('--version', action='store_true', help='引継書一覧')
    p.add_argument('--search', metavar='REGEX', help='正規表現で全文検索')
    p.add_argument('--preset', choices=['phi-sqrt2'], help='プリセット検索')
    p.add_argument('--ask', metavar='QUESTION', help='DB検索+LLM回答 (RAG)')
    p.add_argument('--next', action='store_true', help='次にやることを回答')
    p.add_argument('--summary', nargs='?', const=1, type=int, metavar='N',
                   help='最新N件を要約 (既定1)')
    p.add_argument('--chat', action='store_true', help='対話モード')
    p.add_argument('--raw', action='store_true', help='LLMを呼ばず検索結果のみ')
    p.add_argument('--show-context', action='store_true',
                   help='LLMへ渡した参照資料も表示')
    p.add_argument('--fib', action='store_true',
                   help='フィボナッチ忘却 (0,1,2,3,5,8,... 番目のみ想起)')
    p.add_argument('--decay', action='store_true',
                   help='--fib 時に古い階層を 1/φ ずつ減衰')
    p.add_argument('--model', default=DEFAULT_MODEL)
    p.add_argument('--url', default=DEFAULT_URL)
    p.add_argument('--budget', type=int, default=6000)
    p.add_argument('--top-k', type=int, default=8)
    p.add_argument('--num-ctx', type=int, default=8192)
    p.add_argument('--context', type=int, default=2, help='検索時の前後表示行数')
    args = p.parse_args()

    conn = psycopg2.connect(**DB)
    cur = conn.cursor()
    try:
        m = detect_columns(cur)
        if args.register_dir:
            cmd_register_dir(cur, conn, m, args.register_dir,
                             args.date, args.dry_run)
            return
        if args.register:
            cmd_register(cur, conn, m, args.register, args.ver, args.force)
            return
        rows = fetch_all(cur, m)
        tiered = fibonacci_select(rows) if args.fib else with_tiers(rows)
        if args.decay and not args.fib:
            print("注意: --decay は --fib と併用してください(単独では無効)。")
        if args.version:
            cmd_versions(rows, fib=args.fib)
        elif args.search:
            cmd_search(rows, args.search, args.context)
        elif args.preset == 'phi-sqrt2':
            cmd_preset_phi_sqrt2(rows, args.context)
        elif args.ask:
            cmd_ask(tiered, args.ask, args)
        elif args.next:
            cmd_next(rows, args)
        elif args.summary is not None:
            cmd_summary(tiered, args.summary, args)
        elif args.chat:
            cmd_chat(tiered, args)
        else:
            p.print_help()
    finally:
        conn.close()


if __name__ == '__main__':
    main()
