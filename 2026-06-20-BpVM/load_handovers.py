#!/usr/bin/env python3
"""
load_handovers.py
溜まった引継書 .md を handover テーブルへ一括投入する(冪等)。

スキーマ前提(検索スクリプトの fetch から): handover(id, version, content, created_at, filename)
  - 既存 version は二重登録しない。--reload 指定時のみ content/filename/created_at を更新。
  - created_at はファイルの mtime を使う(一括投入で created_at が潰れるのを防ぐ)。
  - version 列が int/text どちらでも version::text=%s で照合。
  - handover テーブルが無ければ作る。実カラムを information_schema で確認して動的に組む。

使い方:
  python3 load_handovers.py --dir /home/moroc/handovers            # 投入
  python3 load_handovers.py --dir . --dry                          # 走査だけ(DBに触らない)
  python3 load_handovers.py --dir . --reload                       # 既存も content 上書き
  python3 load_handovers.py --dir . --pattern 'B13_引継書_v*.md'   # パターン変更
"""
import argparse, re, os, glob, sys
from datetime import datetime, timezone

DB = dict(host='localhost', dbname='b13db', user='moroc', password='b13research2026')
VER_RE = re.compile(r'v(\d+)', re.I)


def parse_version(fn):
    m = VER_RE.search(fn)
    return int(m.group(1)) if m else None


def scan(directory, pattern):
    out = []
    for p in sorted(glob.glob(os.path.join(directory, pattern))):
        if not os.path.isfile(p):
            continue
        fn = os.path.basename(p)
        ver = parse_version(fn)
        with open(p, encoding='utf-8') as f:
            content = f.read()
        mtime = datetime.fromtimestamp(os.path.getmtime(p), tz=timezone.utc)
        out.append(dict(path=p, filename=fn, version=ver, content=content, mtime=mtime))
    # version 昇順(None は末尾)で安定化
    out.sort(key=lambda r: (r['version'] is None, r['version'] or 0, r['filename']))
    return out


def detect_columns(cur):
    cur.execute("""SELECT column_name FROM information_schema.columns
                   WHERE table_name='handover'""")
    return {r[0] for r in cur.fetchall()}


def ensure_table(cur):
    cols = detect_columns(cur)
    if not cols:
        cur.execute("""
            CREATE TABLE handover (
                id          SERIAL PRIMARY KEY,
                version     INTEGER,
                filename    TEXT,
                content     TEXT NOT NULL,
                created_at  TIMESTAMPTZ DEFAULT now()
            )""")
        print("handover テーブルを新規作成しました。")
        cols = detect_columns(cur)
    return cols


def exists_id(cur, cols, rec):
    """同一引継書が既にあるか。version 優先、無ければ filename で照合。"""
    if rec['version'] is not None and 'version' in cols:
        cur.execute("SELECT id FROM handover WHERE version::text=%s", (str(rec['version']),))
    elif 'filename' in cols:
        cur.execute("SELECT id FROM handover WHERE filename=%s", (rec['filename'],))
    else:
        return None
    row = cur.fetchone()
    return row[0] if row else None


def upsert(cur, cols, rec, reload_existing):
    hid = exists_id(cur, cols, rec)
    has_created = 'created_at' in cols
    has_ver = 'version' in cols
    has_fn = 'filename' in cols

    if hid is not None:
        if not reload_existing:
            return 'skip'
        sets, vals = ["content=%s"], [rec['content']]
        if has_fn:      sets.append("filename=%s");   vals.append(rec['filename'])
        if has_created: sets.append("created_at=%s"); vals.append(rec['mtime'])
        vals.append(hid)
        cur.execute(f"UPDATE handover SET {', '.join(sets)} WHERE id=%s", vals)
        return 'update'

    fields, ph, vals = ['content'], ['%s'], [rec['content']]
    if has_ver:     fields.append('version');    ph.append('%s'); vals.append(rec['version'])
    if has_fn:      fields.append('filename');   ph.append('%s'); vals.append(rec['filename'])
    if has_created: fields.append('created_at'); ph.append('%s'); vals.append(rec['mtime'])
    cur.execute(f"INSERT INTO handover ({', '.join(fields)}) VALUES ({', '.join(ph)})", vals)
    return 'insert'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', default='.', help='引継書 .md のあるディレクトリ')
    ap.add_argument('--pattern', default='B13_引継書_v*.md', help='ファイル名パターン')
    ap.add_argument('--reload', action='store_true', help='既存 version も content を上書き')
    ap.add_argument('--dry', action='store_true', help='走査だけ。DBに触らない')
    args = ap.parse_args()

    recs = scan(args.dir, args.pattern)
    if not recs:
        print(f"該当ファイルなし: {os.path.join(args.dir, args.pattern)}")
        sys.exit(1)

    print(f"検出 {len(recs)} 件:")
    for r in recs:
        v = f"v{r['version']}" if r['version'] is not None else "v?"
        print(f"  {v:>6} | {len(r['content']):>7}字 | {r['mtime']:%Y-%m-%d %H:%M} | {r['filename']}")

    if args.dry:
        novers = [r['filename'] for r in recs if r['version'] is None]
        if novers:
            print("\n[警告] version を抽出できないファイル(filename で照合・投入されます):")
            for fn in novers: print(f"  - {fn}")
        print("\n--dry のため DB には触れていません。")
        return

    import psycopg2  # 接続時のみ必要(--dry では不要)
    conn = psycopg2.connect(**DB)
    cur = conn.cursor()
    try:
        cols = ensure_table(cur)
        tally = dict(insert=0, update=0, skip=0)
        for r in recs:
            tally[upsert(cur, cols, r, args.reload)] += 1
        conn.commit()
        print(f"\n結果: 新規 {tally['insert']} / 更新 {tally['update']} / スキップ既存 {tally['skip']}")
        # 検証
        cur.execute("SELECT count(*) FROM handover")
        total = cur.fetchone()[0]
        order = "version DESC" if 'version' in cols else "created_at DESC"
        cur.execute(f"SELECT version, filename FROM handover ORDER BY {order} LIMIT 1"
                    if 'version' in cols else
                    f"SELECT filename FROM handover ORDER BY {order} LIMIT 1")
        latest = cur.fetchone()
        print(f"テーブル総数: {total} / 最新: {latest}")
    finally:
        conn.close()


if __name__ == '__main__':
    main()
