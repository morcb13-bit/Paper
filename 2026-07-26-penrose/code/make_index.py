#  置いてあるものの目次を作る
#      使い方： python3 make_index.py <ディレクトリ> [出力.md]
#      .py  … 冒頭の連続したコメント行（監督の流儀なら「検定名／対象を一行で」がここにある）
#      .svg … <title> ／ 無ければ最初の <text> ／ viewBox
#      .html… <title> ／ 無ければ最初の <h1>
#      .png … 幅×高さ（IHDR を直接読む。外部の library を使わない）
#      .md  … 最初の見出し行
#      その他 … 大きさだけ
#      中身は書き換えない。読むだけ。

import os, re, sys, zlib, datetime

def head_comment(p, n=3):
    out = []
    try:
        with open(p, encoding='utf-8', errors='replace') as f:
            for line in f:
                s = line.rstrip('\n')
                if s.startswith('#'):
                    t = s.lstrip('#').strip()
                    if t:
                        out.append(t)
                    if len(out) >= n:
                        break
                elif out or s.strip():
                    break
    except OSError:
        pass
    return out

def tag(p, name):
    try:
        s = open(p, encoding='utf-8', errors='replace').read(20000)
    except OSError:
        return None
    m = re.search(r'<%s[^>]*>(.*?)</%s>' % (name, name), s, re.S | re.I)
    return re.sub(r'\s+', ' ', m.group(1)).strip() if m else None

def viewbox(p):
    try:
        s = open(p, encoding='utf-8', errors='replace').read(4000)
    except OSError:
        return None
    m = re.search(r'viewBox="([^"]+)"', s)
    return m.group(1) if m else None

def png_size(p):
    try:
        with open(p, 'rb') as f:
            d = f.read(33)
        if d[:8] != b'\x89PNG\r\n\x1a\n' or d[12:16] != b'IHDR':
            return None
        w = int.from_bytes(d[16:20], 'big'); h = int.from_bytes(d[20:24], 'big')
        return "%d×%d" % (w, h)
    except OSError:
        return None

def md_head(p):
    try:
        for line in open(p, encoding='utf-8', errors='replace'):
            if line.startswith('#'):
                return line.lstrip('#').strip()
    except OSError:
        pass
    return None

def describe(p):
    e = os.path.splitext(p)[1].lower()
    if e == '.py':
        h = head_comment(p)
        return " ／ ".join(h) if h else "（冒頭にコメントなし）"
    if e == '.svg':
        return tag(p, 'title') or (tag(p, 'text') and "文字：" + tag(p, 'text')) or \
               (viewbox(p) and "viewBox " + viewbox(p)) or ""
    if e in ('.html', '.htm'):
        return tag(p, 'title') or tag(p, 'h1') or ""
    if e == '.png':
        return png_size(p) or ""
    if e in ('.md', '.markdown'):
        return md_head(p) or ""
    return ""

def main(root, out=None):
    rows = []
    for dp, dns, fns in os.walk(root):
        dns[:] = [d for d in dns if d not in ('.git', '__pycache__', 'node_modules')]
        for fn in sorted(fns):
            p = os.path.join(dp, fn)
            try:
                st = os.stat(p)
            except OSError:
                continue
            rows.append((os.path.relpath(p, root), st.st_size,
                         datetime.date.fromtimestamp(st.st_mtime).isoformat(),
                         describe(p).replace('|', '/')))
    lines = ["# 置いてあるもの（%d件・%s 時点）" % (len(rows), datetime.date.today().isoformat()), ""]
    cur = None
    for rel, size, mt, desc in sorted(rows):
        d = os.path.dirname(rel) or "."
        if d != cur:
            cur = d
            lines += ["", "## %s" % d, "", "| ファイル | 大きさ | 更新 | 中身 |", "|---|---|---|---|"]
        lines.append("| `%s` | %s | %s | %s |"
                     % (os.path.basename(rel),
                        ("%.0fK" % (size / 1024)) if size >= 1024 else "%dB" % size,
                        mt, desc[:200]))
    text = "\n".join(lines) + "\n"
    if out:
        open(out, 'w', encoding='utf-8').write(text)
        print("→ %s（%d件）" % (out, len(rows)))
    else:
        print(text)

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
