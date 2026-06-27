"""log_first_real.py
初回実観測ログを「そのまま」固定するためのロガー。評価器ではない。
"""

import sys
import re
import io
from a_prep import run_bundle


def _safe_name(target):
    """target をファイル名に使える形へ（値の改変はしない・命名のみ）。"""
    s = re.sub(r"[^0-9A-Za-z_.-]+", "_", str(target)).strip("_")
    return s or "unnamed"


class _Tee:
    """画面出力を保ちつつ buffer にも複写するだけ。内容は一切加工しない。"""

    def __init__(self, *streams):
        self._streams = streams

    def write(self, data):
        for s in self._streams:
            s.write(data)
        return len(data)

    def flush(self):
        for s in self._streams:
            s.flush()


def log_first_real(signals, target="(unnamed target)", n_surr=200):
    """run_bundle を一度だけ通し、生ログを first_real_log に複写する。"""
    path = f"first_real_log_{_safe_name(target)}.txt"

    buf = io.StringIO()
    old = sys.stdout

    try:
        sys.stdout = _Tee(old, buf)
        result = run_bundle(signals, target=target, n_surr=n_surr)
    finally:
        sys.stdout = old

    captured = buf.getvalue()

    with open(path, "w", encoding="utf-8") as f:
        f.write(captured)
        f.write("\n# ---- raw return dict（run_bundle の返り値・改変なし）----\n")
        f.write(f"# target = {result['target']}\n")
        f.write(f"# K = {result['K']}\n")
        f.write(f"# bundle = {result['bundle']}\n")
        f.write(f"# all_blind = {result['all_blind']}\n")

        for i, r in enumerate(result["rows"]):
            f.write(
                f"# obs{i}: verdict={r['verdict']} "
                f"detected_s={r['detected_s']} "
                f"real_d={r['real_d']} "
                f"gap99={r['gap99']} "
                f"guard_reason={r['guard_reason']!r}\n"
            )

        f.write(f"# prior = {result['prior']}\n")

    print(f"\n[logger] 生ログを {path} に固定した（判定なし・解釈なし）。")
    return result


if __name__ == "__main__":
    print("log_first_real.py はロガー。実観測 bundle を渡して呼ぶ:")
    print("  from log_first_real import log_first_real")
    print('  log_first_real([obs0, obs1, obs2], target="...")')
    print("評価器ではない。生ログを固定するだけ。")
