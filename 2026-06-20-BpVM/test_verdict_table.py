"""
test_verdict_table.py  —  v160 凍結ガード: _verdict の 2x2 全象限を pin。
ε=EPS_SEP は排他（gap99==ε は「分離なし」）であることを境界で固定する。
実機 engine が触らない SIGNAL_WEAK 隅も含めて 4 象限すべてを回帰固定。
"""
from engine_reversal_eval import _verdict, EPS_SEP, REAL_OK

cases = [
    # (real_d, gap99, expected_code, 説明)
    (1.00, 0.46,  "ADVANCE",              "real健在 × 分離あり"),
    (1.00, 0.00,  "BLIND_SUSPECT_REJECT", "real健在 × 分離なし（信号盲）"),
    (0.50, 0.30,  "SIGNAL_WEAK",          "real低下 × 分離あり（微弱だが本物）"),
    (0.30, 0.00,  "BLIND_DROP",           "real低下 × 分離なし（盲目）"),
    # 境界: ε は排他、REAL_OK は包含
    (REAL_OK,        EPS_SEP,        "BLIND_SUSPECT_REJECT", "gap99==ε は分離なし"),
    (REAL_OK,        EPS_SEP + 1e-9, "ADVANCE",              "gap99>ε で分離"),
    (REAL_OK - 1e-9, 0.30,           "SIGNAL_WEAK",          "real_d<REAL_OK"),
    (REAL_OK,        0.30,           "ADVANCE",              "real_d==REAL_OK は健在"),
]

print("=" * 72)
print(f"verdict 2x2 table  (EPS_SEP={EPS_SEP}, REAL_OK={REAL_OK})")
print("=" * 72)
fails = 0
for real_d, gap99, want, desc in cases:
    code, label = _verdict(real_d, gap99)
    ok = code == want
    fails += (not ok)
    mark = "✓" if ok else "✗"
    print(f"  [{mark}] real_d={real_d:<6} gap99={gap99:<10} → {code:<22} ({desc})")
    if not ok:
        print(f"        期待 {want} / 実 {code}")
print("-" * 72)
print("ALL PASS" if fails == 0 else f"{fails} FAIL")
