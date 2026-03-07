"""
cone_fft_invertible.py — cone-FFT v29 cos+sin 完全可逆パイプライン
（引継書v30対応・parseval_core.py依存版）
"""
import numpy as np
from parseval_core import (
    build_fractal_integer_signal,
    build_fractal_integer_signal_inv,
)

BASE = 3120

def _split_TU(sig):
    return sig[0::2], sig[1::2]

def _merge_TU(T, U):
    out = [0] * (len(T) + len(U))
    for i, v in enumerate(T): out[2*i]   = v
    for i, v in enumerate(U): out[2*i+1] = v
    return out

def cone_fft_forward(f, B=1):
    N = len(f)
    assert N == BASE
    F = np.fft.fft(np.array(f, dtype=np.int64))
    cos_proj = [round(BASE * F[k].real) for k in range(N)]
    sin_proj = [round(BASE * F[k].imag) for k in range(N)]
    T_cos, U_cos = build_fractal_integer_signal(B, *_split_TU(cos_proj))
    T_sin, U_sin = build_fractal_integer_signal(B, *_split_TU(sin_proj))
    return {"cos_proj": cos_proj, "sin_proj": sin_proj,
            "T_cos": T_cos, "U_cos": U_cos,
            "T_sin": T_sin, "U_sin": U_sin, "B": B}

def cone_fft_inverse(T_cos, U_cos, T_sin, U_sin, B=1):
    N = len(T_cos) + len(U_cos)
    assert N == BASE
    T_cos0, U_cos0 = build_fractal_integer_signal_inv(B, T_cos, U_cos)
    T_sin0, U_sin0 = build_fractal_integer_signal_inv(B, T_sin, U_sin)
    cos_proj = _merge_TU(T_cos0, U_cos0)
    sin_proj = _merge_TU(T_sin0, U_sin0)
    F_rec = np.array([complex(cos_proj[k], sin_proj[k]) / BASE for k in range(N)])
    return [round(v) for v in np.fft.ifft(F_rec).real]

def cone_fft_roundtrip(f, B=1):
    res = cone_fft_forward(f, B)
    return cone_fft_inverse(res["T_cos"], res["U_cos"],
                            res["T_sin"], res["U_sin"], B)
