import argparse
import hashlib
import re
import zlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import mne
import mne_bids
import numpy as np
import pandas as pd

SFREQ = 512
N_SAMPLES = 10240
VHDR_EXPECTED = {
    "DataFormat": "BINARY",
    "DataOrientation": "MULTIPLEXED",
    "NumberOfChannels": "2",
    "SamplingInterval": "1953.125",
    "BinaryFormat": "IEEE_FLOAT_32",
    "Ch1": "x,,1,n/a",
    "Ch2": "y,,1,n/a",
}
# Fig. 1 caption: panel -> (acq, S, Ux, Uy, B)
CAPTION = {
    "a": ("nonfocal", 0, 0, 0, 0),
    "b": ("nonfocal", 0, 0, 1, 1),
    "c": ("nonfocal", 1, 1, 1, 1),
    "d": ("focal", 0, 0, 0, 0),
    "e": ("focal", 0, 1, 1, 1),
    "f": ("focal", 1, 0, 0, 0),
}
OUTCOMES = [
    "stationarity_test_S",
    "randomness_test_Ux",
    "randomness_test_Uy",
    "nonlinear_independence_test_B",
]
TOKEN = re.compile(
    r"\((?:\\.|[^\\)])*\)|[-+]?(?:\d+\.?\d*|\.\d+)|/[^\s/\[\]()<>{}]+|[A-Za-z'\"*]+|[\[\]]"
)


def figure_stream(pdf: bytes) -> str:
    streams = []
    for match in re.finditer(rb"/FlateDecode[^>]*>>\s*stream\r?\n", pdf):
        try:
            streams.append(zlib.decompressobj().decompress(pdf[match.end() :]))
        except zlib.error:
            continue
    return max(streams, key=lambda s: s.count(b" l\n")).decode("latin1")


def parse_paths_and_labels(stream: str) -> tuple[list, list]:
    paths, labels, stack = [], [], []
    path, rect, clip, tm = None, None, None, (0.0, 0.0)
    for tok in TOKEN.findall(stream):
        if tok[0] in "-+.0123456789(/":
            stack.append(tok[1:-1] if tok[0] == "(" else tok)
            continue
        if tok == "m":
            path = [(float(stack[-2]), float(stack[-1]))]
            paths.append((path, clip))
        elif tok in ("l", "c"):
            path.append((float(stack[-2]), float(stack[-1])))
        elif tok == "re":
            rect = tuple(float(v) for v in stack[-4:])
        elif tok == "W":
            clip = rect
        elif tok == "Q":
            clip = None
        elif tok == "Tm":
            tm = (float(stack[-2]), float(stack[-1]))
        elif tok == "Tj" and stack:
            labels.append((stack[-1], tm))
        stack = []
    return paths, labels


def digitize_fig1(pdf: bytes) -> dict:
    paths, labels = parse_paths_and_labels(figure_stream(pdf))
    clipped: dict = {}
    for p, clip in paths:
        if clip is not None:
            chunks = clipped.setdefault(clip, [])
            chunks.extend(p[1:] if chunks and chunks[-1] == p[0] else p)
    traces = sorted(
        [(np.array(p), clip) for clip, p in clipped.items() if len(p) > 1000],
        key=lambda t: -t[1][1],
    )
    assert len(traces) == 12, f"expected 12 traces (6 panels x 2), got {len(traces)}"
    ticks = [np.array(p) for p, _ in paths if len(p) == 2]
    hticks = [p for p in ticks if p[0, 1] == p[1, 1] and abs(p[1, 0] - p[0, 0]) < 5]
    vticks = [p for p in ticks if p[0, 0] == p[1, 0] and abs(p[1, 1] - p[0, 1]) < 5]
    bottom = min(clip[1] for _, clip in traces)
    xlabels = sorted((pos[0], txt) for txt, pos in labels if bottom - 15 < pos[1] < bottom - 2)
    assert [t for _, t in xlabels] == [str(s) for s in range(0, 21, 2)], xlabels
    out = {}
    for i, (vertices, (cx, cy, cw, ch)) in enumerate(traces):
        lo, hi = cy - 1, cy + ch + 1
        xs = sorted({round(p[0, 0], 2) for p in vticks if lo <= p[0, 1] <= hi})
        assert len(xs) == 11, f"trace {i}: expected 11 time ticks, got {len(xs)}"
        ys = sorted({p[0, 1] for p in hticks if abs(p[0, 0] - cx) < 1 and lo <= p[0, 1] <= hi})
        ylabels = [
            (pos[1], txt)
            for txt, pos in labels
            if re.fullmatch(r"\D?\d+", txt) and pos[0] < cx and lo <= pos[1] + 3 <= hi
        ]
        tick_text = {y: min(ylabels, key=lambda lab: abs(lab[0] - y))[1] for y in ys}
        (y0,) = [y for y, txt in tick_text.items() if txt == "0"]
        uv = np.array(
            [int(re.sub(r"\D", "", txt)) * np.sign(y - y0) for y, txt in tick_text.items()]
        )
        pt = np.array([y - y0 for y in tick_text])
        uv_per_pt = float(pt @ uv / (pt @ pt))
        assert np.allclose(pt * uv_per_pt, uv, atol=0.5), (i, tick_text)
        out[("abcdef"[i // 2], "xy"[i % 2])] = {
            "t": (vertices[:, 0] - xs[0]) / (xs[-1] - xs[0]) * 20.0,
            "uv": (vertices[:, 1] - y0) * uv_per_pt,
            "ticks": "/".join(str(int(v)) for v in sorted(uv)),
        }
    return out


def markdown(df: pd.DataFrame) -> str:
    rows = [list(df.columns), ["---"] * len(df.columns), *df.astype(str).values.tolist()]
    return "\n".join("| " + " | ".join(map(str, r)) + " |" for r in rows)


def run_path(root: Path, acq: str, index: int, ext: str) -> Path:
    name = f"sub-pooled_task-interictal_acq-{acq}_run-{index:04d}_ieeg{ext}"
    return root / "sub-pooled" / "ieeg" / name


def load_class(root: Path, acq: str) -> np.ndarray:
    return np.stack(
        [
            np.fromfile(run_path(root, acq, i, ".eeg"), "<f4").reshape(N_SAMPLES, 2)
            for i in range(1, 3751)
        ]
    )


def match_panels(root: Path, fig: dict, scans: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for acq in ("nonfocal", "focal"):
        data = load_class(root, acq)
        for panel, (panel_acq, *outcomes) in CAPTION.items():
            if panel_acq != acq:
                continue
            r = {}
            for ch, col in (("x", 0), ("y", 1)):
                trace = fig[(panel, ch)]
                k = np.clip(np.round(trace["t"] * SFREQ).astype(int), 0, N_SAMPLES - 1)
                cand = data[:, k, col].astype(float)
                cand -= cand.mean(1, keepdims=True)
                v = trace["uv"] - trace["uv"].mean()
                r[ch] = cand @ v / (np.linalg.norm(cand, axis=1) * np.linalg.norm(v))
            score = (r["x"] + r["y"]) / 2
            best, second = np.argsort(score)[::-1][:2]
            pub = scans[(scans.acq == acq) & (scans.pair_index == best + 1)][OUTCOMES]
            row = {"panel": panel, "acq": acq, "pair_index": best + 1}
            for ch, col in (("x", 0), ("y", 1)):
                trace = fig[(panel, ch)]
                k = np.clip(np.round(trace["t"] * SFREQ).astype(int), 0, N_SAMPLES - 1)
                stored = data[best, :, col].astype(float)
                slope, _ = np.polyfit(stored[k], trace["uv"], 1)
                row[f"r_{ch}"] = round(float(r[ch][best]), 4)
                row[f"slope_{ch}"] = round(float(slope), 4)
                row[f"range_ratio_{ch}"] = round(float(np.ptp(trace["uv"]) / np.ptp(stored)), 4)
                row[f"ticks_uV_{ch}"] = trace["ticks"]
            row["next_best_mean_r"] = round(float(score[second]), 4)
            row["caption_outcomes_match"] = tuple(pub.iloc[0]) == tuple(outcomes)
            rows.append(row)
    return pd.DataFrame(rows).sort_values("panel")


def check_run(args: tuple) -> dict:
    root, clone, acq, index = args
    res = {"acq": acq}
    for ext in (".eeg", ".vhdr", ".vmrk"):
        tracked = run_path(clone, acq, index, ext).read_bytes()
        blob = run_path(root, acq, index, ext).read_bytes()
        key = re.fullmatch(rb"/annex/objects/SHA256E-s(\d+)--([0-9a-f]{64})\.\w+\n?", tracked)
        if key is None:
            res[f"{ext[1:]}_matches_head"] = blob == tracked
        else:
            size, digest = int(key[1]), key[2].decode()
            res[f"{ext[1:]}_matches_head"] = (
                len(blob) == size and hashlib.sha256(blob).hexdigest() == digest
            )
    vhdr = dict(
        line.split("=", 1)
        for line in run_path(root, acq, index, ".vhdr").read_text().splitlines()
        if "=" in line and not line.startswith(";")
    )
    res["vhdr_ok"] = all(vhdr.get(k) == v for k, v in VHDR_EXPECTED.items())
    stored = np.fromfile(run_path(root, acq, index, ".eeg"), "<f4").reshape(N_SAMPLES, 2).T
    bp = mne_bids.BIDSPath(
        root=root,
        subject="pooled",
        task="interictal",
        acquisition=acq,
        run=f"{index:04d}",
        datatype="ieeg",
        suffix="ieeg",
        extension=".vhdr",
    )
    raw = mne_bids.read_raw_bids(bp, verbose="error")
    res["read_equals_stored"] = bool(np.array_equal(raw.get_data(), stored.astype(float)))
    res["types"] = ",".join(raw.get_channel_types())
    res["units"] = ",".join(str(c["unit"]) for c in raw.info["chs"])
    res["orig_units"] = ",".join(raw._orig_units.values())
    return res


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path, help="BIDS root: head sidecars + binaries")
    parser.add_argument("clone", type=Path, help="git working tree at head (annex pointers)")
    parser.add_argument("pdf", type=Path, help="published PDF, hdl:10230/43557")
    parser.add_argument("--jobs", type=int, default=1)
    args = parser.parse_args()
    pdf = args.pdf.read_bytes()
    print(f"mne {mne.__version__}, mne_bids {mne_bids.__version__}, numpy {np.__version__}")
    print(f"pdf md5 {hashlib.md5(pdf).hexdigest()} sha256 {hashlib.sha256(pdf).hexdigest()}")
    channels = args.root / "sub-pooled/ieeg/sub-pooled_task-interictal_channels.tsv"
    channels = pd.read_csv(channels, sep="\t", keep_default_na=False)
    print(markdown(channels[["name", "type", "units"]]))
    scans = pd.read_csv(args.root / "sub-pooled/sub-pooled_scans.tsv", sep="\t")
    fig = digitize_fig1(pdf)
    print(f"Fig. 1 traces digitized: {len(fig)}")
    print(markdown(match_panels(args.root, fig, scans)))
    runs = [
        (args.root, args.clone, acq, i) for acq in ("focal", "nonfocal") for i in range(1, 3751)
    ]
    with ProcessPoolExecutor(args.jobs) as pool:
        checks = pd.DataFrame(pool.map(check_run, runs, chunksize=50))
    summary = checks.groupby("acq").agg(
        runs=("acq", "size"),
        eeg_matches_head=("eeg_matches_head", "sum"),
        vhdr_matches_head=("vhdr_matches_head", "sum"),
        vmrk_matches_head=("vmrk_matches_head", "sum"),
        vhdr_ok=("vhdr_ok", "sum"),
        read_equals_stored=("read_equals_stored", "sum"),
        types=("types", lambda s: "; ".join(sorted(set(s)))),
        units=("units", lambda s: "; ".join(sorted(set(s)))),
        orig_units=("orig_units", lambda s: "; ".join(sorted(set(s)))),
    )
    print(markdown(summary.reset_index()))


if __name__ == "__main__":
    main()
