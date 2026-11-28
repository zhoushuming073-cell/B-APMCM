"""Nested LOHO sweep of causal delayed-gas feedback for Q3.

Reuses the repository's frozen V5 implementation and changes only the score
start and feedback interval/gain grid. Run from repository root.
"""
from __future__ import annotations
import argparse
import importlib.util
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[2]


def load_v5(cache: Path, scratch: Path):
    # The frozen V5 module parses its own CLI at import. Supply its exact args,
    # then reuse its functions without invoking its main().
    sys.path.insert(0, str(ROOT / "src"))
    old_argv = sys.argv[:]
    sys.argv = [
        "v5_causal_spectral.py", "--cache", str(cache),
        "--out", str(scratch), "--lag-scheme", "q2_fixed"
    ]
    spec = importlib.util.spec_from_file_location(
        "q3_v5_audit_module", ROOT / "src/apmcm/q3/v5_causal_spectral.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    sys.argv = old_argv
    return module


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--cache", type=Path, default=Path("data/interim/q3_feedback_cache"))
    p.add_argument("--out", type=Path, default=Path("results/q3/feedback_audit/run"))
    p.add_argument("--heat", default="all", help="all or comma-separated heat IDs")
    p.add_argument("--gas", choices=["CO2", "CO", "both"], default="both")
    p.add_argument("--score-start-s", type=float, default=120.0)
    p.add_argument("--intervals", default="5,15,30,60,120,180")
    p.add_argument("--gains", default="0.25,0.5,0.75,1.0")
    args = p.parse_args()
    if args.score_start_s < 60:
        raise ValueError(
            "score-start-s must be >=60 s so every candidate lead has a "
            "nonnegative causal input index"
        )
    args.out.mkdir(parents=True, exist_ok=True)
    scratch = args.out / "_v5_scratch"
    m = load_v5(args.cache, scratch)

    # V5 defaults to a 120 s score start. Replace only the score index; all
    # candidate construction, lag selection, and nested model selection stay intact.
    for h in m.hs:
        m.D[h]["idx"] = np.arange(
            int(round(args.score_start_s * 2)), len(m.D[h]["t"])
        )
    for fn in (m.q2_stats, m.lead, m.features, m.template, m.prior, m.candidates):
        fn.cache_clear()

    intervals = [int(x) for x in args.intervals.split(",")]
    gains = [float(x) for x in args.gains.split(",")]
    heats = m.hs if args.heat == "all" else tuple(
        x for x in args.heat.split(",") if x
    )
    gas_ids = [0, 1] if args.gas == "both" else [0 if args.gas == "CO2" else 1]
    rows = []

    with threadpool_limits(limits=1):
        for h in heats:
            for g in gas_ids:
                train = tuple(j for j in m.hs if j != h)
                inner_rows = []
                for val in train:
                    inner_train = tuple(j for j in train if j != val)
                    pred, d = m.candidates(inner_train, val, g)
                    y = m.D[val]["y"][m.D[val]["idx"], g]
                    for name, v in pred.items():
                        fam = m.family(name)
                        inner_rows.append(dict(
                            model=name, family=fam, inner_heat=m.letters[val],
                            **m.score(y, v)
                        ))
                        if fam == "fusion":
                            for interval in intervals:
                                for gain in gains:
                                    fb_name = f"{name}__fb{interval}_{gain}"
                                    fb_pred = m.feedback(v, y, d, interval, gain)
                                    inner_rows.append(dict(
                                        model=fb_name, family=f"feedback{interval}",
                                        inner_heat=m.letters[val],
                                        **m.score(y, fb_pred)
                                    ))

                inner = pd.DataFrame(inner_rows)
                chosen = (
                    inner.groupby(["family", "model"]).r2.mean().reset_index()
                    .sort_values("r2", ascending=False)
                    .drop_duplicates("family")
                )
                pred, d = m.candidates(train, h, g)
                y = m.D[h]["y"][m.D[h]["idx"], g]
                for z in chosen.itertuples():
                    name, fam = z.model, z.family
                    if fam.startswith("feedback"):
                        base, suffix = name.rsplit("__fb", 1)
                        interval_s, gain_s = suffix.split("_")
                        v = m.feedback(
                            pred[base], y, d, int(interval_s), float(gain_s)
                        )
                    else:
                        v = pred[name]
                    rows.append(dict(
                        score_start_s=args.score_start_s,
                        heat=m.letters[h], heat_id=h,
                        gas=("CO2", "CO")[g], family=fam, model=name,
                        lead_s=d, inner_r2=z.r2, n=len(y),
                        **m.score(y, v)
                    ))
                pd.DataFrame(rows).to_csv(
                    args.out / "fold_metrics_partial.csv", index=False
                )
                print("DONE", m.letters[h], ("CO2", "CO")[g], flush=True)

    result = pd.DataFrame(rows).sort_values(["gas", "family", "heat"])
    result.to_csv(
        args.out / "fold_metrics.csv", index=False, float_format="%.8g"
    )
    print(result.groupby(["gas", "family"]).agg(
        heats=("heat", "nunique"),
        macro_r2=("r2", "mean"),
        macro_rmse=("rmse", "mean"),
        min_r2=("r2", "min"),
    ).to_string())


if __name__ == "__main__":
    main()
