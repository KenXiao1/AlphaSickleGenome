"""CLI entry point for Phase 1: AlphaGenome 2 BCL11A +58 kb enhancer variant reproduction.

Usage:
  Full reproduction (uses local cache if available; queries AlphaGenome 2 API if missing):
    uv run --env-file .env python src/reproduce_bcl11a.py

  Force fresh API queries (re-populates cache):
    uv run --env-file .env python src/reproduce_bcl11a.py --refresh-cache

  Regenerate all figures offline strictly from saved processed outputs (no API key needed):
    uv run python src/reproduce_bcl11a.py --plots-only
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from bcl11a_repro.pipeline import (  # noqa: E402
    regenerate_figures_from_processed,
    run_reproduction_pipeline,
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Reproduce AlphaGenome 2 Fig. 4j-k BCL11A +58 kb enhancer variant predictions."
    )
    parser.add_argument(
        "--refresh-cache",
        action="store_true",
        help="Ignore cached API responses and re-query live AlphaGenome 2 endpoints.",
    )
    parser.add_argument(
        "--offline",
        action="store_true",
        help="Run full validation, cross-validation, and plotting from cached raw files without API key.",
    )
    parser.add_argument(
        "--plots-only",
        action="store_true",
        help="Regenerate Figures A-D strictly from saved files in data/processed/bcl11a_ag2/ without API key.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    dotenv_path = ROOT / ".env"
    cache_dir = ROOT / "data" / "cache" / "bcl11a_ag2"
    raw_export_dir = ROOT / "data" / "raw" / "bcl11a_ag2"
    processed_dir = ROOT / "data" / "processed" / "bcl11a_ag2"
    figures_dir = ROOT / "docs" / "figures" / "bcl11a_ag2"

    if args.plots_only:
        saved = regenerate_figures_from_processed(
            processed_dir=processed_dir,
            figures_dir=figures_dir,
        )
        print(
            json.dumps(
                {
                    "mode": "plots_only",
                    "figures": {k: {fmt: str(p) for fmt, p in v.items()} for k, v in saved.items()},
                },
                indent=2,
            )
        )
        return 0

    summary = run_reproduction_pipeline(
        dotenv_path=dotenv_path,
        cache_dir=cache_dir,
        raw_export_dir=raw_export_dir,
        processed_dir=processed_dir,
        figures_dir=figures_dir,
        offline=args.offline,
        refresh_cache=args.refresh_cache,
        data_provenance_label="REAL_ALPHAGENOME_2_API",
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
