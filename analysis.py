from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class DatasetSchema:
    columns: tuple[str, ...] = (
        "Key",
        "Cell_Type",
        "Species",
        "Functional_State",
        "Diff",
        "Vm_mean_mV",
        "Uncert_mV",
        "Method",
        "Ref",
        "Kingdom",
    )


class UncertaintyParser:
    _num = re.compile(r"[-+]?\d*\.?\d+")

    @classmethod
    def to_float(cls, value: object) -> Optional[float]:
        if value is None:
            return None
        if isinstance(value, (int, float)) and not pd.isna(value):
            return float(value)
        s = str(value).strip()
        if not s:
            return None
        if s.lower() in {"range", "rng", "varies"}:
            return None
        m = cls._num.search(s)
        return float(m.group()) if m else None

    @staticmethod
    def is_sem(value: object) -> bool:
        if value is None:
            return False
        return "se" in str(value).lower()


class BioelectricData:
    def __init__(self, csv_path: Path) -> None:
        self.csv_path = csv_path
        self.schema = DatasetSchema()
        self.df: pd.DataFrame = pd.DataFrame()

    def ensure_exists(self) -> None:
        if self.csv_path.exists():
            return
        df = pd.DataFrame({c: pd.Series(dtype="string")
                          for c in self.schema.columns})
        df.to_csv(self.csv_path, index=False)

    def load(self) -> pd.DataFrame:
        self.ensure_exists()
        df = pd.read_csv(self.csv_path, dtype=str, engine="python").fillna("")
        missing = [c for c in self.schema.columns if c not in df.columns]
        if missing:
            raise ValueError(f"dataset.csv missing columns: {missing}")

        df["Diff"] = pd.to_numeric(df["Diff"], errors="coerce")
        df["Vm_mean_mV"] = pd.to_numeric(df["Vm_mean_mV"], errors="coerce")
        df["Uncert_mV_num"] = df["Uncert_mV"].apply(UncertaintyParser.to_float)

        self.df = df
        return df

    def validate(self) -> None:
        if self.df.empty:
            raise ValueError("dataset is empty; load() first")
        if self.df["Vm_mean_mV"].isna().all():
            raise ValueError("Vm_mean_mV has no numeric values")


class StateMapper:
    @staticmethod
    def for_summary(state: str) -> str:
        s = str(state).strip()
        if s in {"Cancerous", "Oncogene_induced"}:
            return "Cancerous/Oncogene-induced"
        return s

    @staticmethod
    def for_plot_legend(state: str) -> str:
        s = str(state).strip()
        if s in {"Cancerous", "Oncogene_induced"}:
            return "Cancerous/oncogene-induced"
        if s in {"Early embryo", "Early cleavage"}:
            return "Early embryo/cleavage"
        return s


class BioelectricSummary:
    def __init__(self, df: pd.DataFrame) -> None:
        self.df = df.copy()

    def stats_by_state(self) -> pd.DataFrame:
        d = self.df.copy()
        d["Vm"] = pd.to_numeric(d["Vm_mean_mV"], errors="coerce")
        d["State_group"] = d["Functional_State"].apply(StateMapper.for_summary)

        out = (
            d.groupby("State_group", dropna=False)
            .agg(n=("Vm", "count"), mean=("Vm", "mean"), sd=("Vm", "std"))
            .reset_index()
            .rename(columns={"State_group": "Functional State"})
        )
        out = out[out["n"] >= 2].copy()

        order = [
            "Excitable",
            "Progenitor",
            "Differentiated",
            "Differentiating",
            "Primordial",
            "Proliferative",
            "Cancerous/Oncogene-induced",
        ]
        out["__order"] = out["Functional State"].apply(
            lambda x: order.index(x) if x in order else 999)
        out = out.sort_values("__order").drop(
            columns="__order").reset_index(drop=True)
        return out

    @staticmethod
    def to_latex_table(stats: pd.DataFrame, label: str = "tab:summary_stats") -> str:
        lines: list[str] = []
        lines.append(r"\begin{table}[h!]")
        lines.append(r"    \centering")
        lines.append(
            r"    \caption{Summary statistics of Vm by functional state. Categories with single entries (Fertilisation, Early embryo, Early cleavage) are omitted as variance cannot be computed.}"
        )
        lines.append(fr"    \label{{{label}}}")
        lines.append(r"    \begin{tabular}{lccc}")
        lines.append(r"        \toprule")
        lines.append(
            r"        \textbf{Functional State}  & \textbf{n} & \textbf{Mean Vm (mV)} & \textbf{SD (mV)} \\")
        lines.append(r"        \midrule")
        for _, row in stats.iterrows():
            state = str(row["Functional State"])
            n = int(row["n"])
            mean = float(row["mean"])
            sd = float(row["sd"])
            lines.append(
                f"        {state:<27} & {n:<10d} & ${mean:.1f}$               & {sd:.1f}             \\\\")
        lines.append(r"        \bottomrule")
        lines.append(r"    \end{tabular}")
        lines.append(r"\end{table}")
        return "\n".join(lines)


class BioelectricPlots:
    def __init__(self, df: pd.DataFrame) -> None:
        self.df = df.copy()

    @staticmethod
    def _state_colors() -> dict[str, str]:
        return {
            "Primordial": "#4C72B0",
            "Fertilization": "#DD8452",
            "Early embryo": "#55A868",
            "Early cleavage": "#55A868",
            "Differentiating": "#8172B3",
            "Progenitor": "#937860",
            "Proliferative": "#C44E52",
            "Differentiated": "#64B5CD",
            "Excitable": "#1B9E77",
            "Cancerous": "#CC79A7",
            "Oncogene_induced": "#CC79A7",
        }

    def spectrum(self, outpath: Path, xlim: tuple[float, float] = (-110.0, 25.0)) -> None:
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches

        df = self.df.copy()
        df = df.sort_values(
            "Vm_mean_mV", ascending=False).reset_index(drop=True)
        df["YLabel"] = (
            df["Cell_Type"].astype(str)
            + " ("
            + df["Species"].astype(str)
            + ") ["
            + df["Ref"].astype(str)
            + "]"
        )

        colors = self._state_colors()
        legend_handles: dict[str, object] = {}

        fig, ax = plt.subplots(figsize=(10, 12))
        for i, row in df.iterrows():
            x = float(row["Vm_mean_mV"])
            err = row.get("Uncert_mV_num", np.nan)
            xerr = None if pd.isna(err) else float(err)

            state = str(row["Functional_State"]).strip()
            color = colors.get(state, "gray")
            legend_label = StateMapper.for_plot_legend(state)

            lw = 0.8 if UncertaintyParser.is_sem(row.get("Uncert_mV")) else 1.6
            ax.errorbar(
                x,
                i,
                xerr=xerr,
                fmt="o",
                color=color,
                capsize=3,
                markersize=6,
                linewidth=lw,
            )

            if legend_label not in legend_handles:
                legend_handles[legend_label] = mpatches.Patch(
                    color=color, label=legend_label)

        ax.set_yticks(np.arange(len(df)))
        ax.set_yticklabels(df["YLabel"], fontsize=8.5)
        ax.invert_yaxis()
        ax.axvline(0, color="black", linestyle="--", alpha=0.35, linewidth=1.5)
        ax.set_xlabel("Resting membrane potential (mV)", fontsize=10)
        ax.set_title("Bioelectric spectrum across cell states",
                     fontsize=12, fontweight="bold")
        ax.grid(axis="x", alpha=0.25)
        ax.set_xlim(*xlim)

        legend_order = [
            "Primordial",
            "Fertilization",
            "Early embryo/cleavage",
            "Differentiating",
            "Progenitor",
            "Proliferative",
            "Differentiated",
            "Excitable",
            "Cancerous/oncogene-induced",
        ]
        handles = [legend_handles[k]
                   for k in legend_order if k in legend_handles]
        ax.legend(handles=handles, loc="lower right",
                  fontsize=8, framealpha=0.9)

        outpath.parent.mkdir(parents=True, exist_ok=True)
        plt.tight_layout()
        plt.savefig(outpath, dpi=300, bbox_inches="tight")
        plt.close(fig)


class BioelectricPipeline:
    def __init__(self, csv_path: Path, out_dir: Path) -> None:
        self.csv_path = csv_path
        self.out_dir = out_dir
        self.data = BioelectricData(csv_path)

    def run(self) -> dict[str, Path]:
        df = self.data.load()
        self.data.validate()

        spectrum_path = self.out_dir / "bioelectricspectrum.png"
        stats_csv = self.out_dir / "summarystats.csv"
        latex_path = self.out_dir / "summarystatstable.tex"

        BioelectricPlots(df).spectrum(spectrum_path)

        summary = BioelectricSummary(df)
        stats = summary.stats_by_state()
        self.out_dir.mkdir(parents=True, exist_ok=True)
        stats.to_csv(stats_csv, index=False)

        latex = summary.to_latex_table(stats)
        latex_path.write_text(latex, encoding="utf-8")

        return {
            "spectrum": spectrum_path,
            "stats_csv": stats_csv,
            "latex_table": latex_path,
        }


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="analysis.py")
    p.add_argument("--csv", type=str, default="dataset.csv")
    p.add_argument("--out", type=str, default=".")
    return p


def main(argv: Optional[Iterable[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    csv_path = Path(args.csv)
    out_dir = Path(args.out)

    outputs = BioelectricPipeline(csv_path=csv_path, out_dir=out_dir).run()
    print(f"Loaded dataset: {csv_path.resolve()}")
    for k, v in outputs.items():
        print(f"{k}: {v.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
