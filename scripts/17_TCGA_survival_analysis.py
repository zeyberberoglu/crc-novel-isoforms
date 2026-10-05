import os
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import norm
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
from lifelines.plotting import add_at_risk_counts
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")

# ==============================
# User-defined paths
# ==============================

OUT = "/path/to/TCGA/survival"

MED = os.path.join(
    OUT,
    "median"
)

OPT = os.path.join(
    OUT,
    "optimal"
)

MIN_DETECT = 0.25

for d in (MED, OPT):
    os.makedirs(
        os.path.join(
            d,
            "km_plots"
        ),
        exist_ok=True
    )

# ==============================
# Load expression and clinical data
# ==============================

expr = pd.read_csv(
    os.path.join(
        OUT,
        "lncRNA_TPM_matrix.csv"
    )
)

clin = pd.read_csv(
    os.path.join(
        OUT,
        "clinical_parsed.csv"
    )
)

# ==============================
# Build survival variables
# ==============================

clin["time"] = clin["days_to_death"].where(
    clin["vital_status"].str.lower() == "dead",
    clin["days_to_last_follow_up"]
)

clin["event"] = (
    clin["vital_status"]
    .str.lower()
    .eq("dead")
    .astype(int)
)

clin["age_years"] = (
    clin["age_at_diagnosis"]
    / 365.25
)

clin = clin[
    (clin["time"].notna())
    & (clin["time"] > 0)
].set_index(
    "submitter_id"
)

# ==============================
# Prepare expression matrix
# ==============================

meta = expr[
    [
        "ens_id",
        "gene_name"
    ]
]

mat = (
    expr
    .set_index("ens_id")
    .drop(
        columns=["gene_name"]
    )
)

mat.columns = [
    c[:12]
    for c in mat.columns
]

common = [
    s
    for s in mat.columns
    if s in clin.index
]

mat = mat[
    common
]

surv = clin.loc[
    common,
    [
        "time",
        "event",
        "age_years"
    ]
]

# ==============================
# Detection filtering
# ==============================

det = (
    mat > 0
).mean(
    axis=1
)

excl = (
    meta
    .set_index("ens_id")
    .loc[
        det[
            det < MIN_DETECT
        ].index
    ]
    .reset_index()
)

excl["detection_frac"] = det.loc[
    excl["ens_id"]
].values

excl.to_csv(
    os.path.join(
        OUT,
        "excluded_low_expression.csv"
    ),
    index=False
)

mat = mat.loc[
    det[
        det >= MIN_DETECT
    ].index
]

# ==============================
# Log-transform TPM
# ==============================

L = np.log2(
    mat.astype(float)
    + 1.0
)

surv.join(
    L.T
).to_csv(
    os.path.join(
        OUT,
        "analysis_matrix.csv"
    )
)

print(
    f"patients: {len(common)} | "
    f"events: {int(surv['event'].sum())} | "
    f"tested: {len(L)} | "
    f"excluded (<{MIN_DETECT:.0%} detected): {len(excl)}"
)

# ==============================
# Max-statistic adjusted p-value
# ==============================

def maxstat_p(M, eps=0.2):

    M = abs(M)

    if M <= 1.0:
        return 1.0

    d = norm.pdf(M)

    p = (
        d
        * (M - 1.0 / M)
        * np.log(
            ((1 - eps) ** 2)
            / (eps ** 2)
        )
        + 4 * d / M
    )

    return float(
        min(
            max(
                p,
                0.0
            ),
            1.0
        )
    )


# ==============================
# Cutpoint methods
# ==============================

def median_cut(x, t, e):

    return (
        float(
            np.median(x)
        ),
        np.nan
    )


def optimal_cut(
    x,
    t,
    e,
    eps=0.2
):

    lo, hi = np.quantile(
        x,
        [
            eps,
            1 - eps
        ]
    )

    best_z = 0.0

    best_c = float(
        np.median(x)
    )

    for c in np.unique(
        x[
            (x >= lo)
            & (x <= hi)
        ]
    ):

        g = x > c

        if (
            g.sum() < 5
            or (~g).sum() < 5
        ):
            continue

        z = np.sqrt(
            logrank_test(
                t[g],
                t[~g],
                e[g],
                e[~g]
            ).test_statistic
        )

        if z > best_z:
            best_z = z
            best_c = c

    return (
        best_c,
        best_z
    )


# ==============================
# Survival analysis
# ==============================

def analyze(
    cut_fn,
    folder,
    corrected
):

    rows = []

    for eid, row in L.iterrows():

        g = meta.loc[
            meta.ens_id == eid,
            "gene_name"
        ].iloc[0]

        df = pd.DataFrame(
            {
                "time": surv["time"].values,
                "event": surv["event"].values,
                "expr": row.values
            }
        ).dropna()

        if (
            df["expr"].nunique() < 3
            or df["event"].sum() < 5
        ):
            continue

        cut, stat = cut_fn(
            df["expr"].values,
            df["time"].values,
            df["event"].values
        )

        high = (
            df["expr"]
            > cut
        )

        d = dict(
            ens_id=eid,
            gene_name=g,
            n=len(df),
            events=int(
                df["event"].sum()
            ),
            cutpoint_log2tpm=float(cut),
            n_high=int(
                high.sum()
            ),
            n_low=int(
                (~high).sum()
            )
        )

        if (
            high.sum() >= 10
            and (~high).sum() >= 10
        ):

            lr = logrank_test(
                df["time"][high],
                df["time"][~high],
                df["event"][high],
                df["event"][~high]
            )

            c2 = CoxPHFitter().fit(
                df.assign(
                    h=high.astype(int)
                )[
                    [
                        "time",
                        "event",
                        "h"
                    ]
                ],
                "time",
                "event"
            )

            lo, hi = np.exp(
                c2
                .confidence_intervals_
                .loc["h"]
            ).values

            d.update(
                HR=float(
                    np.exp(
                        c2.params_["h"]
                    )
                ),
                HR_lo=lo,
                HR_hi=hi,
                logrank_p=lr.p_value,
                p_report=(
                    maxstat_p(stat)
                    if corrected
                    else lr.p_value
                )
            )

        else:

            d.update(
                HR=np.nan,
                HR_lo=np.nan,
                HR_hi=np.nan,
                logrank_p=np.nan,
                p_report=np.nan
            )

        rows.append(d)

    res = pd.DataFrame(
        rows
    )

    ok = res[
        "p_report"
    ].notna()

    res.loc[
        ok,
        "padj"
    ] = multipletests(
        res.loc[
            ok,
            "p_report"
        ],
        method="fdr_bh"
    )[1]

    res = res.sort_values(
        "p_report"
    )

    res.to_csv(
        os.path.join(
            folder,
            "survival_results.csv"
        ),
        index=False
    )

    return res


# ==============================
# Kaplan-Meier plots
# ==============================

def km(
    eid,
    g,
    cut_fn,
    folder
):

    df = pd.DataFrame(
        {
            "time": surv["time"].values,
            "event": surv["event"].values,
            "expr": L.loc[eid].values
        }
    ).dropna()

    cut, _ = cut_fn(
        df["expr"].values,
        df["time"].values,
        df["event"].values
    )

    high = (
        df["expr"]
        > cut
    )

    if (
        high.sum() < 10
        or (~high).sum() < 10
    ):
        return

    lr = logrank_test(
        df["time"][high],
        df["time"][~high],
        df["event"][high],
        df["event"][~high]
    )

    yrs = (
        df["time"]
        / 365.25
    )

    c = CoxPHFitter().fit(
        pd.DataFrame(
            {
                "time": yrs,
                "event": df["event"].values,
                "h": high.astype(int)
            }
        ),
        "time",
        "event"
    )

    hr = float(
        np.exp(
            c.params_["h"]
        )
    )

    p = lr.p_value

    fig, ax = plt.subplots(
        figsize=(4.3, 4.2),
        dpi=300
    )

    fits = []

    for sel, lab, col in [
        (
            high,
            "High",
            "#C0392B"
        ),
        (
            ~high,
            "Low",
            "#2C3E50"
        )
    ]:

        k = KaplanMeierFitter().fit(
            yrs[sel],
            df["event"].values[
                sel.values
            ],
            label=(
                f"{lab} "
                f"(n={int(sel.sum())})"
            )
        )

        k.plot_survival_function(
            ax=ax,
            ci_show=False,
            color=col,
            linewidth=1.8
        )

        fits.append(k)

    ax.set_xlabel(
        "Time (years)"
    )

    ax.set_ylabel(
        "Overall survival"
    )

    ax.set_title(
        g,
        style="italic"
    )

    ax.set_ylim(
        0,
        1.02
    )

    ax.spines[
        [
            "top",
            "right"
        ]
    ].set_visible(
        False
    )

    ax.legend(
        frameon=False,
        loc="lower left",
        fontsize=9
    )

    ax.text(
        0.97,
        0.96,
        f"HR = {hr:.2f}\n"
        "log-rank P = "
        + (
            f"{p:.1e}"
            if p < 0.001
            else f"{p:.3f}"
        ),
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9
    )

    add_at_risk_counts(
        *fits,
        ax=ax,
        rows_to_show=[
            "At risk"
        ]
    )

    plt.tight_layout()

    s = "".join(
        ch
        if ch.isalnum()
        else "_"
        for ch in str(g)
    )

    for ext in (
        "pdf",
        "png"
    ):

        fig.savefig(
            os.path.join(
                folder,
                "km_plots",
                f"KM_{s}.{ext}"
            ),
            bbox_inches="tight"
        )

    plt.close(fig)


# ==============================
# Run median and optimal analyses
# ==============================

for cut_fn, folder, name, corrected in [

    (
        median_cut,
        MED,
        "median",
        False
    ),

    (
        optimal_cut,
        OPT,
        "optimal",
        True
    )

]:

    res = analyze(
        cut_fn,
        folder,
        corrected
    )

    for _, r in res.iterrows():

        km(
            r["ens_id"],
            r["gene_name"],
            cut_fn,
            folder
        )

    print(
        f"{name}: "
        f"tested={len(res)} | "
        f"padj<0.05: "
        f"{int((res['padj'] < 0.05).sum())}"
    )
