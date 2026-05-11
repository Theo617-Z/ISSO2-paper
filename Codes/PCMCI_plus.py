from pathlib import Path

import numpy as np
import pandas as pd
from tigramite import data_processing as pp
from tigramite.independence_tests.parcorr import ParCorr
from tigramite.pcmci import PCMCI


BASE_DIR = Path(__file__).resolve().parents[1]
DATA_FILE = BASE_DIR / "Data" / "analysis_data.csv"


def get_data(data_all, var_names):
    return data_all[var_names].to_numpy()


def get_selected_links(var_names, tau_min, tau_max):
    season_idx = np.where(np.array(var_names) == r"$Season$")[0]
    season_idx = int(season_idx[0]) if len(season_idx) else None
    temp_idx = np.where(np.array(var_names) == r"$Temp.$")[0]
    temp_idx = int(temp_idx[0]) if len(temp_idx) else None
    humid_idx = np.where(np.array(var_names) == r"$Humid.$")[0]
    humid_idx = int(humid_idx[0]) if len(humid_idx) else None

    selected_links = {}
    for idx, var in enumerate(var_names):
        if var in {r"$IS$", r"$HEM$"}:
            selected_links[idx] = [
                (other_idx, -tau)
                for other_idx, other_var in enumerate(var_names)
                for tau in range(tau_min, tau_max + 1)
                if other_var != r"$Season$"
            ]
        elif var == r"$Humid.$":
            selected_links[idx] = [(idx, -tau) for tau in range(tau_min, tau_max + 1)]
            if temp_idx is not None:
                selected_links[idx] = [(temp_idx, -tau) for tau in range(max(1, tau_min), tau_max + 1)]
        elif var == r"$Temp.$":
            selected_links[idx] = [(idx, -tau) for tau in range(tau_min, tau_max + 1)]
            if humid_idx is not None:
                selected_links[idx] = [(humid_idx, -tau) for tau in range(max(1, tau_min), tau_max + 1)]
        elif var == r"$Season$":
            selected_links[idx] = [(idx, -tau) for tau in range(tau_min, tau_max + 1)]
        elif var == r"$SO_{2}$":
            selected_links[idx] = [
                (other_idx, -tau)
                for other_idx, other_var in enumerate(var_names)
                for tau in range(tau_min, tau_max + 1)
                if other_var not in {r"$PM_{2.5}$", r"$IS$", r"$HEM$", r"$Season$"}
            ]
        else:
            selected_links[idx] = [
                (other_idx, -tau)
                for other_idx, other_var in enumerate(var_names)
                for tau in range(tau_min, tau_max + 1)
                if other_var not in {r"$IS$", r"$HEM$", r"$Season$"}
            ]

        if var != r"$Season$" and season_idx is not None:
            selected_links[idx].append((season_idx, -tau_min))

    return selected_links


def selected_links_to_link_assumptions(selected_links, n_vars):
    link_assumptions = {j: {} for j in range(n_vars)}
    for j, links in selected_links.items():
        for i, neg_tau in links:
            tau = -neg_tau
            link_assumptions[j][(i, 0 if tau == 0 else -tau)] = "o?o" if tau == 0 else "-?>"
    return link_assumptions


def clean_var_label(x):
    return {
        r"$IS$": "IS",
        r"$HEM$": "HEM",
        r"$SO_{2}$": "SO2",
        r"$O_{3}$": "O3",
        r"$PM_{2.5}$": "PM2.5",
        r"$Temp.$": "Temp",
        r"$Humid.$": "Humid",
    }.get(x, x)


def export_pcmci_long(results, var_names, tau_max, out_file):
    rows = []
    for i, source in enumerate(var_names):
        for j, target in enumerate(var_names):
            for tau in range(tau_max + 1):
                graph_value = str(results["graph"][i, j, tau]).strip()
                rows.append(
                    {
                        "source_idx": i,
                        "target_idx": j,
                        "source": source,
                        "target": target,
                        "source_clean": clean_var_label(source),
                        "target_clean": clean_var_label(target),
                        "lag": tau,
                        "effect": results["val_matrix"][i, j, tau],
                        "p_value": results["p_matrix"][i, j, tau],
                        "graph": graph_value,
                        "is_link": graph_value != "",
                        "sig_005": results["p_matrix"][i, j, tau] < 0.05,
                        "sig_010": results["p_matrix"][i, j, tau] < 0.10,
                    }
                )
    out_df = pd.DataFrame(rows)
    out_df.to_csv(out_file, index=False)
    return out_df


def run_pcmci_plus(data_all, var_names, tau_min=0, tau_max=3, pc_alpha=0.05):
    dataframe = pp.DataFrame(get_data(data_all, var_names), var_names=var_names, missing_flag=999.0)
    pcmci = PCMCI(dataframe=dataframe, cond_ind_test=ParCorr(), verbosity=0)
    selected_links = get_selected_links(var_names, tau_min, tau_max)
    link_assumptions = selected_links_to_link_assumptions(selected_links, len(var_names))
    return pcmci.run_pcmciplus(
        tau_min=tau_min,
        tau_max=tau_max,
        pc_alpha=pc_alpha,
        link_assumptions=link_assumptions,
    )


if __name__ == "__main__":
    data_all = pd.read_csv(DATA_FILE)
    data_all.rename(
        columns={
            "isc": r"$IS$",
            "hem": r"$HEM$",
            "o3": r"$O_{3}$",
            "fsp": r"$PM_{2.5}$",
            "so2": r"$SO_{2}$",
            "temp": r"$Temp.$",
            "rh": r"$Humid.$",
        },
        inplace=True,
    )

    var_names = [r"$IS$", r"$SO_{2}$", r"$O_{3}$", r"$PM_{2.5}$", r"$Temp.$", r"$Humid.$"]
    results = run_pcmci_plus(data_all, var_names)
    export_pcmci_long(results, var_names, tau_max=3, out_file=BASE_DIR / "pcmci_long_Overall_Year.csv")
    print("Wrote", BASE_DIR / "pcmci_long_Overall_Year.csv")
