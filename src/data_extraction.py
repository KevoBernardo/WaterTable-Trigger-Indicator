import os
import numpy as np
import pandas as pd
from plxscripting.easy import get_equivalent


def get_material_indices(g_i, materials):
    mat_names = [str(mat.Name) for mat in g_i.Materials]
    return [i for i, name in enumerate(mat_names) if name in materials]


def extract_phase_results(go_phase, g_o, num_nodes_per_elem):
    x_nodes = np.fromstring(
        g_o.getresults(go_phase, g_o.ResultTypes.Soil.X, "node").echo()[1:-1],
        dtype=float, sep=',')
    y_nodes = np.fromstring(
        g_o.getresults(go_phase, g_o.ResultTypes.Soil.Y, "node").echo()[1:-1],
        dtype=float, sep=',')
    mat_nodes = np.fromstring(
        g_o.getresults(go_phase, g_o.ResultTypes.Soil.MaterialIndex, "node").echo()[1:-1],
        dtype=int, sep=',')

    n_nodes = len(x_nodes)
    element_id = np.arange(n_nodes) // num_nodes_per_elem
    local_node = np.arange(n_nodes) % num_nodes_per_elem

    df_nodes = pd.DataFrame({
        "element_id": element_id,
        "local_node": local_node,
        "x_node": x_nodes,
        "y_node": y_nodes,
        "mat_id": mat_nodes
    })

    x_sp = np.fromstring(
        g_o.getresults(go_phase, g_o.ResultTypes.Soil.X, "stresspoint").echo()[1:-1],
        dtype=float, sep=',')
    y_sp = np.fromstring(
        g_o.getresults(go_phase, g_o.ResultTypes.Soil.Y, "stresspoint").echo()[1:-1],
        dtype=float, sep=',')
    sat_sp = np.fromstring(
        g_o.getresults(go_phase, g_o.ResultTypes.Soil.EffSaturationAsPercentage, "stresspoint").echo()[1:-1],
        dtype=float, sep=',')
    pactive_sp = np.fromstring(
        g_o.getresults(go_phase, g_o.ResultTypes.Soil.PActive, "stresspoint").echo()[1:-1],
        dtype=float, sep=',')
    psteady_sp = np.fromstring(
        g_o.getresults(go_phase, g_o.ResultTypes.Soil.PSteady, "stresspoint").echo()[1:-1],
        dtype=float, sep=',')
    # Try to extract both MeanEffStress and SigyyE if available in this Plaxis version
    try:
        mean_eff = np.fromstring(
            g_o.getresults(go_phase, getattr(g_o.ResultTypes.Soil, "MeanEffStress"), "stresspoint").echo()[1:-1],
            dtype=float, sep=',')
    except Exception:
        mean_eff = np.full_like(psteady_sp, np.nan, dtype=float)

    try:
        sigyye = np.fromstring(
            g_o.getresults(go_phase, getattr(g_o.ResultTypes.Soil, "SigyyE"), "stresspoint").echo()[1:-1],
            dtype=float, sep=',')
    except Exception:
        sigyye = np.full_like(psteady_sp, np.nan, dtype=float)
    mat_sp = np.fromstring(
        g_o.getresults(go_phase, g_o.ResultTypes.Soil.MaterialIndex, "stresspoint").echo()[1:-1],
        dtype=int, sep=',')

    df_sp = pd.DataFrame({
        "x_sp": x_sp,
        "y_sp": y_sp,
        "saturation": sat_sp,
        "pactive": pactive_sp,
        "psteady": psteady_sp,
        "MeanEffStress": mean_eff,
        "SigyyE": sigyye,
        "mat_id": mat_sp
    })

    return df_nodes, df_sp


def save_phase_csvs(df_nodes, df_sp, phase_name, phase_id, outdir_model):
    nodes_fname = os.path.join(outdir_model, f"{phase_name}_ID{phase_id}_nodes.csv")
    sp_fname = os.path.join(outdir_model, f"{phase_name}_ID{phase_id}_stresspoints.csv")

    df_nodes.to_csv(nodes_fname, index=False)
    df_sp.to_csv(sp_fname, index=False)

    return nodes_fname, sp_fname


def extract_and_save(
    g_i,
    g_o,
    phases_to_extract,
    materials,
    num_nodes_per_elem,
    outdir_model,
):
    sel_mat_indices = get_material_indices(g_i, materials)
    os.makedirs(outdir_model, exist_ok=True)

    saved_info = []
    for phase_name in phases_to_extract:
        go_phase = None
        for phase in g_i.Phases:
            if phase.Name == phase_name:
                go_phase = get_equivalent(phase, g_o)
                phase_id = phase.Identification
                break

        if go_phase is None:
            print(f'Phase {phase_name} not found')
            continue

        df_nodes, df_sp = extract_phase_results(go_phase, g_o, num_nodes_per_elem)
        df_nodes = df_nodes[df_nodes["mat_id"].isin(sel_mat_indices)]
        df_sp = df_sp[df_sp["mat_id"].isin(sel_mat_indices)]

        nodes_fname, sp_fname = save_phase_csvs(
            df_nodes, df_sp, phase_name, phase_id, outdir_model
        )

        saved_info.append({
            "phase_name": phase_name,
            "phase_id": phase_id,
            "nodes_file": nodes_fname,
            "stresspoints_file": sp_fname
        })
        print(f"Saved nodes + stresspoints for {phase_id} [{phase_name}] in {outdir_model}")

    print("\nEXTRACCIÓN COMPLETA.")
    return saved_info
