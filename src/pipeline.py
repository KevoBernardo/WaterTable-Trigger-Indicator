from multiprocessing.util import debug
import os
import numpy as np
import pandas as pd
from .analysis import (
    build_domain,
    saturated_areas,
    hydraulic_energy,
    hydraulic_stress_state,
    ru_method,
    make_grid,
    mask_domain,
    interpolate_field,
    mask_field,
)
from .plotting import plot_saturation, plot_pactive, plot_ru, plot_debug_contour


def list_phases_from_csv(dir_raw):
    phases = sorted(
        {f.split("_ID")[0] for f in os.listdir(dir_raw) if f.endswith("_nodes.csv")}
    )
    return phases


def get_phase_files(dir_raw, phase_name):
    nodes_file = [f for f in os.listdir(dir_raw) if f.startswith(phase_name) and f.endswith("_nodes.csv")][0]
    return os.path.join(dir_raw, nodes_file), os.path.join(dir_raw, nodes_file.replace("_nodes.csv", "_stresspoints.csv"))


def extract_phase_id(nodes_path):
    file_name = os.path.basename(nodes_path)
    parts = file_name.split("_ID")
    return parts[1].split("_")[0] if len(parts) > 1 else "unknown"


def process_phase_files(phase_name, nodes_path, sp_path, config, outdir_figs, ref_df_sp=None):
    df_nodes = pd.read_csv(nodes_path)
    df_sp = pd.read_csv(sp_path)

    df_nodes = df_nodes[
        (df_nodes["x_node"] >= config["xmin"])
        & (df_nodes["x_node"] <= config["xmax"])
    ]
    df_sp = df_sp[
        (df_sp["x_sp"] >= config["xmin"])
        & (df_sp["x_sp"] <= config["xmax"])
    ]

    if df_nodes.empty or df_sp.empty:
        return None

    domain = build_domain(df_nodes, config["num_nodes_per_elem"])

    process_method = config.get("process_method", "all")
    if process_method not in {
        "saturated_area",
        "hydraulic_energy",
        "effective_stress_energy",
        "ru",
        "all",
    }:
        raise ValueError(f"Unknown process_method: {process_method}")

    area_sat = area_unsat = total = saturation_pct = None
    hydraulic_energy_value = hydraulic_debug = None
    effective_stress_energy_value = effective_stress_debug = None
    ru_index = ru_debug = None
    ru_path = None

    if process_method in {"saturated_area", "all"}:
        area_sat, area_unsat = saturated_areas(domain, df_sp)
        total = area_sat + area_unsat
        saturation_pct = 100 * area_sat / total if total > 0 else None

    H = config.get("dam_height")
    if H is None:
        H = df_nodes["y_node"].max() - df_nodes["y_node"].min()

    if process_method in {"hydraulic_energy", "all"}:
        hydraulic_energy_value, hydraulic_debug = hydraulic_energy(
            domain,
            df_sp,
            gamma_w=config.get("gamma_w", 10.0),
            H=H,
        )

    stress_field = config.get("stress_field", "SigyyE")

    if process_method in {"effective_stress_energy", "all"}:
        effective_stress_energy_value, effective_stress_debug = hydraulic_stress_state(
            domain,
            df_sp,
            stress_field=stress_field,
            gamma_w=config.get("gamma_w", 10.0),
        )

    if process_method in {"ru", "all"}:
        reference_phase = config.get("reference_phase")
        if reference_phase is None:
            raise ValueError(
                "reference_phase must be set in config for Ru calculation"
            )

        if ref_df_sp is None:
            raise ValueError(
                "Reference phase stresspoint data must be provided for Ru calculation"
            )

        include_positive_pressure_diff = config.get(
            "include_positive_pressure_diff",
            config.get("include_negative_pressure_diff", False),
        )

        ru_index, ru_debug = ru_method(
            domain,
            df_sp,
            ref_df_sp,
            stress_field=stress_field,
            include_positive_pressure_diff=include_positive_pressure_diff,
        )

    if process_method == "saturated_area":
        selected_metric = "saturation_pct"
        selected_value = saturation_pct

    elif process_method == "hydraulic_energy":
        selected_metric = "hydraulic_energy"
        selected_value = hydraulic_energy_value

    elif process_method == "effective_stress_energy":
        selected_metric = "effective_stress_energy"
        selected_value = effective_stress_energy_value

    elif process_method == "ru":
        selected_metric = "ru_index"
        selected_value = ru_index

    else:
        selected_metric = "all"
        selected_value = None

    grid_x, grid_y = make_grid(
        config["xmin"],
        config["xmax"],
        df_nodes["y_node"].min(),
        df_nodes["y_node"].max(),
        resolution=config.get("grid_resolution", 500),
    )

    grid_sat = interpolate_field(df_sp, "saturation", grid_x, grid_y)
    grid_pactive = interpolate_field(df_sp, "pactive", grid_x, grid_y)

    grid_stress = None
    ref_grid_pactive = None
    ref_grid_psteady = None
    pressure_diff_grid = None
    ru_index_grid = None
    ru_index_grid_masked = None

    if config["plot_water_table"]:
        grid_psteady = interpolate_field(df_sp, "psteady", grid_x, grid_y)

    mask = mask_domain(domain, grid_x, grid_y)

    grid_sat_masked = mask_field(grid_sat, mask)
    grid_sat_masked = np.clip(grid_sat_masked, 0, 100)

    grid_pactive_masked = mask_field(grid_pactive, mask)

    water_table_args = {
        "enabled": config["plot_water_table"],
        "grid_psteady_masked": (
            mask_field(grid_psteady, mask)
            if config["plot_water_table"]
            else None
        ),
        "label": config["water_table_label"],
        "color": config["water_table_colour"],
    }

    phase_id = extract_phase_id(nodes_path)

    sat_path = plot_saturation(
        grid_x,
        grid_y,
        grid_sat_masked,
        config["xmin"],
        config["xmax"],
        water_table_args,
        outdir_figs,
        phase_id,
        phase_name,
    )

    pactive_path = plot_pactive(
        grid_x,
        grid_y,
        grid_pactive_masked,
        config["xmin"],
        config["xmax"],
        water_table_args,
        outdir_figs,
        phase_name,
        phase_id,
    )

    if process_method in {"ru", "all"} and ref_df_sp is not None:

        grid_stress = interpolate_field(
            df_sp,
            stress_field,
            grid_x,
            grid_y,
        )

        ref_grid_pactive = interpolate_field(
            ref_df_sp,
            "pactive",
            grid_x,
            grid_y,
        )

        ref_grid_psteady = interpolate_field(
            ref_df_sp,
            "psteady",
            grid_x,
            grid_y,
        )

        pressure_diff_grid = grid_pactive - ref_grid_pactive

        ru_index_grid = np.full_like(
            pressure_diff_grid,
            np.nan,
            dtype=float,
        )

        ru_plot_mask = (
            np.isfinite(pressure_diff_grid)
            & np.isfinite(grid_stress)
            & (grid_stress != 0)
        )

        with np.errstate(divide="ignore", invalid="ignore"):
            ru_index_grid[ru_plot_mask] = (
                pressure_diff_grid[ru_plot_mask]
                / grid_stress[ru_plot_mask]
            )
            
        if config.get("mask_negative_ru_for_plot", True):
            ru_index_grid = np.where(
                ru_index_grid >= 0,
                ru_index_grid,
                np.nan,
            )
            
        if config.get("apply_ru_mask", True):
            ru_index_grid_masked = mask_field(
                ru_index_grid,
                mask,
            )
        else:
            ru_index_grid_masked = ru_index_grid

        reference_water_table_args = {
            "grid_psteady_masked": mask_field(
                ref_grid_psteady,
                mask,
            ),
            "label": f"Reference {config['water_table_label'].lower()}",
            "color": config.get(
                "reference_water_table_colour",
                "black",
            ),
        }

        water_table_args["label"] = (
            f"Current {config['water_table_label'].lower()}"
        )

        ru_path = plot_ru(
            grid_x,
            grid_y,
            ru_index_grid_masked,
            config["xmin"],
            config["xmax"],
            water_table_args,
            outdir_figs,
            phase_name,
            phase_id,
            reference_water_table_args=reference_water_table_args,
            ru_min=config.get("ru_min", 0.0),
            ru_max=config.get("ru_max", None),
        )

    if config.get("print_debug_info", False):

        plot_debug_contour(
            grid_x,
            grid_y,
            grid_sat_masked,
            config["xmin"],
            config["xmax"],
            "Saturation (%)",
            title=f"{phase_name} {phase_id} - Saturation",
            cmap="viridis",
            center=None,
        )

        plot_debug_contour(
            grid_x,
            grid_y,
            grid_pactive_masked,
            config["xmin"],
            config["xmax"],
            "Pactive [kPa]",
            title=f"{phase_name} {phase_id} - Pactive",
            cmap="coolwarm",
            center=0.0,
        )

        if process_method in {
            "effective_stress_energy",
            "all",
            "ru",
        }:

            if grid_stress is None:
                grid_stress = interpolate_field(
                    df_sp,
                    stress_field,
                    grid_x,
                    grid_y,
                )

            grid_stress_masked = mask_field(
                grid_stress,
                mask,
            )

            plot_debug_contour(
                grid_x,
                grid_y,
                grid_stress_masked,
                config["xmin"],
                config["xmax"],
                f"{stress_field} [kPa]",
                title=f"{phase_name} {phase_id} - {stress_field}",
                cmap="coolwarm",
                center=0.0,
            )

        if process_method in {"ru", "all"} and ref_df_sp is not None:

            ref_grid_pactive_masked = mask_field(
                ref_grid_pactive,
                mask,
            )

            plot_debug_contour(
                grid_x,
                grid_y,
                ref_grid_pactive_masked,
                config["xmin"],
                config["xmax"],
                "Reference phase Pactive [kPa]",
                title=f"{phase_name} {phase_id} - Reference phase Pactive",
                cmap="coolwarm",
                center=0.0,
            )

            pressure_diff_grid_masked = mask_field(
                pressure_diff_grid,
                mask,
            )

            plot_debug_contour(
                grid_x,
                grid_y,
                pressure_diff_grid_masked,
                config["xmin"],
                config["xmax"],
                "Pactive - Reference phase Pactive [kPa]",
                title=f"{phase_name} {phase_id} - Pactive diff",
                cmap="coolwarm",
                center=0.0,
            )

            plot_debug_contour(
                grid_x,
                grid_y,
                ru_index_grid_masked,
                config["xmin"],
                config["xmax"],
                "Ru index [-]",
                title=f"{phase_name} {phase_id} - Ru index",
                cmap="viridis",
                center=0.0,
                vmin=None,
                vmax=None,
            )

    return {
        "phase_name": phase_name,
        "phase_id": phase_id,
        "process_method": process_method,
        "selected_metric": selected_metric,
        "selected_value": selected_value,
        "saturation_pct": saturation_pct,
        "hydraulic_energy": hydraulic_energy_value,
        "effective_stress_energy": effective_stress_energy_value,
        "ru_index": ru_index,
        "hydraulic_debug": hydraulic_debug,
        "effective_stress_debug": effective_stress_debug,
        "ru_debug": ru_debug,
        "saturation_fig": sat_path,
        "pactive_fig": pactive_path,
        "ru_fig": ru_path,
    }


def process_all_phases(dir_raw, config, outdir_figs):
    phase_names = config.get("phases_to_plot")
    if not phase_names:
        phase_names = list_phases_from_csv(dir_raw)

    reference_phase = config.get("reference_phase")
    ref_df_sp = None
    if reference_phase:
        _, ref_sp_path = get_phase_files(dir_raw, reference_phase)
        ref_df_sp = pd.read_csv(ref_sp_path)

    results = []
    for phase_name in phase_names:
        nodes_path, sp_path = get_phase_files(dir_raw, phase_name)
        phase_result = process_phase_files(
            phase_name,
            nodes_path,
            sp_path,
            config,
            outdir_figs,
            ref_df_sp=ref_df_sp,
        )
        if phase_result is not None:
            results.append(phase_result)

    print_results(results, config)
    return results


def print_results(results, config):
    print_debug_info = config.get("print_debug_info", False)

    for phase_result in results:
        print(f"Phase: {phase_result['phase_name']}, ID {phase_result['phase_id']}")
        # print(f"  Method: {phase_result['process_method']}")

        if phase_result['process_method'] == 'saturated_area':
            print(f"  Saturated area: {phase_result['saturation_pct']:.2f}%")
        elif phase_result['process_method'] == 'hydraulic_energy':
            print(f"  Hydraulic energy: {phase_result['hydraulic_energy']:.4f}")
            if print_debug_info:
                debug = phase_result['hydraulic_debug']
                print(
                    f"    Debug: numerator={debug['negative_integral']:.2f}, denominator={debug['denominator']:.2f} "
                    f"(γw={debug['gamma_w']:.1f}, A={debug['area']:.2f})"
                )
        elif phase_result['process_method'] == 'effective_stress_energy':
            print(f"  Effective stress energy: {phase_result['effective_stress_energy']:.4f}")
            if print_debug_info:
                debug = phase_result['effective_stress_debug']
                print(
                    f"    Debug: numerator={debug['stress_integral']:.2f}, denominator={debug['denominator']:.2f} "
                    f"(stress_field={debug['stress_field']}, γw={debug['gamma_w']:.1f}, A={debug['area']:.2f})"
                )
        elif phase_result['process_method'] == 'ru':
            print(f"  Ru index: {phase_result['ru_index']:.4f}")
            if print_debug_info:
                debug = phase_result['ru_debug']
                print(
                    f"    Debug: integral={debug['ru_integral']:.2f}, integrated_area={debug['integrated_area']:.2f}, "
                    f"domain_area={debug['area']:.2f}, stress_field={debug['stress_field']}, "
                    f"stress_denominator={debug['stress_denominator']}, "
                    f"include_positive_pressure_diff={debug['include_positive_pressure_diff']}, "
                    f"ignored_area={debug['ignored_area']:.2f}"
                )
        else:  # all
            if phase_result['saturation_pct'] is not None:
                print(f"  Saturated area: {phase_result['saturation_pct']:.2f}%")
            print(f"  Hydraulic energy: {phase_result['hydraulic_energy']:.4f}")
            if print_debug_info:
                debug = phase_result['hydraulic_debug']
                print(
                    f"    Debug: integral={debug['negative_integral']:.2f}, denominator={debug['denominator']:.2f} "
                    f"(γw={debug['gamma_w']:.1f}, A={debug['area']:.2f})"
                )
            print(f"  Effective stress energy: {phase_result['effective_stress_energy']:.4f}")
            if print_debug_info:
                debug2 = phase_result['effective_stress_debug']
                print(
                    f"    Debug: integral={debug2['stress_integral']:.2f}, denominator={debug2['denominator']:.2f} "
                    f"(stress_field={debug2['stress_field']}, γw={debug2['gamma_w']:.1f}, A={debug2['area']:.2f})"
                )
            if phase_result['ru_index'] is not None:
                print(f"  Ru index: {phase_result['ru_index']:.4f}")
                if print_debug_info:
                    debug3 = phase_result['ru_debug']
                    print(
                        f"    Debug: integral={debug3['ru_integral']:.2f}, integrated_area={debug3['integrated_area']:.2f} "
                        f"(stress_field={debug3['stress_field']}, A={debug3['area']:.2f}, "
                        f"stress_denominator={debug3['stress_denominator']}, "
                        f"include_positive_pressure_diff={debug3['include_positive_pressure_diff']}, "
                        f"ignored_area={debug3['ignored_area']:.2f})"
                    )

        # print(f"  Saturation figure: {phase_result['saturation_fig']}")
        # print(f"  Pactive figure: {phase_result['pactive_fig']}")
        print()
