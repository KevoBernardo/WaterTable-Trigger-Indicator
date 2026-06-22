import os
import warnings
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm
import numpy as np


def save_figure(fig, outdir_figs, file_name, dpi=300):
    os.makedirs(outdir_figs, exist_ok=True)
    fig_path = os.path.join(outdir_figs, file_name)
    fig.savefig(fig_path, dpi=dpi)
    return fig_path


def draw_water_table(ax, grid_x, grid_y, grid_psteady_masked, label, color):
    ax.contour(
        grid_x, grid_y, grid_psteady_masked,
        levels=[0],
        colors=[color],
        linewidths=1.5,
        linestyles='solid',
        zorder=10
    )
    proxy_line = Line2D([0], [0], color=color, linewidth=1.5, linestyle='solid')
    ax.legend([proxy_line], [label], loc="lower right")


def plot_saturation(grid_x, grid_y, grid_sat_masked, xmin, xmax, water_table_args, outdir_figs, phase_id, phase_name):
    fig, ax = plt.subplots(figsize=(10, 5))
    cp = ax.contourf(
        grid_x, grid_y, grid_sat_masked,
        levels=21,
        cmap="viridis",
        vmin=0, vmax=100
    )
    fig.colorbar(cp, ax=ax, label="Saturation [-] (%)")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_xlim(xmin, xmax)
    ax.axis("equal")

    if water_table_args["enabled"]:
        draw_water_table(
            ax,
            grid_x,
            grid_y,
            water_table_args["grid_psteady_masked"],
            water_table_args["label"],
            water_table_args["color"]
        )

    fig_name = f"{phase_id}_{phase_name}_saturation.png"
    fig_path = save_figure(fig, outdir_figs, fig_name)
    plt.show()
    return fig_path


def plot_debug_contour(grid_x, grid_y, field_masked, xmin, xmax, label, title=None, cmap="coolwarm", center=0.0, vmin=None, vmax=None):
    if np.all(np.isnan(field_masked)) or not np.any(np.isfinite(field_masked)):
        warnings.warn(f"Skipping debug contour for '{label}' because the field contains no finite values.")
        return

    fig, ax = plt.subplots(figsize=(10, 5))
    finite_masked = np.where(np.isfinite(field_masked), field_masked, np.nan)

    # Definir límites del colormap
    if vmin is None:
        vmin = np.nanmin(finite_masked)
    if vmax is None:
        vmax = np.nanmax(finite_masked)

    if np.isclose(vmin, vmax):
        warnings.warn(f"Debug field '{label}' has constant values; contour will use a single level.")
        levels = [vmin]
        norm = None
    else:
        # CLAVE: recortar el rango de color a [vmin, vmax]
        levels = np.linspace(vmin, vmax, 21)

        if center is not None and vmin < center < vmax:
            norm = TwoSlopeNorm(vmin=vmin, vcenter=center, vmax=vmax)
        else:
            norm = None

    cp = ax.contourf(
        grid_x,
        grid_y,
        finite_masked,
        levels=levels,
        cmap=cmap,
        norm=norm,
        extend="both",
    )

    fig.colorbar(cp, ax=ax, label=label)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_xlim(xmin, xmax)
    ax.axis("equal")

    if title:
        ax.set_title(title)

    plt.show()


def plot_pactive(grid_x, grid_y, grid_pactive_masked, xmin, xmax, water_table_args, outdir_figs, phase_name, phase_id):
    fig, ax = plt.subplots(figsize=(10, 5))
    vmin = np.nanmin(grid_pactive_masked)
    vmax = np.nanmax(grid_pactive_masked)
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)

    cp = ax.contourf(
        grid_x, grid_y, grid_pactive_masked,
        levels=21,
        cmap="coolwarm",
        norm=norm
    )
    fig.colorbar(cp, ax=ax, label="Pactive [kPa]")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_xlim(xmin, xmax)
    ax.axis("equal")

    if water_table_args["enabled"]:
        draw_water_table(
            ax,
            grid_x,
            grid_y,
            water_table_args["grid_psteady_masked"],
            water_table_args["label"],
            water_table_args["color"]
        )

    fig_name = f"{phase_name}_{phase_id}_pactive.png"
    fig_path = save_figure(fig, outdir_figs, fig_name)
    # plt.close(fig)
    plt.show()
    return fig_path
