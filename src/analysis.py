import numpy as np
from scipy.interpolate import griddata
from scipy.spatial import Voronoi
from shapely.geometry import Polygon
from shapely.ops import unary_union
from matplotlib.path import Path

def build_domain(df_nodes, num_nodes_per_elem):
    element_polys = []
    for elem_id, group in df_nodes.groupby("element_id"):
        if len(group) != num_nodes_per_elem:
            continue

        group = group.sort_values("local_node")
        coords = list(zip(group["x_node"], group["y_node"]))

        poly = Polygon(coords)
        if not poly.is_valid:
            poly = poly.convex_hull

        if poly.area > 0:
            element_polys.append(poly)

    return unary_union(element_polys)


def _classified_areas(domain, df_sp, classifier):
    points = df_sp[["x_sp", "y_sp"]].values
    vor = Voronoi(points)

    area_matching = 0.0
    area_not_matching = 0.0
    vor_polys = []

    for i, region_index in enumerate(vor.point_region):
        vertices = vor.regions[region_index]
        if -1 in vertices or len(vertices) == 0:
            continue

        region_coords = [vor.vertices[v] for v in vertices]
        cell = Polygon(region_coords)
        if not cell.is_valid or cell.area == 0:
            continue

        clipped = cell.intersection(domain)
        if clipped.is_empty:
            continue

        vor_polys.append(clipped)
        area = clipped.area
        if classifier(df_sp.iloc[i]):
            area_matching += area
        else:
            area_not_matching += area

    return area_matching, area_not_matching


def suction_areas(domain, df_sp):
    return _classified_areas(domain, df_sp, lambda row: row["pactive"] <= 0)


def saturated_areas(domain, df_sp):
    return _classified_areas(domain, df_sp, lambda row: row["saturation"] == 100)


def hydraulic_energy(domain, df_sp, gamma_w=10.0, H=None):
    points = df_sp[["x_sp", "y_sp"]].values
    vor = Voronoi(points)

    if H is None:
        H = df_sp["y_sp"].max() - df_sp["y_sp"].min()

    if H <= 0:
        raise ValueError("H must be positive for hydraulic energy calculation")

    pwp_negative_integral = 0.0
    for i, region_index in enumerate(vor.point_region):
        vertices = vor.regions[region_index]
        if -1 in vertices or len(vertices) == 0:
            continue

        region_coords = [vor.vertices[v] for v in vertices]
        cell = Polygon(region_coords)
        if not cell.is_valid or cell.area == 0:
            continue

        clipped = cell.intersection(domain)
        if clipped.is_empty:
            continue

        pactive = df_sp.iloc[i]["pactive"]
        if pactive < 0:
            pwp_negative_integral += pactive * clipped.area

    area = domain.area
    if area <= 0:
        raise ValueError("Domain area must be positive for hydraulic energy calculation")

    denominator = gamma_w * H * area
    energy = -1 * pwp_negative_integral / denominator
    debug = {
        "negative_integral": pwp_negative_integral,
        "gamma_w": gamma_w,
        "H": H,
        "area": area,
        "denominator": denominator,
        "hydraulic_energy": energy,
    }
    return energy, debug


def hydraulic_stress_state(domain, df_sp, stress_field="SigyyE", gamma_w=10.0):
    points = df_sp[["x_sp", "y_sp"]].values
    vor = Voronoi(points)

    if stress_field not in {"SigyyE", "MeanEffStress"}:
        raise ValueError("stress_field must be 'SigyyE' or 'MeanEffStress'")
    if stress_field not in df_sp.columns:
        raise ValueError(f"Stress field '{stress_field}' not found in stresspoints")

    stress_integral = 0.0
    for i, region_index in enumerate(vor.point_region):
        vertices = vor.regions[region_index]
        if -1 in vertices or len(vertices) == 0:
            continue

        region_coords = [vor.vertices[v] for v in vertices]
        cell = Polygon(region_coords)
        if not cell.is_valid or cell.area == 0:
            continue

        clipped = cell.intersection(domain)
        if clipped.is_empty:
            continue

        stress_value = df_sp.iloc[i][stress_field]
        stress_integral += stress_value * clipped.area

    area = domain.area
    if area <= 0:
        raise ValueError("Domain area must be positive for effective stress energy calculation")

    denominator = gamma_w * area # TODO The denominator should be the integral of SigyyE over the domain, not gamma_w * area. This is a placeholder until the correct calculation is implemented.
    energy = -1 * stress_integral / denominator
    debug = {
        "stress_integral": stress_integral,
        "stress_field": stress_field,
        "gamma_w": gamma_w,
        "area": area,
        "denominator": denominator,
        "hydraulic_stress_state_index": energy,
    }
    return energy, debug

def ru_method(domain, df_sp, df_ref_sp, stress_field="SigyyE"):
    import numpy as np
    from scipy.spatial import Voronoi, KDTree
    from shapely.geometry import Polygon

    # ---------------- CHECKS ----------------
    if "pactive" not in df_ref_sp.columns:
        raise ValueError("Reference stresspoints must contain pactive values")

    if stress_field not in {"SigyyE", "MeanEffStress"}:
        raise ValueError("stress_field must be 'SigyyE' or 'MeanEffStress'")

    if stress_field not in df_sp.columns:
        raise ValueError(f"Stress field '{stress_field}' not found in stresspoints")

    # ---------------- ALIGN BY COORDS ----------------
    df_current = df_sp.copy()
    df_ref = df_ref_sp.copy()

    df_current["x_r"] = df_current["x_sp"].round(8)
    df_current["y_r"] = df_current["y_sp"].round(8)

    df_ref["x_r"] = df_ref["x_sp"].round(8)
    df_ref["y_r"] = df_ref["y_sp"].round(8)

    df_merged = df_current.merge(
        df_ref[["x_r", "y_r", "pactive"]],
        on=["x_r", "y_r"],
        how="inner",
        suffixes=("", "_ref")
    )

    print("\n[DEBUG] Matching summary:")
    print(" current   :", len(df_current))
    print(" reference :", len(df_ref))
    print(" matched   :", len(df_merged))

    if len(df_merged) != len(df_current):
        raise ValueError("Not all stresspoints were matched between phases")

    # ---------------- DISTANCE CHECK (KDTree) ----------------
    points = df_merged[["x_sp", "y_sp"]].values
    ref_points = df_ref[["x_r", "y_r"]].values

    tree = KDTree(ref_points)
    dist, ind = tree.query(points)

    print("\n[DEBUG] Distance statistics (current vs nearest ref):")
    print(" min  :", dist.min())
    print(" mean :", dist.mean())
    print(" max  :", dist.max())

    np.random.seed(0)
    idx = np.random.choice(len(points), min(5, len(points)), replace=False)

    print("\n[DEBUG] Sample distance check:")
    for i in idx:
        print(f"i={i}")
        print(" current :", points[i])
        print(" nearest :", ref_points[ind[i]])
        print(" dist    :", dist[i])
        print()

    # ---------------- VORONOI ----------------
    df_merged = df_merged.reset_index(drop=True).sort_values(by=["x_sp", "y_sp"])
    points = df_merged[["x_sp", "y_sp"]].values
    vor = Voronoi(points)

    ru_integral = 0.0

    for i, region_index in enumerate(vor.point_region):

        vertices = vor.regions[region_index]

        if -1 in vertices or len(vertices) == 0:
            continue

        region_coords = [vor.vertices[v] for v in vertices]
        cell = Polygon(region_coords)

        if not cell.is_valid or cell.area == 0:
            continue

        clipped = cell.intersection(domain)

        if clipped.is_empty:
            continue

        stress_value = df_merged.iloc[i][stress_field]

        if stress_value == 0:
            continue

        pressure_diff = (
            df_merged.iloc[i]["pactive"]
            - df_merged.iloc[i]["pactive_ref"]
        )

        ratio = pressure_diff / stress_value

        ru_integral += ratio * clipped.area

    debug = {
        "ru_integral": ru_integral,
        "stress_field": stress_field,
        "area": domain.area,
        "n_points": len(df_merged),
        "max_point_distance": float(dist.max()),
        "mean_point_distance": float(dist.mean()),
    }

    return ru_integral, debug


# def ru_method(domain, df_sp, df_ref_sp, stress_field="SigyyE"):
#     points = df_sp[["x_sp", "y_sp"]].values
#     vor = Voronoi(points)

#     if "pactive" not in df_ref_sp.columns:
#         raise ValueError("Reference stresspoints must contain pactive values for Ru calculation")

#     if stress_field not in {"SigyyE", "MeanEffStress"}:
#         raise ValueError("stress_field must be 'SigyyE' or 'MeanEffStress'")

#     if stress_field not in df_sp.columns:
#         raise ValueError(f"Stress field '{stress_field}' not found in stresspoints")

#     ref_points = df_ref_sp[["x_sp", "y_sp"]].values
#     ref_pactive_values = df_ref_sp["pactive"].values

#     ref_pactive = griddata(
#         points=ref_points,
#         values=ref_pactive_values,
#         xi=points,
#         method='linear'
#     )

#     ref_pactive_nearest = griddata(
#         points=ref_points,
#         values=ref_pactive_values,
#         xi=points,
#         method='nearest'
#     )

#     ref_pactive = np.where(np.isnan(ref_pactive), ref_pactive_nearest, ref_pactive)

#     if np.any(np.isnan(ref_pactive)):
#         raise ValueError("Reference phase pactive interpolation failed for some stresspoints")

#     ru_integral = 0.0

#     for i, region_index in enumerate(vor.point_region):
#         vertices = vor.regions[region_index]

#         if -1 in vertices or len(vertices) == 0:
#             continue

#         region_coords = [vor.vertices[v] for v in vertices]
#         cell = Polygon(region_coords)

#         if not cell.is_valid or cell.area == 0:
#             continue

#         clipped = cell.intersection(domain)

#         if clipped.is_empty:
#             continue

#         stress_value = df_sp.iloc[i][stress_field]
#         pressure_diff = df_sp.iloc[i]["pactive"] - ref_pactive[i]

#         if stress_value == 0:
#             continue  # evita singularidad

#         ratio = pressure_diff / stress_value

#         ru_integral += ratio * clipped.area

#     debug = {
#         "ru_integral": ru_integral,
#         "stress_field": stress_field,
#         "area": domain.area,
#     }

#     return ru_integral, debug


def make_grid(xmin, xmax, ymin, ymax, resolution=500):
    return np.mgrid[xmin:xmax:resolution*1j, ymin:ymax:resolution*1j]


def mask_domain(domain, grid_x, grid_y):
    grid_points = np.vstack((grid_x.ravel(), grid_y.ravel())).T
    mask = np.zeros(grid_x.shape, dtype=bool)

    for geom in getattr(domain, "geoms", [domain]):
        path = Path(np.array(geom.exterior.coords))
        mask |= path.contains_points(grid_points).reshape(grid_x.shape)

    return mask


def interpolate_field(df_sp, field_name, grid_x, grid_y, method='linear'):
    return griddata(
        points=df_sp[["x_sp", "y_sp"]].values,
        values=df_sp[field_name].values,
        xi=(grid_x, grid_y),
        method=method
    )


def mask_field(field_grid, mask):
    return np.where(mask, field_grid, np.nan)
