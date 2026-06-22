# WaterTable-Trigger-Indicator

A Python tool to compute spatially integrated indicators for water table triggering assessment in Tailings Storage Facilities (TSFs), using numerical outputs from PLAXIS stress-point data.

---

## Overview

This repository implements a family of spatial indicators designed to quantify hydraulic and mechanical triggering conditions associated with water table evolution in geotechnical systems.

All indicators are computed through post-processing of finite element results (e.g., PLAXIS). The spatial integration is performed using a Voronoi discretization of the stress-point cloud, where each numerical point represents a physically consistent area contribution within a user-defined domain (e.g., a vulnerable region of a tailings dam).

---

## General Methodology

All indicators share the same spatial framework:

\[
I = \int_A f(\mathbf{x}) \, dA
\]

where:

- \(A\): spatial domain of interest (vulnerable region)
- \(f(\mathbf{x})\): field function derived from PLAXIS outputs
- Discretization: Voronoi tessellation of stress points

Each stress point contributes proportionally to its associated Voronoi cell area.

---

## 1. $\mathbf{r_u}$ Indicator

\[
I_{WT} = \int_A \frac{\Delta p}{\sigma} \, dA
\]

where:

- \(\Delta p = p_{\text{active}} - p_{\text{reference}}\)
- \(\sigma\): selected stress field (e.g. SigyyE, MeanEffStress)

### Interpretation

Measures spatially accumulated pore pressure changes normalized by stress state, highlighting regions approaching potential hydraulic triggering conditions.

---

## 2. Suction Areas Indicator

### Definition

Partitions the domain based on pore pressure sign:

- Suction (negative pressure): \(p_{\text{active}} < 0\)
- Non-suction (positive or zero pressure): \(p_{\text{active}} \ge 0\)

### Outputs

- \(A_{suction}\): total area under suction
- \(A_{non-suction}\): remaining area

### Interpretation

Quantifies the spatial extent of matric suction conditions, relevant for partially saturated soil behavior and stability.

---

## 3. Saturation Areas Indicator

### Definition

Classifies the domain based on full saturation condition:

- Fully saturated: \(p_{\text{active}} \le 0\) (100% saturation assumption)
- Not fully saturated: \(p_{\text{active}} > 0\)

### Outputs

- \(A_{sat}\): fully saturated area
- \(A_{unsat}\): non-saturated or partially saturated area

### Interpretation

Provides a geometric measure of saturation extent within the domain, useful for tracking water table advancement or retreat.

---

## 4. Hydraulic Energy Indicator

\[
E_h = - \frac{\int_A p_{\text{active}}^{-} \, dA}{\gamma_w H A}
\]

where:

- \(p_{\text{active}}^{-}\): negative pore pressure contribution
- \(\gamma_w\): unit weight of water
- \(H\): characteristic hydraulic head
- \(A\): domain area

---

## 5. Hydraulic Stress-State Indicator

\[
I_{\sigma} = - \frac{\int_A \sigma \, dA}{\gamma_w A}
\]

*(current implementation uses a simplified denominator; refinement pending)*

---

## Features

- Reads PLAXIS stress-point outputs
- Handles multiple simulation phases
- Matches stress points by spatial coordinates
- Builds Voronoi tessellation for spatial integration
- Computes multiple hydraulic and mechanical indicators
- Includes debugging tools for consistency checks

---

## Applications

- Tailings dam stability assessment
- Water table triggering analysis
- Saturation and suction mapping
- Liquefaction-related screening metrics
- Geotechnical numerical model post-processing

---

## Requirements

- numpy
- scipy
- shapely
- matplotlib (optional, debugging)

---

## Usage Example

```python
from ru_module import ru_method

Iwt, debug = ru_method(domain, df_sp, df_ref_sp, stress_field="SigyyE")
print(Iwt)