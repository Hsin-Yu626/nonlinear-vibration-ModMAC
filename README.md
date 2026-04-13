# Nonlinear Vibration Analysis: Closely Spaced Modes with Cubic Stiffness

This repository contains MATLAB implementations developed for a Master of Engineering capstone at the University of Cincinnati (2026), extending the work of Liao [1] on a symmetric 5-DOF lumped-parameter system with hardening cubic stiffness nonlinearity.

## Research Summary

When natural frequencies are closely spaced and nonlinear stiffness is present, linear superposition fails. This project investigates:

- **Nonlinear frequency response** via the Harmonic Balance Method (HBM) and shooting method with arc-length continuation
- **Stability classification** via Floquet analysis (fold / Neimark–Sacker / period-doubling bifurcations)
- **Modal participation tracking** via a modified Modal Assurance Criterion (ModMAC)
- **Beating phenomenon** confirming simultaneous excitation of closely spaced modes under nonlinear forcing

Key finding: at 6 N excitation, seven coexisting periodic solutions exist at 45.8 Hz (consistent with Bézout's theorem), and the beat frequency (0.36 Hz) is smaller than the linear modal spacing (0.487 Hz), confirming that hardening nonlinearity pulls the modes closer together.

---

## Repository Structure

```
├── LiaoCase_ShootingV2.m     % Shooting method + Floquet stability (up/down sweep)
├── ModMAC_PostProcess.m      % ModMAC post-processing and visualization
└── README.md
```

---

## Dependencies

This code requires the **NLvib toolbox** (Krack & Gross, Springer 2019):

> Download: https://www.ifb.uni-stuttgart.de/en/research/software/nlvib/

After downloading, set the path in `LiaoCase_ShootingV2.m`:
```matlab
addpath('path/to/NLvib/SRC');
addpath('path/to/NLvib/SRC/MechanicalSystems');
```

**MATLAB version:** Tested on R2023b. Requires core MATLAB only (no additional toolboxes beyond those bundled with NLvib).

---

## Usage

### 1. Run the shooting method

```matlab
% In LiaoCase_ShootingV2.m, set parameters at the top:
F_amp = 6;      % excitation amplitude [N]
ds    = 0.0005; % arc-length step size

% Then run — results are saved to workspace and a .mat file
run('LiaoCase_ShootingV2.m')
```

This computes:
- `X_up`, `X_dn` — solution branches (up/down sweep)
- `a_rms_up`, `a_rms_dn` — RMS response amplitudes
- `stable_up`, `stable_dn` — Floquet stability flags
- `BPs_up`, `BPs_dn` — bifurcation point structs (type, frequency, amplitude)

### 2. Run ModMAC post-processing

```matlab
% Requires oscillator, X_up, X_dn, stable_up, stable_dn in workspace
run('ModMAC_PostProcess.m')
```

Generates:
- **Figure 1:** FRF coloured by dominant ModMAC mode + heatmap
- **Figure 2:** ModMAC curves per mode (solid = stable, dashed = unstable)

---

## System Description

The 5-DOF symmetric chain model (from Liao [1]) has masses `m = [150, 5, 5, 5, 5]` kg with cubic stiffness nonlinearity (`γ = 9×10¹⁰ N/m³`) at nodes 2 and 3. Modes 3 and 4 are closely spaced with a linear frequency gap of Δf = 0.487 Hz.

| Mode | Frequency (Hz) | Description |
|------|---------------|-------------|
| 1 | 16.900 | In-phase translation |
| 2 | 17.037 | Out-of-phase translation |
| 3 | 44.116 | In-phase bending |
| 4 | 44.603 | Out-of-phase bending |
| 5 | 51.299 | Higher-order bending |

---

## ModMAC Definition

The modified Modal Assurance Criterion tracks the correlation between the operational deflection shape (ODS) and each linear mode shape:

$$\text{ModMAC}_{ij} = \frac{|\{\Phi_B\}^H \mathbf{M} \{\phi_j\}|^2}{(\{\phi_j\}^T \mathbf{M} \{\phi_j\})(\{\Phi_B\}^H \mathbf{M} \{\Phi_B\})}$$

where `Φ_B = q(0) + i·q̇(0)` is the complex ODS extracted from the shooting solution, and `φ_j` is the j-th mass-normalized linear mode shape.

---

## References

[1] W.-T. Liao, "Vibration Characteristics of Symmetric Structures with Nonlinear Elements," M.S. thesis, Dept. Mech. Eng., National Chung Hsing University, Taiwan, August 2025.

[2] M. Krack and J. Gross, *Harmonic Balance for Nonlinear Vibration Problems*, Springer, 2019.

[3] M. Peeters et al., "Nonlinear normal modes, Part II: Toward a practical computation using numerical continuation techniques," *Mechanical Systems and Signal Processing*, vol. 23, 2009.

---

## Author

**Hsin-Yu (George) Tseng**  
M.Eng. Mechanical Engineering, University of Cincinnati  
M.S. Mechanical Engineering, National Chung Hsing University
