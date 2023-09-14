# NMSM_Forward_Dynamics
Neuromusculoskeletal Forward Dynamics Framework with Continuous Muscle Wrapping

## Simulations

A simple but representative neuromusculoskeletal model (NMSM)  was simulated under two conditions: (I) without an external wrench applied to the upper-limb and (II) with a time-varying external wrench applied on the wrist due to the interaction with a robot.

For condition I, the results obtained between the proposed framework (using the methodology presented in Section III) and OpenSim (Muscle Analysis Tool) are compared.

- Simulation_I.m: A forward upper-limb NMSM is driven only by excitation signals with no external wrench applied
- Simulation_II.m: A forward upper-limb NMSM is driven only by excitation signals with a time-varying external wrench applied on the wrist

> **Note:** The file OpenSim_Sim/arm22_test.cpp was used to generate the muscle analysis in OpenSim.
> The NMSM can be found in spatial_v2/models/arm22_model.m for use in Matlab, and in OpenSim_Sim/arm22_dG.osim for use in OpenSim.

**Equation (13) should read as follows**
$$\mathbf{R}_{w} &= \mathbf{l}^{{(i+1)}^{T}}\begin{bmatrix}
-{^{0}\mathbf{E}_{\lambda_{P^{*}}}}\left[{^{(\lambda_{Q})}\bm{r}_{Q^{i}}}\times\right] & {^{0}\mathbf{E}_{\lambda_{P^{*}}}}
\end{bmatrix}{^{(\lambda_{P^{*}})}\mathbf{J}}_{\lambda_{P^{*}}} \\
&\quad
- \mathbf{l}^{{(i)}^{T}}\begin{bmatrix}
-{^{0}\mathbf{E}_{\lambda_{Q^{*}}}}\left[{^{(\lambda_{P})}\bm{r}_{P^{i}}}\times\right] & {^{0}\mathbf{E}_{\lambda_{Q^{*}}}}
\end{bmatrix}{^{(\lambda_{Q^{*}})}\mathbf{J}_{\lambda_{Q^{*}}}} \\ &\quad
+ \mathbf{l}^{{(i)}^{T}}\begin{bmatrix}
-{^{0}\mathbf{E}_{\lambda_{S}}}\left[{^{(\lambda_{S})}\bm{r}_{S^{i}}}\times\right] & {^{0}\mathbf{E}_{\lambda_{S}}}
\end{bmatrix}{^{(\lambda_{S})}\mathbf{J}}_{\lambda_{S}} \\ &\quad
- \mathbf{l}^{{(i+1)}^{T}}\begin{bmatrix}
-{^{0}\mathbf{E}_{\lambda_{S}}}\left[{^{(\lambda_{S})}\bm{r}_{S^{i}}}\times\right] & {^{0}\mathbf{E}_{\lambda_{S}}}
\end{bmatrix}{^{(\lambda_{S})}\mathbf{J}}_{\lambda_{S}}$$
