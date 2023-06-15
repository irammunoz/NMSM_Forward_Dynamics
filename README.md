# NMSM_Forward_Dynamics
Neuromusculoskeletal Forward Dynamics Framework with Continuous Muscle Wrapping

## Simulations

A simple but representative neuromusculoskeletal model (NMSM)  was simulated under two conditions: (I) without an external wrench applied to the upper-limb and (II) with a time-varying external wrench applied on the wrist due to the interaction with a robot.

For condition I, the results obtained between the proposed framework (using the methodology presented in Section III) and OpenSim (Muscle Analysis Tool) are compared.

- Simulation_I.m: A forward upper-limb NMSM is driven only by excitation signals with no external wrench applied
- Simulation_II.m: A forward upper-limb NMSM is driven only by excitation signals with a time-varying external wrench applied on the wrist

> **Note:** The file OpenSim_Sim/arm22_test.cpp was used to generate the muscle analysis in OpenSim.
