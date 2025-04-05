# Consensus-Based MPC for Multi-Robot Formation Control

## Introduction
This repository contains the MATLAB implementation of a **Consensus-Based Model Predictive Control (MPC)** algorithm for multi-robot formation control. The code simulates a distributed multi-robot system where robots maintain a desired formation while following a leader robot's trajectory. The implementation incorporates **control barrier functions (CBFs)** for obstacle avoidance and ensures consensus among robots to achieve desired trajectories.

The code simulates a system of three robots following a virtual leader while avoiding obstacles. Each robot operates under constraints defined by its kinematics, dynamics, and inter-robot communication network.

## System Overview
- **Leader-follower architecture:** The leader robot generates a trajectory, and the followers maintain a predefined formation relative to the leader.
- **Distributed control:** Each robot computes its control inputs locally using consensus-based MPC.
- **Obstacle avoidance:** Control barrier functions ensure safety by maintaining a minimum distance from obstacles.
- **Simulation visualization:** The code generates plots to visualize robot trajectories, velocities, formation errors, wheel torques, and other key parameters.

### Diagram
Below is an example of the multi-robot communication network:

![Multi-Robot Formation Control]

## Installation Instructions
To run the simulation and visualize the results, follow these steps:

1. **Install MATLAB**: Ensure that MATLAB is installed on your system. This code was tested on MATLAB R2023b but should work on other versions.
2. **Install CasADi**:
   - Download the CasADi toolbox from [CasADi official website](https://web.casadi.org/).
   - Add the CasADi toolbox to your MATLAB path.
3. **Clone the Repository**:
   ```
   bash
   git clone https://github.com/your-username/Consensus-MPC-Formation-Control.git
   cd Consensus-MPC-Formation-Control
   ```
4. **Run the code**:
   - Open MATLAB.
   - Navigate to the folder containing the code file, e.g., Consensus_MPC_formation_control_straight.m.
   - Run the script:
       ```run('Consensus_MPC_formation_control_straight.m')```
## How to Run the Code
1. Open the script, e.g., `Consensus_MPC_formation_control_straight.m`, in MATLAB.
2. Modify parameters if necessary:
   - `activate_formation_control`: Set to `true` to enable formation control with consensus.
   - **Obstacle Position**: Adjust `p_obs` to change the obstacle's position.
   - **Simulation time**: Modify `t0`, `tf`, and `dt` for the desired simulation duration and timestep.
3. Run the script and view the generated plots.

The script will produce the following visualizations:
- Trajectories of robots and the leader.
- Linear and angular velocities.
- Formation errors in x and y directions.
- Wheel angular velocities and torques.

## Technical Details
### System Dynamics and Kinematics
- **Kinematics**: The robots' motion is governed by differential drive kinematics, which includes states for position, orientation, and wheel velocities.
- **Dynamics**:  The motor dynamics for each wheel are modeled, including resistances, inertias, and damping coefficients.
- **Consensus Algorithm**: The robots use a consensus-based approach to minimize formation errors relative to their neighbors and the leader.

### Control Methodology
- **Model Predictive Control (MPC)**: The MPC formulation optimizes the robots' trajectories over a finite time horizon.
- **Control Barrier Functions (CBF)**: CBFs ensure that robots avoid collisions with obstacles by enforcing safety constraints.

### Parameters
The following key parameters are defined in the script:
- **Robot Parameters**
  - e.g., `R, L, J_m, m`: Physical properties of the robot and its wheels.
- **Formation COntorl**
  - e.g., `Delta`: Desired relative positions between robots in the formation.
- **Simulation Settings**
  - e.g., `tf, dt`: Final simulation time and timestep.

### Key Equations
Control inputs are derived by solving optimization problems formulated as nonlinear programs (NLPs) using CasADi. The objective function includes:
- Formation error minimization.
- Penalization of control input rates.
- Obstacle avoidance constraints.

## Known Issues and Future Improvements
### Known Issues
1. **Obstacle position**: The obstacle position is hardcoded and must be manually modified in the script.
2. **Leader trajectory**: The leader's trajectory is predefined and does not dynamically adapt to the environment.
3. **Computation time**: The simulation can be slow for large numbers of robots due to the computational complexity of solving NLPs for each robot.

### Future Improvements
1. **Dynamic obstacles**: Incorporate dynamic obstacle avoidance using real-time updates.
2. **Scalability**: Optimize the code for larger robot teams by reducing computational overhead.
3. **Real-world implementation**: Extend the simulation to work with hardware robots.

## Example Results
Sample plots generated by the simulation:
1. **Robot Trajectories**:
   ![Robot Trajectories](consensus_robot_trajectories_straight.pdf)
2. **Left Wheel Velocities**:

3. **Formation errors**:

## Repository Link
You can access the repository at: https://github.com/FHL-08/Consensus_based_formation_NMPC.
