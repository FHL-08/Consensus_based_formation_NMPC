clear
close all
clc
import casadi.*

t0 = 0; % start time
tf = 40; % end time
dt = 0.1; % step size
k = t0:dt:tf-dt;
ACTIVATE_CONSENSUS = 'active';

% NMPC finite horizon lengths
N = 10; % trajectory following horizon length
T = length(k) + N;
t = t0:dt:(T-1)*dt;

U_0 = 10; % Source Voltage
R = 0.5; % Armature Resistance
R_z = 2e-4; % Internal Resistance of voltage source
J_m = 1e-6; % Inertia of Motor Shaft
L = 100e-3; % Armature Inductance
N_g = 1; % Gearbox Transmission Ratio
l_L = 0.125; % Distance between left wheel and centre of wheel axle
l_R = l_L; % Distance between right wheel and centre of wheel axle
r = 0.05; % wheel d_safe
r_G = r/N_g;
m_p = 4.385;
m_w = 0.0575;
m = m_p + 2*m_w; % Mass of robot
l_T = 0.11; % distance between centre of wheel shaft and centre of gravity
J_zp = 0.05; % Intertia of the platform about the centroid
J_zw = 0.002; % Intertia of the wheels about the centroid
J_B = J_zp + 2*J_zw + m_w*(2*l_T^2 + l_L^2 + l_R^2); % Intertia about the centre of rotation
K_t = 5.305e-2; % Torque Constant
K_e = K_t; % Torque Constant
k_r = 4e-5; % Damping coefficient of Motor Shaft
k_v = 0.25; % Damping coefficient of linear motion
k_w = 0.15; % Dampin coefficient of rotational motion

alpha = (r_G^2)/(l_L + l_R);
scaled_radius = (alpha/r_G);
a_L = k_r + k_v*l_R*alpha;
a_R = k_r + k_v*l_L*alpha;
c_L = k_r*l_L + k_w*alpha;
c_R = k_r*l_R + k_w*alpha;

b_L = J_m + m*l_R*alpha;
b_R = J_m + m*l_L*alpha;
d_L = J_m*l_L + J_B*alpha;
d_R = J_m*l_R + J_B*alpha;

RL = -(1/L)*[R+R_z R_z;
             R_z R+R_z];
KL = -(K_e/L)*eye(2);
ac_mat = [a_L a_R;
          -c_L c_R];
bd_mat = [d_R -b_R;
          d_L b_L];
LR = [1 1;
      -l_L l_R];
beta = 1/det(bd_mat);

%%% Matrix Definition for low-level controller
A = [RL KL;
     (K_t*beta)*bd_mat*LR -beta*bd_mat*ac_mat];
B = (U_0/L)*eye(4, 2);
C = r_G*[zeros(2) (LR'\eye(2))];
electr_sys = ss(A, B, C, zeros(2));

%%% System Parameters
x = MX.sym('x');
y = MX.sym('y');
theta = MX.sym('theta');

omega = MX.sym('omega');
vel = MX.sym('vel', 2);
Phi = MX.sym('Phi', 2);
tau = MX.sym('tau', 2);

dim_vel = length(Phi);
dim_tau = length(tau);

%%% Robot Kinematics part 1
q = [x; y; theta]; % states
dim_state = length(q); % dimension of state vector

R_q = [cos(q(3)) -l_T*sin(q(3));
       sin(q(3)) l_T*cos(q(3));
       0 1];

R_q_func = Function('R_q_func', {q}, {R_q});
dq = R_q*vel;
dq_func = Function('dq_func', {q, vel}, {dq});

% Nonlinear solver. Produces q(k + 1) given q(k) and u
q_ode = struct('x', q, 'p', vel, 'ode', dq);
q_solver = integrator('q_solver', 'cvodes', q_ode, t0, dt);

% Defines the ith standard basis vector e
function e = std_basis(n, i)
    e = zeros(n,1);
    e(i) = 1;
end

%%% Robot Kinematics part 2
M_bD = [l_R l_L; l_T -l_T];
bar_M_bD = [l_T l_R; l_T -l_L];
M_permute = [0 1; 1 0];
unit_heading_vec = [cos(q(3)); sin(q(3))];
S_y_q = scaled_radius*[M_bD (M_permute')*M_bD*M_permute]*kron(unit_heading_vec, eye(2));
inv_S_y_q = (1/(r*l_T))*[bar_M_bD (M_permute')*bar_M_bD*M_permute]*kron(unit_heading_vec, eye(2));
S_q = [S_y_q; scaled_radius*[1 -1]];

% convert wheel velocities to robot cartesian velocities
S_y_q_func = Function('S_y_q_func', {theta}, {S_y_q});
% inverse transformation: converts cartesian velocities to wheel velocities
inv_S_y_q_func = Function('inv_S_y_q_func', {theta}, {inv_S_y_q});

skew_mat = [0 -1;
            1  0];

% angular velocity
omega_q = scaled_radius*[1 -1]*Phi;

% skew symmetric matrix S(w)
S_w = omega_q*skew_mat;
S_omega = omega*skew_mat;

dS_y_q = S_w*S_y_q;
S_w_func = Function('S_w_func', {Phi, theta}, {S_w});
S_omega_func = Function('S_omega_func', {omega}, {S_omega});

%%% Defintion of unprojected system matrices
% M_q = diag([m; m; J_B]);
% C_q = zeros(3, 1);
% B_q = (N_g/r_G)*[cos(q(3)), cos(q(3));
%                  sin(q(3)), sin(q(3));
%                  l_R      , -l_L];

det_S_y = -(scaled_radius^2)*(l_R + l_L)*(l_T);
bar_M_q = m*(scaled_radius^2)*[(l_R^2 + l_T^2), (l_L*l_R - l_T^2);
                             (l_R*l_L - l_T^2), (l_L^2 + l_T^2)]...
          + J_B*(scaled_radius^2)*[1, -1; -1, 1];
inv_bar_M_q = bar_M_q\eye(2);
bar_C_q = m*det_S_y*S_w;
bar_B_q_coeff = N_g;
bar_B_q = bar_B_q_coeff*eye(dim_tau);

dphi = inv_bar_M_q*(-bar_C_q*Phi + bar_B_q*tau);
dphi_func = Function('dphi_func', {Phi, tau}, {dphi});

x_state = [q; Phi];
dx_state = [S_q*Phi; dphi];
global_robot_acceleration = dS_y_q*Phi + S_y_q*dphi;

% Nonlinear solver. Produces x_state(k + 1) given x_state(k) and tau
state_ode = struct('x', x_state, 'p', tau, 'ode', dx_state);
state_solver = integrator('state_solver', 'cvodes', state_ode, t0, dt);
dx_state_func = Function('dx_state_func', {x_state, tau}, {dx_state});
global_robot_acceleration_func = Function('global_robot_acceleration_func',...
                                            {x_state, tau}, {global_robot_acceleration});

%%% Define Saturation limits
volts_ub = U_0; % Volts
volts_lb = -U_0; % Volts
forward_vel = -C*(A\B)*[volts_ub; volts_ub];
wheel_vel_ub = full(inv_S_y_q_func(0)*forward_vel);
wheel_vel_lb = -wheel_vel_ub;
max_torque = 0.25;
torque_ub = max_torque*ones(2, 1);
torque_lb = -max_torque*ones(2, 1);

%%% Define Communication Network
mrs_size = 3; % number of robots
adjacency = [0 1 1;
             1 0 1;
             1 1 0]; % adjacency matrix (consensus)
vertex_degrees = sum(adjacency, 2);
leader_neighbour = [1 0 0]; % leader neighbourhood vector (consensus)
switch ACTIVATE_CONSENSUS
    case 'inactive'
        adjacency = zeros(mrs_size);  % adjacency matrix (leader-follower)
        vertex_degrees = sum(adjacency, 2); 
        leader_neighbour = ones(1, mrs_size); % leader neighbourhood (leader-follower)
    case 'active'
        % do nothing
end
degree = diag(vertex_degrees);
laplacian = degree - adjacency;
leader_adjacency = diag(leader_neighbour);
modified_laplacian = laplacian + leader_adjacency;
g_dis = modified_laplacian + adjacency;

%%% trajectory error variables
leader_state = MX.sym('leader_state', dim_state);
dleader_state = MX.sym('dleader_state', dim_state);
ddleader_state = MX.sym('ddleader_state', dim_state);
g_i_dis_var = MX.sym('g_i_dis_var');
trust_rad_var = MX.sym('trust_rad_var', 2);
basis_var = MX.sym('xi_var', mrs_size);
xi_var = MX.sym('xi_var', dim_state-1);
e_i = MX.sym('e_i', dim_state-1);
de_i = MX.sym('de_i', dim_state-1);
e_i_var = MX.sym('e_i_var', dim_state-1);
de_i_var = MX.sym('de_i_var', dim_state-1);
h_i_dis_var = MX.sym('h_i_dis_var', dim_state-1);
dh_i_dis_var = MX.sym('dh_i_dis_var', dim_state-1);
nu_var = MX.sym('nu_var', dim_tau);
P_j_var = MX.sym('P_j_var', dim_state-1, mrs_size);
dP_j_var = MX.sym('dP_j_var', dim_state-1, mrs_size);
ddP_j_var = MX.sym('ddP_j_var', dim_state-1, mrs_size);
Delta_var = MX.sym('Delta_var', dim_state-1, mrs_size);
dDelta_var = MX.sym('dDelta_var', dim_state-1, mrs_size);
ddDelta_var = MX.sym('dDelta_var', dim_state-1, mrs_size);

%%% trajectory error matrices
R_theta = [cos(q(3)), -sin(q(3));
           sin(q(3)),  cos(q(3))];
R_theta_func = Function('R_theta_func', {theta}, {R_theta});

calculate_torque = bar_B_q\(bar_M_q*nu_var + bar_C_q*Phi);
calculate_torque_func = Function('calculate_torque_func', {x_state, nu_var}, {calculate_torque});

Lambda = eye(2, 3);
copy_over_network = ones(1, mrs_size);
Leader_position_matrix = Lambda*leader_state*copy_over_network;
Leader_velocity_matrix = Lambda*dleader_state*copy_over_network;
Leader_acceleration_matrix = Lambda*ddleader_state*copy_over_network;

tracking_error = ((P_j_var - Delta_var)*modified_laplacian...
                      - Leader_position_matrix*leader_adjacency)*basis_var;
tracking_error_rate = -(dP_j_var*adjacency + dDelta_var*modified_laplacian...
                           + Leader_velocity_matrix*leader_adjacency)*basis_var...
                           + g_i_dis_var*S_y_q*Phi;
eta = [tracking_error; tracking_error_rate];
eta_func = Function('eta_func', {P_j_var, dP_j_var, Delta_var, dDelta_var,...
                                 leader_state, dleader_state,...
                                 x_state, g_i_dis_var, basis_var},...
                     {eta});

tilde_n_i = g_i_dis_var*S_y_q;
inv_tilde_n_i = inv_S_y_q/g_i_dis_var;
tilde_m_i = g_i_dis_var*dS_y_q*Phi;
tilde_l_i = -(ddP_j_var*adjacency + ddDelta_var*modified_laplacian...
            + Leader_acceleration_matrix*leader_adjacency)*basis_var;

tracking_error_accel = tilde_l_i + tilde_m_i + tilde_n_i*nu_var;
deta = [tracking_error_rate; tracking_error_accel];
deta_func = Function('deta_func', {dP_j_var, ddP_j_var, dDelta_var, ddDelta_var,...
                                   dleader_state, ddleader_state,...
                                   x_state, g_i_dis_var, basis_var, nu_var},...
                      {deta});

convert_xi_to_nu = tilde_n_i\(xi_var - tilde_l_i - tilde_m_i);
convert_xi_to_nu_func = Function('convert_xi_to_nu_func',...
                                 {ddP_j_var, ddDelta_var, ddleader_state,...
                                  x_state, g_i_dis_var, basis_var, xi_var},...
                                 {convert_xi_to_nu});

%%% formation parameters
Delta = [-1.5, -3.0, -3.0;
            0,  1.5, -1.5];

%%% leader trajectory (virtual signal generation)
start = [1.5; 1.2; 0];
vel_start = [0.25; 0; 0];
p_l = start;
dp_l = vel_start;
ddp_l = zeros(dim_state, 1);
Delta_history = Delta;
dDelta_history = zeros(dim_state-1, mrs_size);
ddDelta_history = zeros(dim_state-1, mrs_size);

% Let A be the point at the centre of the shaft connecting the two driven
% wheels.
%%% define velocities of virtual leader over trajectory
v_l = zeros(2, T);
% linear velocity of A relative ot fixed centre
v_l(1, :) = 0.25*ones(1, T);
% angular velocity of A about point remote centre
v_l(2, :) =  zeros(1, T);

q_l = start; % iniital state of virtual robot
for l=2:T
    q_l_solution = q_solver('x0', q_l, 'p', v_l(:, l-1));
    q_l = full(q_l_solution.xf);
    dq_l = R_q_func(q_l)*v_l(:, l-1);
    p_l(:, l) = q_l;
    dp_l(:, l) = full(dq_l);
    ddp_l(:, l) = zeros(dim_state, 1);
    Delta_history(:, :, l) = full(R_theta_func(q_l(3))*Delta);
    dDelta_history(:, :, l) = full(-S_omega_func(dq_l(3))*R_theta_func(q_l(3))*Delta);
    ddDelta_history(:, :, l) = full(-S_omega_func(0)*R_theta_func(q_l(3))*Delta...
                                    + (S_omega_func(dq_l(3))^2)*R_theta_func(q_l(3))*Delta);
end

%%% initial conditions
q_k = zeros(dim_state, mrs_size); % empty vector for states

% Initial Robot possitions and orientations
p1_0 = [0; 0.7; pi];
p2_0 = [-2.5; 2; pi];
p3_0 = [-1.8; -1.2; pi];

q_next_k = [p1_0 p2_0 p3_0]; % initial system state
dq_next_k = zeros(dim_state, mrs_size); % initial system velocity
robot_accelerations = zeros(dim_state-1, mrs_size);

%%% position of obstacle
tolerance = 0.25;
obstacle_radius = 0.25;
d_safe = obstacle_radius + tolerance;
obstacle_encounter_time = double(int32(10/dt));
p_obs = p_l(1:2, obstacle_encounter_time) - [0; 0.1];

%%% Definition of Discrete-time Control Barrier Function
% extract position information of robot
pos = (e_i + de_i*dt - (h_i_dis_var + dh_i_dis_var*dt))/g_i_dis_var;
pos_next = (e_i + 2*de_i*dt + xi_var*(dt^2) - (h_i_dis_var + 2*dh_i_dis_var*dt))/g_i_dis_var;
gamma_CBF = MX.sym('gamma_CBF'); 

% larger gamma values corresponds with moving closer to obstacle.
% smaller gamma values corresponds with avoiding the obstacle earlier.
% 0 < gamma_val <= 1
min_gamma = 1e-9;
max_gamma = 1 - min_gamma;
desired_gamma_CBF = 0.2;

% vector from obstacle position to robot position
r_vec = (pos - p_obs);
r_vec_next = (pos_next - p_obs);

% Definitions of CBFs
h_p = (r_vec')*r_vec;
h_p_next = (r_vec_next')*r_vec_next;

% Definitions of CBF constraints (CBF interpolation)
cbf_constraint = h_p_next - (1 - gamma_CBF)*h_p - gamma_CBF*(d_safe)^2;
cbf_constraint_func = Function('cbf_constraint_func',...
                               {e_i, de_i, h_i_dis_var, dh_i_dis_var, g_i_dis_var,...
                                gamma_CBF, xi_var}, {cbf_constraint});

% Define Optimisation Problem Bounds
Phi_e_i = inv_tilde_n_i*(de_i - dh_i_dis_var);
next_theta = theta + S_q(3, :)*Phi_e_i*dt;
next_theta_func = Function('next_theta_func', {de_i, theta, dh_i_dis_var, g_i_dis_var},...
                                              {next_theta});
z_var = 1/(g_i_dis_var*r*l_T);
X_var = 0.5*m*l_T*r^2;
Y_var = scaled_radius*r*(J_B + m*l_T^2);
R_var = norm([X_var; Y_var]);
alpha_var = atan2(Y_var, X_var);
lambda = (l_L + l_R)*J_B*r*(z_var*scaled_radius)^2;
H_1 = lambda*[-sin(2*theta) cos(2*theta); cos(2*theta) sin(2*theta)];
H_2 = -H_1;
torque_bound = 0.5*[((de_i - dh_i_dis_var)')*H_1*(de_i - dh_i_dis_var);
                    ((de_i - dh_i_dis_var)')*H_2*(de_i - dh_i_dis_var)] +...
               z_var*R_var*[cos(theta + alpha_var) sin(theta + alpha_var);
                            cos(theta - alpha_var) sin(theta - alpha_var)]*xi_var;
torque_bound_func = Function('torque_bound_func', {de_i, xi_var, theta, dh_i_dis_var, g_i_dis_var},...
                                                  {torque_bound});

%%% cost weighting matrices (tune)
CBF_weight = 100;
penalise_vel = 50;
Q = diag([15, 15, 10, 10]);
R = 25*eye(dim_state-1);
R_horizon = sparse(kron(eye(N), R));
H = [zeros(dim_state-1) eye(dim_state-1);
     zeros((dim_state-1), 2*(dim_state-1))];
G = [zeros(dim_state-1); eye(dim_state-1)];
[P,~,~] = icare(H, G, Q, R, [], [], []);

lmin_q = min(eig(Q));
gamma = min(eig(Q))/max(eig(P));
form_beta_val = 25;

% discrete consensus error model
auxiliary_sys = ss(H, G, [], []);
auxiliary_sysd = c2d(auxiliary_sys, dt);
dt_H = auxiliary_sysd.A;
dt_G = auxiliary_sysd.B;

%%% memory allocation
consensus_opt_output = repmat([repmat(zeros(dim_state-1, 1), N, 1)], 1, mrs_size);
consensus_opt_params = MX.sym('consensus_opt_params', length(consensus_opt_output));

v_i = repmat(zeros(dim_vel, 1), 1, mrs_size);
wheel_vel_i = repmat(zeros(dim_vel, 1), 1, mrs_size);
tau_i = repmat(zeros(dim_tau, 1), 1, mrs_size);
q_i_library = repmat(zeros(dim_state, length(k)), 1, 1, mrs_size);
formation_error = repmat(zeros(dim_state-1, length(k)), 1, 1, mrs_size);
vel_library = repmat(zeros(dim_vel, length(k)), 1, 1, mrs_size);
tau_library = repmat(zeros(dim_tau, length(k)), 1, 1, mrs_size);
wheeled_vel_library = repmat(zeros(dim_vel, length(k)), 1, 1, mrs_size);

%%% Define IPOPT options
opts = struct;
opts.ipopt.tol = 1e-8; % Tolerance for optimization
opts.ipopt.constr_viol_tol = 1e-8; % Constraint violation tolerance
opts.ipopt.acceptable_tol = 1e-6;
opts.ipopt.acceptable_constr_viol_tol = 1e-6;
opts.ipopt.max_iter = 120;  % Max number of iterations
opts.ipopt.print_level = 0; % Suppress output (0 for no output)
opts.print_time = false;
opts.ipopt.sb = 'yes'; % Suppress the IPOPT startup banner

p_diffeomorphism_k_init = MX.sym('p_diffeomorphism_k_init', 4);
p_h_i_dis = MX.sym('p_h_i_dis', 2);
p_dh_i_dis = MX.sym('p_dh_i_dis', 2);
p_angle_init = MX.sym('p_angle_init');
p_g_i_dis = MX.sym('p_g_i_dis');
p_relative_average_velocity = MX.sym('p_relative_average_velocity', 2);

% Initialise symbolic states and cost using the parameters
p_diffeomorphism_k_next = p_diffeomorphism_k_init;
h_i_dis_next = p_h_i_dis;
angle = p_angle_init;
consensus_stage_cost = (p_diffeomorphism_k_next')*Q*p_diffeomorphism_k_next;

formation_lb = MX(N*(dim_state-1), 1);
formation_ub = MX(N*(dim_state-1), 1);
CBF_constraint_vec = MX(N*(dim_state-1), 1);

for z=1:N
    opt_var_z = consensus_opt_params((z-1)*(dim_state-1)+1:z*(dim_state-1));
    xi_bound = torque_bound_func(p_diffeomorphism_k_next(3:4), opt_var_z,...
                                 angle, p_dh_i_dis, p_g_i_dis);
    CBF_constraint_z = cbf_constraint_func(p_diffeomorphism_k_next(1:2),...
                                           p_diffeomorphism_k_next(3:4),...
                                           h_i_dis_next, p_dh_i_dis,...
                                           p_g_i_dis, desired_gamma_CBF,...
                                           opt_var_z);

    formation_lb((z-1)*(dim_state-1)+1:z*(dim_state-1)) = xi_bound - torque_lb;
    formation_ub((z-1)*(dim_state-1)+1:z*(dim_state-1)) = torque_ub - xi_bound;
    CBF_constraint_vec((z-1)*(dim_state-1)+1:z*(dim_state-1)) = CBF_constraint_z;   

    angle = next_theta_func(p_diffeomorphism_k_next(3:4), angle, p_dh_i_dis, p_g_i_dis);
    p_diffeomorphism_k_next = dt_H*p_diffeomorphism_k_next + dt_G*opt_var_z;
    h_i_dis_next = h_i_dis_next + p_dh_i_dis*dt;
    if (z == N)
        consensus_stage_cost = consensus_stage_cost...
                          + form_beta_val*(p_diffeomorphism_k_next')*P...
                            *p_diffeomorphism_k_next;
    else
        consensus_stage_cost = consensus_stage_cost + (p_diffeomorphism_k_next')*Q...
                            *p_diffeomorphism_k_next...
                          + penalise_vel*(((p_diffeomorphism_k_next(3:4) - p_dh_i_dis)/p_g_i_dis...
                                - p_relative_average_velocity)')*...
                                ((p_diffeomorphism_k_next(3:4) - p_dh_i_dis)/p_g_i_dis...
                                - p_relative_average_velocity);
    end
end

%%% Define NLP problem
formation_sys_objective = consensus_stage_cost...
                +  0.5*(consensus_opt_params(1:N*(dim_state-1))')*R_horizon...
                    *consensus_opt_params(1:N*(dim_state-1));
% 'g' vector is formulated so that all constraints are 'g >= 0'
formation_sys_constraint = [formation_lb;
                            formation_ub;
                            CBF_constraint_vec];

%%% optimisation parameter bounds
g_vec = [formation_lb; formation_ub; CBF_constraint_vec];
formation_lbg = zeros(size(g_vec, 1), 1);
formation_ubg = inf(size(g_vec, 1), 1);

% Assemble all symbolic parameters into a single vector
all_params_sym = [p_diffeomorphism_k_init; p_h_i_dis; p_dh_i_dis; p_angle_init; p_g_i_dis;...
                  p_relative_average_velocity];

% Create the NLP solver object ONCE
consensus_nlp = struct('x', consensus_opt_params, 'p', all_params_sym, ...
                       'f', formation_sys_objective, 'g', formation_sys_constraint);
consensus_nlpsolve = nlpsol('solver', 'ipopt', consensus_nlp, opts);

%% Simulate Distributed Multirobot System Control
for j=1:length(k)
    fprintf('\b\b\b\b\b\b\b%3.2f%%', (j/length(k))*100); % Simulation Progress

    leader_distribution = Lambda*p_l(:, j)*copy_over_network;
    leader_velocity_distribution = Lambda*dp_l(:, j)*copy_over_network;
    p_j_array_database = q_next_k(1:dim_state-1, :);
    p_jdot_database = dq_next_k(1:dim_state-1, :);
    P_jddot_database = robot_accelerations;
    Delta_j = Delta_history(:, :, j);
    dDelta_j = dDelta_history(:, :, j);
    ddDelta_j = ddDelta_history(:, :, j);

    global_formation_errors = (p_j_array_database - Delta_j)*modified_laplacian...
                                - leader_distribution*leader_adjacency;
    velocity_sums = p_jdot_database*adjacency + leader_velocity_distribution*leader_adjacency;

    for i=1:mrs_size
        q_k(:, i) = q_next_k(:, i);
        measured_q = q_k(:, i);
        x_state_k = [measured_q; wheel_vel_i(:, i)];
        angle_val = x_state_k(3);
        p_i = p_j_array_database(:, i);

        q_i_library(:, j, i) = measured_q;
        vel_library(:, j, i) = v_i(:, i);
        tau_library(:, j, i) = tau_i(:, i);
        wheeled_vel_library(:, j, i) = wheel_vel_i(:, i);

        leader_error = leader_neighbour(i)*((p_i - p_l(1:dim_state-1, j)) - Delta_j(:, i));
        global_formation_error = global_formation_errors(:, i);
        % Evaluate local formation error excluding trajectory error
        formation_error(:, j, i) = abs(global_formation_error - leader_error);

        %%% Consensus Errors
        g_i_dis_val = g_dis(i,i);
        ith_basis_vector = std_basis(mrs_size, i);
        h_i_dis_val = -(p_j_array_database*adjacency + Delta_j*modified_laplacian...
                        + leader_distribution*leader_adjacency)*ith_basis_vector;
        dh_i_dis_val = -(p_jdot_database*adjacency + dDelta_j*modified_laplacian...
                        + leader_velocity_distribution*leader_adjacency)*ith_basis_vector;
        diffeomorphism_k_val = eta_func(p_j_array_database, p_jdot_database,...
                                     Delta_j, dDelta_j, p_l(:, j), dp_l(:, j),...
                                     x_state_k, g_i_dis_val, ith_basis_vector);
        relative_average_velocity_val = velocity_sums*ith_basis_vector/g_i_dis_val;

        % Pack numerical data into a vector in the SAME ORDER as 'all_params_sym'
        numerical_param_values = [diffeomorphism_k_val; h_i_dis_val; dh_i_dis_val;...
                                  angle_val; g_i_dis_val; relative_average_velocity_val];
        
        % Call the pre-built solver with the numerical data
        consensus_horizon_variables = consensus_nlpsolve('x0', consensus_opt_output(:, i), ...
                                                         'p', numerical_param_values, ...
                                                         'lbg', formation_lbg,...
                                                         'ubg', formation_ubg);

        consensus_opt_output(:, i) = full(consensus_horizon_variables.x);
        optimal_xi = consensus_opt_output(1:dim_state-1, i);
        optimal_nu = convert_xi_to_nu_func(P_jddot_database, ddDelta_j, ddp_l(:, j),...
                                           x_state_k, g_i_dis_val, ith_basis_vector,...
                                           optimal_xi);
        tau_i(:, i) = full(calculate_torque_func(x_state_k, optimal_nu));

        %%% update states
        state_solution_f = state_solver('x0', x_state_k, 'p', tau_i(:, i));
        x_state_next = full(state_solution_f.xf);
        dx_state_next = full(dx_state_func(x_state_k, tau_i(:, i)));
        robot_i_acceleration = full(global_robot_acceleration_func(x_state_k, tau_i(:, i)));

        q_next_k(:, i) = x_state_next(1:dim_state);
        dq_next_k(:, i) = dx_state_next(1:dim_state);
        robot_accelerations(:, i) = robot_i_acceleration;
        wheel_vel_i(:, i) = x_state_next(end-dim_vel+1:end);
        v_i(:, i) = [norm(dq_next_k(1:2, i)); dq_next_k(3, i)];
        consensus_opt_output(:, i) = [consensus_opt_output(2:end, i); consensus_opt_output(end, i)];
    end
end

save('Straight_Trajectory_Values.mat', 'formation_error', 'q_i_library', 'p_obs', 'obstacle_radius',...
        'k', 'p_l', 'vel_library', 'v_l');

%% IAE calculations
% IAE_1 = sum(formation_error(:, :, 1)*dt, 2);
% IAE_2 = sum(formation_error(:, :, 2)*dt, 2);
% IAE_3 = sum(formation_error(:, :, 3)*dt, 2);

%% Plot Results
%%% trajectories
figure(1);
hold on;
plot(q_i_library(1, :, 1), q_i_library(2, :, 1), 'r', 'LineWidth', 1.5);
plot(q_i_library(1, :, 2), q_i_library(2, :, 2), 'b', 'LineWidth', 1.5);
plot(q_i_library(1, :, 3), q_i_library(2, :, 3), 'g', 'LineWidth', 1.5);
plot(p_l(1, 1:length(k)), p_l(2, 1:length(k)), 'c--', 'LineWidth', 1.5);
rectangle('Position', [p_obs(1)-obstacle_radius,...
                       p_obs(2)-obstacle_radius,...
                       2*obstacle_radius, 2*obstacle_radius],...
                      'Curvature', [1, 1], 'FaceColor',[0 0 0]);
cord1x = reshape(q_i_library(1, 1, 1:3), [], 1);
cord1y = reshape(q_i_library(2, 1, 1:3), [], 1);
plot([cord1x; cord1x(1); p_l(1, 1)], [cord1y; cord1y(1); p_l(2, 1)], '.-',...
    'Color', [0.5, 0.5, 0.5, 0.5],'MarkerSize', 10, 'LineWidth', 2);
plot(nan, nan, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 10);
time = 170;
cord1x = reshape(q_i_library(1, time, 1:3), [], 1);
cord1y = reshape(q_i_library(2, time, 1:3), [], 1);
plot([cord1x; cord1x(1); p_l(1, time)], [cord1y; cord1y(1); p_l(2, int32(2*length(k)/5))], '.-',...
    'Color', [0.5, 0.5, 0.5, 0.5],'MarkerSize', 10, 'LineWidth', 2);
cord1x = reshape(q_i_library(1, int32(3.5*length(k)/5), 1:3), [], 1);
cord1y = reshape(q_i_library(2, int32(3.5*length(k)/5), 1:3), [], 1);
plot([cord1x; cord1x(1); p_l(1, int32(3.5*length(k)/5))], [cord1y; cord1y(1); p_l(2, int32(3.5*length(k)/5))], '.-',...
    'Color', [0.5, 0.5, 0.5, 0.5],'MarkerSize', 10, 'LineWidth', 2);
cord1x = reshape(q_i_library(1, length(k), 1:3), [], 1);
cord1y = reshape(q_i_library(2, length(k), 1:3), [], 1);
plot([cord1x; cord1x(1); p_l(1, length(k))], [cord1y; cord1y(1); p_l(2, length(k))], '.-',...
    'Color', [0.5, 0.5, 0.5, 0.5],'MarkerSize', 10, 'LineWidth', 2);
hold off;
legend('show', 'Location', 'northeast');
legend('Robot 1', 'Robot 2', 'Robot 3', 'Leader', 'formation tracker', 'Obstacle');
title("Robots' Trajectories", 'fontsize', 15);
xlabel('X coordinate [m]', 'fontSize', 12);
xlim([-3 12]);
ylabel('Y coordinate [m]', 'fontSize', 12);
ylim([-1.5 6]);
pbaspect([diff(xlim) diff(ylim) 1]);
box on
grid on
grid minor
set(gca, 'fontsize', 12);
set(gca, 'TickDir', 'out');

%{
%%% linear velocities
figure(2);
hold on;
plot(k, vel_library(1, :, 1), 'r', 'LineWidth', 1.5);
plot(k, vel_library(1, :, 2), 'b', 'LineWidth', 1.5);
plot(k, vel_library(1, :, 3), 'g', 'LineWidth', 1.5);
plot(k, v_l(1, 1:length(k)), 'c--', 'LineWidth', 1.5);
hold off;
title('Linear Velocities', 'fontsize', 15);
xlabel('time [s]', 'fontSize', 12);
ylabel('Linear velocities v_i [m/s]', 'fontSize', 12);
legend('show', 'Location', 'northeast');
legend('Robot 1', 'Robot 2', 'Robot 3', 'Leader');
box on
grid on
grid minor
set(gca, 'fontsize', 12);
set(gca, 'TickDir', 'out');

%%% angular velocities
figure(3);
hold on;
plot(k, vel_library(2, :, 1), 'r', 'LineWidth', 1.5);
plot(k, vel_library(2, :, 2), 'b', 'LineWidth', 1.5);
plot(k, vel_library(2, :, 3), 'g', 'LineWidth', 1.5);
plot(k, v_l(2, 1:length(k)), 'c--', 'LineWidth', 1.5);
hold off;
title('Angular Velocities', 'fontsize', 15);
xlabel('time [s]', 'fontSize', 12);
ylabel('Angular velocities $\omega_i$ [rad/s]', 'FontSize', 12, 'Interpreter', 'latex');
legend('show', 'Location', 'northeast');
legend('Robot 1', 'Robot 2', 'Robot 3', 'Leader');
box on
grid on
grid minor
set(gca, 'fontsize', 12);
set(gca, 'TickDir', 'out');

%%% X formation errors
figure(4);
hold on;
plot(k, formation_error(1, :, 1), 'r', 'LineWidth', 1.5);
plot(k, formation_error(1, :, 2), 'b', 'LineWidth', 1.5);
plot(k, formation_error(1, :, 3), 'g', 'LineWidth', 1.5);
hold off;
title('Absolute longitudinal axis formation error (x)', 'fontsize', 15);
xlabel('time [s]', 'fontSize', 12);
ylabel('x-position error', 'fontSize', 12);
legend('show', 'Location', 'northeast');
legend('Robot 1', 'Robot 2', 'Robot 3');
box on
grid on
grid minor
set(gca, 'fontsize', 12);
set(gca, 'TickDir', 'out');

%%% Y formation errors
figure(5);
hold on;
plot(k, formation_error(2, :, 1), 'r', 'LineWidth', 1.5);
plot(k, formation_error(2, :, 2), 'b', 'LineWidth', 1.5);
plot(k, formation_error(2, :, 3), 'g', 'LineWidth', 1.5);
hold off;
title('Absolute lateral axis formation error (y)', 'fontsize', 15);
xlabel('time [s]', 'fontSize', 12);
ylabel('y-position error', 'fontSize', 12);
legend('show', 'Location', 'northeast');
legend('Robot 1', 'Robot 2', 'Robot 3');
box on
grid on
grid minor
set(gca, 'fontsize', 12);
set(gca, 'TickDir', 'out');

%%% right wheel angular velocities
figure(6);
hold on;
plot(k, wheeled_vel_library(1, :, 1), 'r', 'LineWidth', 1.5);
plot(k, wheeled_vel_library(1, :, 2), 'b', 'LineWidth', 1.5);
plot(k, wheeled_vel_library(1, :, 3), 'g', 'LineWidth', 1.5);
hold off;
title('Right Wheel Angular Velocities', 'fontsize', 15);
xlabel('time [s]', 'fontSize', 12);
ylabel('$\dot{\Phi}_r$ [rad/s]', 'FontSize', 12, 'Interpreter', 'latex');
legend('show', 'Location', 'northeast');
legend('Robot 1', 'Robot 2', 'Robot 3');
box on
grid on
grid minor
set(gca, 'fontsize', 12);
set(gca, 'TickDir', 'out');

%%% left wheel angular velocities
figure(7);
hold on;
plot(k, wheeled_vel_library(2, :, 1), 'r', 'LineWidth', 1.5);
plot(k, wheeled_vel_library(2, :, 2), 'b', 'LineWidth', 1.5);
plot(k, wheeled_vel_library(2, :, 3), 'g', 'LineWidth', 1.5);
hold off;
title('Left Wheel Angular Velocities', 'fontsize', 15);
xlabel('time [s]', 'fontSize', 12);
ylabel('$\dot{\Phi}_l$ [rad/s]', 'FontSize', 12, 'Interpreter', 'latex');
legend('show', 'Location', 'northeast');
legend('Robot 1', 'Robot 2', 'Robot 3');
box on
grid on
grid minor
set(gca, 'fontsize', 12);
set(gca, 'TickDir', 'out');

%%% right wheel torque
figure(8);
hold on;
plot(k, tau_library(1, :, 1), 'r', 'LineWidth', 1.5);
plot(k, tau_library(1, :, 2), 'b', 'LineWidth', 1.5);
plot(k, tau_library(1, :, 3), 'g', 'LineWidth', 1.5);
hold off;
title('Right Wheel Torque', 'fontsize', 15);
xlabel('time [s]', 'fontSize', 12);
ylabel('$\tau_r\ [N\cdot m]$', 'FontSize', 12, 'Interpreter', 'latex');
legend('show', 'Location', 'northeast');
legend('Robot 1', 'Robot 2', 'Robot 3');
box on
grid on
grid minor
set(gca, 'fontsize', 12);

%%% left wheel torque
figure(9);
hold on;
plot(k, tau_library(2, :, 1), 'r', 'LineWidth', 1.5);
plot(k, tau_library(2, :, 2), 'b', 'LineWidth', 1.5);
plot(k, tau_library(2, :, 3), 'g', 'LineWidth', 1.5);
hold off;
title('Left Wheel Torque', 'fontsize', 15);
xlabel('time [s]', 'fontSize', 12);
ylabel('$\tau_l\ [N\cdot m]$', 'FontSize', 12, 'Interpreter', 'latex');
legend('show', 'Location', 'northeast');
legend('Robot 1', 'Robot 2', 'Robot 3');
box on
grid on
grid minor
set(gca, 'fontsize', 12);
set(gca, 'TickDir', 'out');
%}