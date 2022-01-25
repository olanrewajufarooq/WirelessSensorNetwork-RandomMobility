close all;
clc;
clear;

%% Simulation Details
disp("Welcome to Wireless Sensor Simulation\n");
disp("............................................................\n");
disp("............................................................\n\n");
 
% Initialization Inputs
max_dimension = 100; % Maximum Dimension of the WSN Plot

initial_energy = 500e-3; % Initial node energy
transceiver_energy = 50e-9; % Energy Required for transmission and receiving of packet
ener_amp = 100e-12; % Amplification Energy
ener_agg = 100e-12; % Aggregation Energy

% Simulation Parameters
n = 100; % Number of nodes

sn = 3; % Number of mobile sink
sn_method = 'random'; % the mobile sink can be either selected randomly 'random' or evenly spaceed 'even'.

rounds = 100; % Number of rounds per simulation
k = 80000; % Bits transmitted per packet

% Clustering Paramters
n_clusters = 8;

% Mobility Parameters
min_dist = 1; % Minimum mobility for sensor nodes (in m)
max_dist = 3; % Maximum mobility for sensor nodes (in m)
sn_min_dist = 2; % Minimum mobility for sink nodes (in m)
sn_max_dist = 5; %Maximum mobility for sink nodes (in m)
min_visit_dist = 1.5; % Minimum distance to affirm visitation by sink nodes (in m)
mob_params = containers.Map({'min_dist', 'max_dist', 'sn_min_dist', 'sn_max_dist', 'min_visit_dist'}, {min_dist, max_dist, sn_min_dist, sn_max_dist, min_visit_dist});

%% Parameters Initialization
[dims, ener] = param_init(max_dimension, initial_energy, transceiver_energy, ener_agg, ener_amp);

%% Initialization of the WSN
[SN, ms_ids] = createWSN(n, sn, sn_method, dims, ener('init'), rounds);

%% Smiluation of the WSN
[SN, round_params, sim_params] = simulation_rounds(rounds, SN, dims, ener, k, ms_ids, n_clusters, mob_params);

%% Data Visualisation Conclusion
plot_simulation(SN, rounds, dims)

%% Lifetime and Stability Periods.

fprintf('\n\nSimulation Summary\n\n')
fprintf('Stability Period: %d secs\n', round(round_params('stability period'), 2))
fprintf('Stability Period Round: %d\n', round_params('stability period round'))
fprintf('Lifetime: %d secs\n', round(round_params('lifetime'), 2))
fprintf('Lifetime Round: %d\n', round_params('lifetime round'))