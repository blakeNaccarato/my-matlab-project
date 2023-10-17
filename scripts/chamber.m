clc, clear, close all

%#ok<*UNRCH> allow unreachable branches of if/then/else statements

%% Common

L_ins_wall = 2.5; % (in)
L_ins_other = 0.5; % (in)
use_washer = false;
L_washer = 0.024; % (in)

k_aluminum = 200; % (W/m-K)
k_ins = 0.05; % (W/m-K)

A_one_leg = 0.4342; % (in^2)
N_legs = 12;
A_legs = N_legs*A_one_leg;

inch_to_meter = 0.0254; % (m/in)
Conductance = @(k, A, L) k*A/L * inch_to_meter; % (W/K)

pool_temp = 100; % (C)
air_temp = 25; % (C)

%% Lid (1, 2, 3)

% Lid (1 to 2)
A_lid = 230; % (in^2)
L_lid = 1.5; % (in)
K_lid = Conductance(k_aluminum, A_lid, L_lid);

% Lid insulation (2 to 3)
L_lid_ins = 0.5; % (in)
K_lid_ins = Conductance(k_ins, A_lid, L_lid_ins);

% Initial thermal network
N_nodes = 3;
thermal_network = ThermalNet(N_nodes, [K_lid, K_lid_ins]);
lid_bound = numel(thermal_network.Nodes);

%% Wall (1, 4, 5)

% Wall (1, 4)
A_wall = 470; % (in^2)
L_wall = 1/16; % (in)
K_wall = Conductance(k_aluminum, A_wall, L_wall);

% Wall insulation (4, 5)
K_wall_ins = Conductance(k_ins, A_wall, L_ins_wall);

% Branch existing network
pool_node = 1;
N_nodes = 2;
thermal_network.branch(pool_node, N_nodes, [K_wall, K_wall_ins])
wall_bound = numel(thermal_network.Nodes);

%% Floor (1, 6, 7)

% Floor (1, 6)
A_floor = 250; % (in^2)
L_floor = 1; % (in)
K_floor = Conductance(k_aluminum, A_floor, L_floor);

% Floor insulation (6, 7)
K_floor_ins = Conductance(k_ins, A_floor, L_ins_other);

% Branch existing network
N_nodes = 2;
thermal_network.branch(pool_node, N_nodes, [K_floor, K_floor_ins])
floor_bound = numel(thermal_network.Nodes);

%% Washer (6, 8)

floor_outside_node = 6;
N_nodes = 1;

if use_washer
    % Washer
    k_washer = 0.5; % (W/m-K) epoxy-fiberglass
    A_washer = pi/4 * (.5^2 -.25^2); % (in^2)
    K_washer = Conductance(k_washer, A_washer, L_washer);
    % Branch existing network
    K_total = N_legs * K_washer;
    thermal_network.branch(floor_outside_node, N_nodes, K_total)
else
    K_legs_instead = Conductance(k_aluminum, A_legs, L_washer);
    % Branch existing network
    thermal_network.branch(floor_outside_node, N_nodes, K_legs_instead)
end

%% Legs (8, 9)

L_legs = 16; % (in)
K_legs = Conductance(k_aluminum, A_legs, L_legs);

% Branch existing network
leg_start_node = 8;
N_nodes = 1;
thermal_network.branch(leg_start_node, N_nodes, K_legs)
leg_bound = numel(thermal_network.Nodes);

%% Leg Intermediate Nodes (8, 10, 11, 12, 9)

% Split the leg link into four, introducing three additional nodes
refine_between = [leg_start_node, leg_bound];
N_splits = 2;
thermal_network.refine(refine_between, N_splits)
% Get the most recently created nodes
new_end_node = numel(thermal_network.Nodes);
leg_refine_nodes = (leg_bound + 1):new_end_node;

% Apply one-fourth of the heat leak through insulation to nodes 8,10,11,12
A_legs_ins = N_legs * 130; % (in^2)
K_legs_ins = Conductance(k_ins, A_legs_ins, L_legs);
branch_nodes = [leg_start_node, leg_refine_nodes];
N_branches = numel(branch_nodes);
N_nodes = 1;
% Loop over branches
for i = 1:N_branches
    thermal_network.branch(branch_nodes(i), N_nodes, K_legs/N_branches);
end

new_end_node = numel(thermal_network.Nodes);
leg_ins_bounds = (new_end_node - N_branches + 1):new_end_node;

%% Boundaries

pool_node = 1;
thermal_network.addbound('temp', pool_node, pool_temp)

air_nodes = [lid_bound, wall_bound, floor_bound, leg_bound, leg_ins_bounds];
N_bounds = numel(air_nodes);
for i = 1:N_bounds
    thermal_network.addbound('temp', air_nodes(i), air_temp)
end

%% Plot
thermal_network.solve()
thermal_network.getfluxes()
figure, thermal_network.plot('node');
% figure, thermal_network.plot('cond');
figure, thermal_network.plot('flux');
