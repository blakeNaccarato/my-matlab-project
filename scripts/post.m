clear

%% Inputs

% conversion factor
inToM = 0.0254; % (m/in)

% geometry inputs
dBot    = 2*inToM; % (m) OD of copper post before taper
dTop    = 1.25*inToM; % (m) OD of copper post after taper
dTap    = mean([dBot dTop]); % (m) average OD of copper post during taper
dCart   = 0.5*inToM; % (m) cartridge hole diameter
tIns = 0.25*inToM; % (m) insulation thickness

% thermal conductivities
kCu  = 400; % (W/m-K)
kSS  = 20; % (W/m-K)
kIns = 1; % (W/m-K)

% heat transfer coefficients
hAir  = 10; % (W/m^2-K)
hWat  = 1e3; % (W/m^2-K)
hNucl = 1e6; % (W/m^2-K)

% boundary conditions
heatIn = 1.5e3; % (W)
poolT  = 100; % (C)
airT   = 25; % (C)

%% Geometry

% link lengths    betw. nodes
dy = [2.00...   % (in)  1 -  2, base, w/ cart holes
      0.25...   % (in)  2 -  3, base, w/o cart holes
      0.40...   % (in)  3 -  4, taper
      0.50...   % (in)  4 -  5, before mount
      0.25...   % (in)  5 -  6, mount
      0.45...   % (in)  6 -  7, after mount to 1st TC hole
      0.65...   % (in)  7 -  8, between 1st and 2nd TC holes
      0.65...   % (in)  8 -  9, between 2nd and 3rd TC holes
      0.60...   % (in)  9 - 10, 3rd TC hole to o-ring
      0.15...   % (in) 10 - 11, o-ring to top of copper post
      0.25]...  % (in) 11 - 12, stainless steel sample
      *inToM;   % convert inch to meter

numL = numel(dy);   % number of links to be defined
d = zeros(1,numL);  % (m) preallocate vector of ODs between nodes
% link diameters    betw. nodes
d(1:5) = [dBot...   % (m) 1 - 2, before taper, cartridge holes to be subtracted
          dBot...   % (m) 2 - 3, before taper, no cartridge holes
          dTap...   % (m) 3 - 4, average taper diameter
          dTop...   % (m) 4 - 5, after taper, before mounting ring
          dBot];    % (m) 5 - 6, mounting ring
d(6:end) = dTop;    % (m) 6 - 11, after mount to top of mounted sample

% link cross-sectional areas
A = pi/4*d.^2;                      % (m^2) area of links
A(1) = A(1) - pi/4*(4*dCart.^2);    % (m^2) subtract cartridge holes from link 1
% link materials
k(1:numL-1) = kCu;  % (W/m-K) copper post
k(numL)     = kSS;  % (W/m-K) stainless sample
% link conductances
K = k.*A./dy;               % (W/K) conductive links in copper post and sample
K(end+1) = hNucl*A(end);    % (W/K) nucleate boiling link to pool
numL = numL + 1;            % update link count to include convective link

%% Thermal Network before Incorporating Losses
% thermal conduction from the base of the post to the boiling pool

numN = numL+1;                  % number of nodes is (number of links + 1)
net = ThermalNet(numL+1, K);    % create network with specified conductances
% add BCs
addbound(net,'heat',1,heatIn)       % apply heat input BC to node 1
addbound(net,'temp',numN,poolT)     % apply fixed temp BC to node 13

%% Mount Losses
% get losses through stainless-steel mounts

mIDs = [1 5];           % nodes where mount plates exist
numM = numel(mIDs);     % number of mounts
% mount link lengths
dyMount = dy(5);        % (m) thickness of mount plates
% mount link cross-sectional areas
mountA = A(mIDs);                       % (m^2) cross-sectional area of mounts
mountA(2) = mountA(2) - pi/4*dTop^2;    % (m^2) exclude copper post core from 2nd mount
% mount thermal conductances
mountK(1) = kSS* mountA(1)/dyMount;     % (W/K) conduction through mount
mountK(2) = hAir*mountA(2);             % (W/K) convective loss to air

for i = 1:numM                      % for each mounting node
    branch(net,mIDs(i),2,mountK)        % apply conductances to a two-node sidechain
    endNode = numel(net.Nodes);         % get the latest node in the net
    addbound(net,'temp',endNode,airT)   % apply ambient temperature BC to it
end

%% Radial Losses to Air
% get radial conductances thru insulation to air

% geometry
airIDs = 1:9;           % nodes with radial losses to air
numA = numel(airIDs);   % number of nodes with radial losses to air
dPost = d(airIDs);      % (m) diameter of post at nodes with losses
dyPost = dy(airIDs);    % (m) length of links above nodes with losses
dIns = dPost + 2*tIns;  % (m) OD at outside of insulation at these nodes
% conductances
airK = zeros(2,numA);                           % (W/K) preallocate radial thermal conductance array
airK(1,:) = kIns*2*pi*dy(1:9)./log(dIns./dPost);   % (W/K) radial thermal conductance through insulation
airK(2,:) = hAir  *pi*dIns.*dyPost;             % (W/K) convective conductance over circumferential area

% create radial branches
for i = 1:numA                          % for each node with radial losses to air
    branch(net, airIDs(i), 2, airK(:,i))    % apply conductances to a two-node sidechain
    endNode = numel(net.Nodes);             % get the latest node in the net
    addbound(net,'temp',endNode,airT)       % apply ambient temperature BC to it
end

%% Radial Losses to Water
% radial losses by natural convection to boiling pool

% geometry
watIDs = 10:11;                 % nodes with radial losses to water
numW = numel(watIDs);           % number of nodes with radial losses to water
dPost = d(watIDs);              % (m) diameter of post at nodes with losses
dyPost = dy(watIDs);            % (m) length of links above nodes with losses
watK = hWat*pi*dPost.*dyPost;   % (W/K) natural convection conductance to pool

% get radial conductances to water
for i = 1:numW                      % for each node with natural convection
    branch(net, watIDs(i), 1, watK(i))  % apply conductance to one-node sidechain
    endNode = numel(net.Nodes);         % get the latest node in the net
    addbound(net,'temp',endNode,poolT)  % apply pool temperature BC to it
end

%% Solve Before Refinement
update(net)
solve(net)                          % solve for unknown temperatures
getfluxes(net)                      % compute fluxes

% figure(3), hF = plot(net,'flux');   % plot nodal network with temperatures at nodes and flux along paths
postPlot;

%% Refine Post

% add more nodes in the cartridge heater section to distribute heat input
numRef = 4;
refine(net,[1 2],numRef)

update(net)
solve(net)                          % solve for unknown temperatures
getfluxes(net)                      % compute fluxes

% figure(3), hF = plot(net,'flux');   % plot nodal network with temperatures at nodes and flux along paths
reduce(net,[1 2])

update(net)
solve(net)                          % solve for unknown temperatures
getfluxes(net)                      % compute fluxes

postPlot;

%% Refine Water 1

nodePairID = [10 11];
refine(net,nodePairID,numRef)

numL = 2^numRef;
nodeEnd = numel(net.Nodes);
nodeStart = nodeEnd - (2^numRef-2);

link = setdiff(net.Nodes([nodePairID(1) nodeStart]).Links);
link = setdiff(link, net.Nodes(nodePairID(1)-1).Links);
newCond = link.Conductance/numL;
link.Conductance = newCond;

for nodeID = nodeStart:nodeEnd
    newNodeID = nodeID + 2^numRef-1;
    branch(net,nodeID,1,newCond)
    addbound(net,'temp',newNodeID,poolT)  % apply pool temperature BC to it
end

update(net)
solve(net)                          % solve for unknown temperatures
getfluxes(net)                      % compute fluxes

% figure(3), hF = plot(net,'flux');   % plot nodal network with temperatures at nodes and flux along paths
reduce2(net,nodePairID)

update(net)
solve(net)                          % solve for unknown temperatures
getfluxes(net)                      % compute fluxes

postPlot;

%% Refine Water 2

nodePairID2 = [11 12];
refine(net,nodePairID2,numRef)

numL = 2^numRef;
nodeEnd = numel(net.Nodes);
nodeStart = nodeEnd - (2^numRef-2);

link = setdiff(net.Nodes([nodePairID2(1) nodeStart]).Links);
link = setdiff(link, net.Nodes(nodePairID2(1)-1).Links);
newCond = link.Conductance/numL;
link.Conductance = newCond;

for nodeID = nodeStart:nodeEnd
    newNodeID = nodeID + 2^numRef-1;
    branch(net,nodeID,1,newCond)
    addbound(net,'temp',newNodeID,poolT)  % apply pool temperature BC to it
end

update(net)
solve(net)                          % solve for unknown temperatures
getfluxes(net)                      % compute fluxes

% figure(3), hF = plot(net,'flux');   % plot nodal network with temperatures at nodes and flux along paths
reduce2(net,nodePairID2)

update(net)
solve(net)                          % solve for unknown temperatures
getfluxes(net)                      % compute fluxes

postPlot;

%% Plot

% figure(1), hN = plot(net,'node');   % plot nodal network
% figure(2), hC = plot(net,'cond');   % plot nodal network with conductance along paths
% figure(3), hF = plot(net,'flux');   % plot nodal network with temperatures at nodes and flux along paths
