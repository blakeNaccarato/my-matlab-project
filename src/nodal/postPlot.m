%% Plot Parameters

alpha = 0.5;
posit = [2000 100 500 750];

DX = 20;
DY = -80;

%% Node Spacing in Meter Space

% preallocate X and Y
X = zeros(1,endNode);
Y = X;
% sidechain spacing
dymin = 0.01;
drmin = sqrt(dymin^2/2);

% find node positions along copper post
Y(1:11) = cumsum([0 dy(1:end-1)]);
LCu = Y(11);
% find node position of stainless steel sample
Y(12) = Y(11) + dy(end);
% find node position of pool
Y(13) = Y(12) + dymin;

% find node position of mount and air side
X(14) = X(1) + drmin;
Y(14) = Y(1) - drmin;
X(15) = X(14) + drmin;
Y(15) = Y(14) - drmin;
% find node position of bolt mount and air side
X(16) = -drmin;
Y(16) = Y(5) + drmin;
X(17) = X(16) - drmin;
Y(17) = Y(16) + drmin;

% find node position of insulation
coreIDs = 1:9;
insIDs = 18:2:34;
Y(insIDs) = Y(coreIDs);
X(insIDs) = dIns/2;
% find node position of air outside insulation
airIDs = 19:2:35;
Y(airIDs) = Y(coreIDs);
X(airIDs) = X(insIDs) + dymin;

% find node position of natural convection
coreIDs = 10:11;
convIDs = 36:37;
Y(convIDs) = Y(coreIDs);
X(convIDs) = dTop;

%% Converting Node Spacing in Meters to Pixels

% known node positions in meters
X0 = 0;
Y0 = 0;
X1 = dBot/2;
Y1 = LCu;

% known node positions in pixels
filename = 'post2.png'; % filename = 'post.png';
Xi0 = 2608; % Xi0 = 363;
Yi0 = 5592; % Yi0 = 629;
Xi1 = 3206; % Xi1 = 451;
Yi1 = 2056; % Yi1 = 110;

% convert node positions from meters to pixels
Xi = Xi0 + X*(Xi1-Xi0)/(X1-X0);
Yi = Yi0 + Y*(Yi1-Yi0)/(Y1-Y0);

%% Plot Nodal Network

% plot background image
f = figure(4); clf
f.Position = posit;
h = imshow(filename); hold on
h.AlphaData = alpha;

% plot directed network graph
update(net)
G = digraph(net.Fluxes);
G = flipedge(G, find(G.Edges.Weight<0));
h = plot(G, 'ArrowSize' , 10 ,...
            'EdgeAlpha' ,  1 ,...
            'EdgeColor' , 'r',...
            'LineWidth' ,  1 ,...
            'Marker'    , 'o',...
            'MarkerSize',  3 ,...
            'NodeColor' , 'r',...
            'NodeLabel' ,  []);
set(gca,'Position',[0 0 1 1])

% move nodes to pixel positions
h.YData = Yi;
h.XData = Xi;

% plot node labels
text(Xi+DX, Yi+DY, num2str((1:37)'),'FontSize',9,'Color','r')

%% Plot Fluxes

% plot background image
f = figure(5); clf
f.Position = posit;
h = imshow(filename); hold on
h.AlphaData = alpha;

% plot directed network graph
update(net)
fluxes = round(net.Fluxes,0);
G = digraph(fluxes);
G = flipedge(G, find(G.Edges.Weight<0));
W = abs(G.Edges.Weight);
h2 = plot(G, 'EdgeLabel' ,  W ,...
            'ArrowSize' , 10 ,...
            'EdgeAlpha' ,  1 ,...
            'EdgeColor' , 'r',...
            'LineWidth' ,  1 ,...
            'Marker'    , 'o',...
            'MarkerSize',  3 ,...
            'NodeColor' , 'r',...
            'NodeLabel' ,  []);
set(gca,'Position',[0 0 1 1])

% move nodes to pixel positions
h2.YData = Yi;
h2.XData = Xi;

% hide fluxes in inconvenient spots
hideFluxIDs = [2 3 5:7 9:14 16 18 20 21 23];
N = numel(hideFluxIDs);
for i = 1:N
    h2.EdgeLabel{hideFluxIDs(i)} = '';
end

% hide nodes in inconvenient spots
hideNodeIDs = [18 20 22 24 26 28 30 32 34];
Xi2 = Xi;
Xi2(hideNodeIDs) = [];
Yi2 = Yi;
Yi2(hideNodeIDs) = [];
IDs = (1:37)';
IDs(hideNodeIDs) = [];

% plot node labels
text(Xi2+DX, Yi2+DY, num2str(IDs),'FontSize',9,'Color','r')

%% Plot Temperatures

% plot background image
f = figure(6); clf
f.Position = posit;
h = imshow(filename); hold on
h.AlphaData = alpha;

% plot directed network graph
update(net)
G = digraph(net.Fluxes);
G = flipedge(G, find(G.Edges.Weight<0));
h3 = plot(G,'ArrowSize' , 10 ,...
            'EdgeAlpha' ,  1 ,...
            'EdgeColor' , 'r',...
            'LineWidth' ,  1 ,...
            'Marker'    , 'o',...
            'MarkerSize',  3 ,...
            'NodeColor' , 'r',...
            'NodeLabel' ,  []);
set(gca,'Position',[0 0 1 1])

% move nodes to pixel positions
h3.YData = Yi;
h3.XData = Xi;

% hide nodes in inconvenient spots
hideNodeIDs = [15:2:35 16];
Xi3 = Xi;
Xi3(hideNodeIDs) = [];
Yi3 = Yi;
Yi3(hideNodeIDs) = [];

% plot node labels
temps = num2str(round([net.Nodes.Temperature]'));
temps = strcat(num2str((1:37)'),':',temps,'°C');
temps(hideNodeIDs,:) = [];
text(Xi3+DX, Yi3+DY, temps,'FontSize',9,'Color','r') % text(Xi3+3,Yi3-8,temps,'FontSize',9,'Color','r')
