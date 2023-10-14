classdef ThermalNet < handle
%THERMALNET A thermal network of nodes linked by conductive paths.
%   THERMALNET generates a network of 10 nodes and the 9 links connecting
%   them.
%
%   THERMALNET(N) generates a network of N nodes and N-1 links connecting
%   them.
%
%   THERMALNET(N,K) generates a network of N nodes and N-1 links connecting
%   them, with a vector K of length N-1 specifying the thermal conductances
%   of the links.

%% Properties

    properties
        Nodes Node              % list of nodes that are in the network
        Links Link              % list of paths that are in the network

        Adjacency logical       % adjacency array used in network plotting
        Conductances double     % conductance array used in solving the network
        Fluxes double           % flux array to be computed after solving for temperatures

        Knowns uint16           % keeps track of boundary conditions set on the network
    end

%% Methods

    methods
        %% Constructor

        function net = ThermalNet(varargin)

            % handle input arguments
            if nargin == 1                  % if arguments are given
                numN = varargin{1};             % the first arg is the number of nodes
            elseif nargin == 2              % if more than one arg
                numN = varargin{1};             % the first arg is the number of nodes
                K = varargin{2};                % the 2nd arg is link conductances
            else                            % else
                numN = 10;                      % use 10 nodes by default
            end

            % build the network
            nodes(numN) = Node;             % generate numN Node objects
            links(numN-1) = Link;           % generate numN-1 Link objects
            net.joinchain(nodes, links)     % connect Nodes with Links
            net.Nodes = nodes;              % store Nodes in the ThermalNet object
            net.Links = links;              % store Links in the ThermalNet object

            % apply conductances, if given
            if nargin == 2
                net.setcond(links,K)
            end
        end

        %% Networking

        function branch(net, rootNodeID, numN, varargin)
        %BRANCH Generate a sidechain on the network.
        %   BRANCH(net, rootNodeID, N) generates a sidechain of N nodes in
        %   net. The sidechain branches off of the node specified by
        %   rootNodeID.
        %
        %   BRANCH(net, rootNodeID, N, K) generates an N-node sidechain off
        %   of rootNodeID, with thermal conductances in the sidechain links
        %   specified by K, an N-element vector.
        %
        %   See also ThermalNet/join, ThermalNet/netjoin.

            rootNode = net.Nodes(rootNodeID);   % get the root node

            newNodes(numN) = Node;              % generate enough nodes for the sidechain
            links(numN) = Link;                 % generate enough links for the sidechain

            nodes = [rootNode newNodes];        % operating on the rootNode and newNodes...
            net.joinchain(nodes, links)         % connect these nodes with links

            net.Nodes = [net.Nodes nodes(2:end)];   % add the new nodes to the network
            net.Links = [net.Links links];          % add the new links to the network

            % apply conductances, if given
            if nargin > 3
                K = varargin{1};
                net.setcond(links,K)
            end
        end

        function join(net, rootNodeID, joinNodeIDs, varargin)
        %JOIN Join any number of nodes in the network to a root node.
        %   JOIN(net, rootNodeID, joinNodeIDs) connects all of the nodes in
        %   net specified by joinNodeIDs to the node specified by
        %   rootNodeID.
        %
        %   See also ThermalNet/branch, ThermalNet/netjoin.

            rootNode = net.Nodes(rootNodeID);       % get the root node
            joinNodes = net.Nodes(joinNodeIDs);     % get the nodes that will be joined to rootNode
            numL = numel(joinNodes);                % the number of links needed for the operation
            links(numL) = Link;                     % generate enough links for the operation

            rootNode.Neighbors = [rootNode.Neighbors joinNodes];    % update the rootNode's list of neighbors
            rootNode.Links = [rootNode.Links links];                % update the rootNode's list of links

            for i = 1:numL                                                  % for each node to be joined
                joinNodes(i).Neighbors = [joinNodes(i).Neighbors rootNode];     % update the node's neighbors
                joinNodes(i).Links = [joinNodes(i).Links links(i)];             % update the node's paths

                % store the nodeIDs in the connecting link in the correct order
                if rootNodeID < joinNodeIDs(i)
                    links(i).Nodes = [rootNode joinNodes(i)];
                else
                    links(i).Nodes = [joinNodes(i) rootNode];
                end

            end

            % apply conductances, if given
            if nargin > 3
                K = varargin{1};
                net.setcond(links,K)
            end

            net.Links = [net.Links links];                          % add the new links to the network
        end

        function netjoin(rootNet, rootNodeID, joinNet, joinNodeID)
        %NETJOIN Join a second network to the existing network.
        %   NETJOIN(rootNet, rootNodeID, joinNet, joinNodeID) adds the
        %   joinNet network to the rootNet network. The node specified by
        %   joinNodeID in joinNet is connected to the node specified by
        %   rootNodeID in rootNet.
        %
        %   See also ThermalNet/branch, ThermalNet/join.

            % add the links from the network to be joined to the root network
            rootNet.Links = [rootNet.Links joinNet.Links];

            % add the nodes from the network to be joined to the root network
            rootNodes = rootNet.Nodes;              % get all nodes in the root network
            joinNodes = joinNet.Nodes;              % get all nodes in the network to be joined
            rootNet.Nodes = [rootNodes joinNodes];

            % get the node addresses for nodes coming from the network to be joined
            address(joinNet)
            oldAddress = [joinNodes.Address];

            % update these to their new addresses in the root network
            address(rootNet)
            newAddress = [joinNodes.Address];

            % update the nodeID of the node to be joined to reflect its new address
            joinNodeID = newAddress(joinNodeID==oldAddress);

            % get the knowns/BCs from the network to be joined
            joinKnowns = joinNet.Knowns;
            numK = size(joinKnowns,1);

            % update the addresses of these knowns/BCs to reflect their position in the root network
            if numK     % only if there are any knowns/BCs in the network to be joined
                newAddress = repmat(newAddress,numK,1);
                joinKnowns(:,1) = newAddress(joinKnowns(:,1)==oldAddress);
                rootNet.Knowns = [rootNet.Knowns; joinKnowns];  % add the knowns to the root network
            end

            % connect the root network and the network to be joined at the specified point
            join(rootNet, rootNodeID, joinNodeID)
        end

        function refine(net, nodePairID, multiple)
            boundNodes = net.Nodes(nodePairID);     % nodes bounding the link to be refined
            numL = 2^multiple;                      % number of links after refinement
            numN = numL-1;                          % number of nodes to be created

            [oldLink, ind1, ind2] = intersect(boundNodes.Links);  % the link shared by the two boundary nodes
            oldCond = oldLink.Conductance;          % the conductance of the original link
            newCond = oldCond*numL;                 % get equivalent conductance for shorter links

            branch(net, nodePairID(1), numN, newCond)   % create refined sidechain off of the first node in the pair
            lastN = numel(net.Nodes);                   % find the ID of the last node in the sidechain
            join(net, lastN, nodePairID(2), newCond)    % join that node to the second node in the pair to be refined

            delete(oldLink)
            boundNodes(1).Links(ind1) = [];
            boundNodes(2).Links(ind2) = [];
            net.Links(~isvalid(net.Links)) = [];

            idx = boundNodes(1).Neighbors == boundNodes(2);
            boundNodes(1).Neighbors(idx) = [];
            idx = boundNodes(2).Neighbors == boundNodes(1);
            boundNodes(2).Neighbors(idx) = [];

            if ~isempty(boundNodes(1).BoundHeat)

                boundNodeHeat = boundNodes(1).BoundHeat/numL;
                rembound(net,nodePairID(1))
                addbound(net, 'heat', nodePairID(1), boundNodeHeat)

                totN = numel(net.Nodes);
                refNodeIDs = (1+totN-numN):totN;
                for i = 1:numN
                    addbound(net, 'heat', refNodeIDs(i), boundNodeHeat)
                end
            end
        end

        function reduce(net, nodePairID)
            figure,
            [~, G] = plot(net,'node');
            close

            nodeIDs = shortestpath(G,nodePairID(1),nodePairID(2));
            nodes = net.Nodes(nodeIDs);

            numN = numel(nodes);
            numL = numN-1;

            idx = net.Knowns(:,1) > nodeIDs(end-1);
            net.Knowns(idx,1) = net.Knowns(idx,1) - (numN-2);

            cond = zeros(1,numL);
            flux = zeros(1,numL);
            for i = 1:numL
                [link ind1 ind2] = intersect(nodes([i i+1]).Links);
                cond(i) = link.Conductance;
                flux(i) = link.Flux;
                delete(link)
                nodes(i).Links(ind1) = [];
                nodes(i+1).Links(ind2) = [];
            end
            net.Links(~isvalid(net.Links)) = [];
            equivCond = sum(1./cond)^-1;

            nodeTemp1 = nodes(1).Temperature;
            rembound(net,nodeIDs(1))
            addbound(net,'temp',nodeIDs(1),nodeTemp1)

            nodeTemp2 = nodes(end).Temperature;
            rembound(net,nodeIDs(end))
            addbound(net,'temp',nodeIDs(end),nodeTemp2)

            for i = 2:numN-1
                rembound(net,nodeIDs(i))
                delete(nodes(i))
            end
            net.Nodes(~isvalid(net.Nodes)) = [];

            join(net,nodePairID(1),nodePairID(2),equivCond)

            link = intersect(nodes([1 end]).Links);
            if sign(flux(1)) ~= sign(flux(end))
                equivFlux = -flux(end);
            else
                equivFlux = flux(end);
            end

            link.Flux = equivFlux;
            link.Conductance = equivFlux/(nodeTemp1-nodeTemp2);

            nodes( 1 ).Neighbors(~isvalid(nodes( 1 ).Neighbors)) = [];
            nodes(end).Neighbors(~isvalid(nodes(end).Neighbors)) = [];

        end

        %% Addressing

        function address(net)
        %ADDRESS Give each node in the network a unique address.
        %   ADDRESS(net) gives each node in net a unique address.
        %
        %   See also ThermalNet/update.

            nodes = net.Nodes;          % get all nodes in the network
            numN = numel(nodes);        % count them
            for i = 1:numN              % for each node
                nodes(i).Address = i;       % assign a sequential address
            end
        end

        function update(net)
        %UPDATE Update adjacency and conductance arrays for the network.
        %   UPDATE(net) Updates adjacency and conductance arrays in net.
        %
        %   See also ThermalNet/address.

            address(net)                    % update node addresses

            numN = numel(net.Nodes);        % get the number of nodes
            links = net.Links;              % get the links in the network
            numL = numel(links);            % get the number of links

            adjacency = false(numN);        % preallocate adjacency
            conductances = zeros(numN);     % preallocate conductances
            for i = 1:numL                  % for each path

                % set adjacency to true for nodes connected to that path
                nodePairID = [links(i).Nodes.Address];
                adjacency(nodePairID(1), nodePairID(2)) = true;

                % set conductances for nodes connected to that path
                conductance = links(i).Conductance;
                conductances(nodePairID(1), nodePairID(2)) = conductance;
                conductances(nodePairID(2), nodePairID(1)) = conductance;
            end

            net.Adjacency = sparse(adjacency);          % update adjacency
            net.Conductances = sparse(conductances);    % update conductances
        end

        %% BCs and Knowns

        function addbound(net, type, nodeID, boundvalue)
        %ADDBOUND Add a boundary condition to the network.
        %   ADDBOUND(net, 'type', nodeID, boundvalue) adds a boundary
        %   condition of 'type' ('temp' or 'heat') and with boundvalue (in
        %   units of C or W) to the node specified by nodeID in net.

            node = net.Nodes(nodeID);   % get the node to be operated on

            if strcmp(type, 'temp')             % if it's a temp BC
                node.Temperature = boundvalue;      % set the node temp
                net.Knowns = [net.Knowns;...
                             nodeID 1];             % add it to the network knowns
            elseif strcmp(type, 'heat')         % if it's a heat BC
                node.BoundHeat = boundvalue;        % set the node heat
                net.Knowns = [net.Knowns;...
                             nodeID 2];             % add it to the network knowns
            end
        end

        function rembound(net, nodeID)
            nodes = net.Nodes;
            nodes(nodeID).Temperature = [];
            nodes(nodeID).BoundHeat = [];
            remIdx = net.Knowns(:,1)==nodeID;
            net.Knowns(remIdx,:) = [];
        end
        %% Solver

        function solve(net)
        %SOLVE Solve for nodal temperatures in the network.
        %   SOLVE(net) solves for nodal temperatures in net.
        %
        %   See also ThermalNet/getfluxes.

%             update(net)                 % update adjacency and conductances
            conds = net.Conductances;   % get the conductances array

            nodes = net.Nodes;          % get all the nodes in the network
            numN = numel(nodes);        % count them

            % get coeffs and rhs in coeffs*T=rhs, before applying BCs
            coeffs = -conds + diag(sum(conds));     % (W/K)
            rhs = zeros(numN,1);                    % (W)

            % apply BCs: remove/adjust rows/cols in coeffs and rhs
            knowns = net.Knowns;        % get the knowns
            numK = size(knowns,1);      % count them
            % remove rows/cols from higher addresses first to avoid readressing
            [~,idx] = sort(knowns(:,1),'descend');
            knowns = knowns(idx,:);
            for i = 1:numK              % for each known

                known = knowns(i,:);        % get the BC nodes and BC types
                nodeID = known(1);          % get the BC nodeID
                node = nodes(nodeID);       % get the BC node

                if known(2) == 1            % if temperature is known
                    temp = node.Temperature;    % get the node temperature
                    nodes(nodeID) = [];         % don't solve this node's temp
                    numN = numN-1;              % reduce unknowns by one

                    % get the conductances from the BC to neighboring nodes
                    K = -coeffs([1:nodeID-1 nodeID+1:end], nodeID);

                    coeffs(nodeID,:) = [];      % remove the row of the BC
                    coeffs(:,nodeID) = [];      % remove the column of the BC
                    rhs(nodeID) = [];           % remove the rhs of the BC

                    % move K*temp to rhs in nodes adjacent to the known temp
                    rhs = rhs + K*temp;

                elseif known(2) == 2                    % if heat is known
                    boundheat = node.BoundHeat;             % get the node heat
                    rhs(nodeID) = rhs(nodeID) + boundheat;  % add the heat to rhs as known
                end
            end

            % solve for nodal temperatures and assign sol'n to each node
            temperatures = coeffs\rhs;
            for i = 1:numN
                nodes(i).Temperature = temperatures(i);
            end
        end

        function getfluxes(net)
        %GETFLUXES Obtain heat fluxes once temperatures have been solved.
        %   GETFLUXES(net) obtains heat fluxes in net once temperatures
        %   have been solved.
        %
        %   See also ThermalNet/solve.

            numN = numel(net.Nodes);    % get the number of nodes
            links = net.Links;          % get all the links in the network
            numL = numel(links);        % count them

            fluxes = zeros(numN);       % preallocate flux array
            for i = 1:numL              % for each link

                % get the nodes attached to that link and their addresses
                nodes = links(i).Nodes;
                nodePairID = [nodes.Address];

                % get the flux in the path, positive from low-to-high address
                deltaT = diff([nodes.Temperature]);
                conductance = links(i).Conductance;
                flux = -conductance*deltaT;

                % assign the flux to the link and the flux array
                links(i).Flux = flux;
                fluxes(nodePairID(1), nodePairID(2)) = flux;
            end

            net.Fluxes = sparse(fluxes);    % store the flux array
        end

        %% Plotting

        function labels = gettemplabels(net)
        %GETTEMPLABELS Get node labels in the form "NodeID=Temp".
        %   GETTEMPLABELS(net) gets node labels for net in the form
        %   "NodeID=Temp".

            nodes = net.Nodes;
            numN = numel(nodes);

            nodeTemps = [nodes.Temperature];

            labels = cell(1, numN);
            for i = 1:numN
               nodeStr = num2str(i);
               valStr = sprintf('%.2f',nodeTemps(i));
               label = [nodeStr '=' valStr];
               labels{i} = label;
            end
        end

        function [h G] = plot(net, type)
        %PLOT Plot the nodal network.
        %   PLOT(net,type) plots the nodal network for net of 'type'
        %   ('node', 'flux', or 'cond'). For 'node', just plot the node
        %   numbers and connectivity. For 'cond', plot the node numbers and
        %   conductances along paths. FOr 'flux', plot the node numbers,
        %   temperatures, and fluxes along paths.

            update(net)
            lw = 3;

            if strcmp(type, 'node')
                adjacency = net.Adjacency;
                adjacency = adjacency + adjacency';
                G = graph(adjacency);
                LW = 0.5*lw;
                h = plot(G);
            elseif strcmp(type,'cond')
                conductances = round(net.Conductances,3);
                G = graph(conductances);
                W = G.Edges.Weight;
                LW = lw*W/max(abs(W));
                h = plot(G, 'EdgeLabel',W);
            elseif strcmp(type,'flux')
                fluxes = round(net.Fluxes,0);
                labels = gettemplabels(net);
                G = digraph(fluxes,labels);
                G = flipedge(G, find(G.Edges.Weight<0));
                W = abs(G.Edges.Weight);
                LW = lw*abs(W)/max(abs(W));
                h = plot(G, 'EdgeLabel',W, 'ArrowSize',15);
            end

            axis off equal
            set(gcf, 'Color','w')
            set(gca, 'Position',[0 0 1 1])
            layout(h,'force')
            h.EdgeColor = 'k';
            h.EdgeAlpha = 1;
            h.LineWidth = LW;
            h.NodeColor = 'k';
            h.Marker = 'o';
        end
    end

%% Static Methods
% These methods operate directly on lists of nodes and links.

    methods (Static)

        function joinchain(nodes, links)
        %JOINCHAIN Join a list of nodes with a list of links.
        %   JOINCHAIN(nodes, links) Attaches links to nodes in series.
        %   Supply N nodes and N-1 links.
        %
        %   See also ThermalNet.setcond.

            numN = numel(nodes);    % get the number of nodes

            % the first node's new neighbor is the second node in the list
            nodes(1).Neighbors =    [nodes(1).Neighbors nodes(2)];
            % the first node's new link is the first link in the list
            nodes(1).Links =        [nodes(1).Links links(1)];

            % for each internal node
            for i = 2:numN-1
                % the node's new neighbors are the nodes before and after
                nodes(i).Neighbors =    [nodes(i).Neighbors nodes([i-1 i+1])];
                % the node's new links are the links before and after
                nodes(i).Links =        [nodes(i).Links links([i-1 i])];
            end

            % the last node's new neighbor is the second-to-last node
            nodes(numN).Neighbors =     [nodes(numN).Neighbors nodes(numN-1)];
            % the last node's new link is the last path in the list
            nodes(numN).Links =         [nodes(numN).Links links(numN-1)];

            % the first path's nodes are the first and second nodes
            links(1).Nodes = nodes(1:2);

            % for each internal path
            for i = 2:numN-2
                % the path's nodes are the next two nodes in the list
                links(i).Nodes = nodes([i i+1]);
            end

            % the last path's nodes are the second-to-last and last nodes
            links(numN-1).Nodes = nodes([numN-1 numN]);
        end

        function setcond(links,K)
        %SETCOND Set the thermal conductances of a list of links.
        %   SETCOND(links,K) Apply thermal conductances K to links. Supply
        %   N links and N conductances, or N links and one conductance.
        %
        %   See also ThermalNet.joinchain.

            numL = numel(links);            % get the number of links
            numK = numel(K);                % get the number of conductances
            if numL == numK                 % if each link is given unique conductance
                for i = 1:numL                  % for each link
                    links(i).Conductance = K(i);    % set its conductance
                end
            else                            % else one conductance should be given
                for i = 1:numL                  % for each link
                    links(i).Conductance = K;       % set it to the single given conductance
                end
            end
        end

    end

end
