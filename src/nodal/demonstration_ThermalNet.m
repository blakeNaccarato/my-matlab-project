% This is a non-functioning demo of the thermal network program w/ function bodies snipped out

classdef ThermalNet < handle
%THERMALNET A thermal network of nodes linked by conductive paths.
    properties
        Nodes Node              % list of nodes that are in the network
        Links Link              % list of paths that are in the network

        Adjacency logical       % adjacency array used in network plotting
        Conductances double     % conductance array used in solving the network
        Fluxes double           % flux array to be computed after solving for temperatures

        Knowns uint16           % keeps track of boundary conditions set on the network
    end

    methods
        function net = ThermalNet(varargin)  % Constructor

        % Networking
        function branch(net, rootNodeID, numN, varargin) % Generate a sidechain
        function join(net, rootNodeID, joinNodeIDs, varargin) % Join nodes to  a root node
        function netjoin(rootNet, rootNodeID, joinNet, joinNodeID) % Join two networks
        function refine(net, nodePairID, multiple) % Create intermediate nodes
        function reduce(net, nodePairID) % Destroy intermediate nodes
        function address(net) % Give each node a unique address
        function update(net) % Update adjacency and conductance arrays

        % Set BC's and solve
        function addbound(net, type, nodeID, boundvalue) % Add a boundary condition
        function solve(net) % Solve for nodal temperatures
        function getfluxes(net) % Obtain heat fluxes

        % Visualize
        function labels = gettemplabels(net) % Get node labels
        function [h G] = plot(net, type) % Plot the network
    end

    % More networking
    methods (Static)
        function joinchain(nodes, links) % Join a list of nodes with a list of links.
        function setcond(links,K) % Set the thermal conductances of a list of links.
    end
end
