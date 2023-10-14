classdef Node < handle
%NODE A node to be used in a thermal network.
    %   NODE generates a node.
    %
    %   NODE(N) generates N nodes.

    properties
        Neighbors Node          % nodes that neighbor this node
        Links Link              % links attached to this node
        Address uint16          % the latest address assigned to this node by a thermal network

        Temperature double      % the temperature of this node
        BoundHeat double        % heat input to this node
    end
end
