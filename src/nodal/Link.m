classdef Link < handle
%LINK A link to be used in a thermal network.
    %   LINK generates a link.
    %
    %   LINK(N) generates N links.

    properties
        Nodes Node                  % nodes attached to this link
        Conductance double = 0      % conductance of this link
        Flux double                 % heat flux through this link
    end
end
