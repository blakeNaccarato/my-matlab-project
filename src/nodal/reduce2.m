function reduce2(net, nodePairID)
    figure,
    [~, G] = plot(net,'node');
    close

    nodeIDs = shortestpath(G,nodePairID(1),nodePairID(2));
    nodes = net.Nodes(nodeIDs);

    numN = numel(nodes);
    numL = numN-1;

    cond = zeros(1,numL);
    flux = zeros(1,numL);

    i = 1;

    rootNode = nodePairID(1)-1;

    rtLink = intersect(nodes(i).Links, net.Nodes(rootNode).Links);
    [link ind1 ind2] = intersect(nodes([i i+1]).Links);

    rtLinkInf = setdiff(nodes(i).Links, [rtLink link]);
    rtInfNode = setdiff(rtLinkInf.Nodes, nodes(i));
    rtLinkInf.Conductance = numL*rtLinkInf.Conductance;
    fluxInf = rtLinkInf.Flux;

    cond(i) = link.Conductance;
    flux(i) = link.Flux;
    delete(link)
    nodes(i).Links(ind1) = [];
    nodes(i+1).Links(ind2) = [];

    for i = 2:numL
        [link ind1 ind2] = intersect(nodes([i i+1]).Links);

        [linkInf indInf] = setdiff(nodes(i).Links, link);
        infNode = setdiff(linkInf.Nodes, nodes(i));
        fluxInf = fluxInf + linkInf.Flux;
        delete(linkInf)
        nodes(i).Links(indInf) = [];
        rembound(net,infNode.Address)
        delete(infNode)

        cond(i) = link.Conductance;
        flux(i) = link.Flux;
        delete(link)
        nodes(i).Links(ind1) = [];
        nodes(i+1).Links(ind2) = [];
    end

    rembound(net,rtInfNode.Address)
    addbound(net,'heat',rtInfNode.Address,-fluxInf)

    net.Links(~isvalid(net.Links)) = [];
    equivCond = sum(1./cond)^-1;

    for i = 2:numN-1
        rembound(net,nodeIDs(i))
        delete(nodes(i))
    end
    net.Nodes(~isvalid(net.Nodes)) = [];

    join(net,nodePairID(1),nodePairID(2),equivCond)

    link = intersect(nodes(1).Links,nodes(end).Links);
    if sign(flux(1)) ~= sign(flux(end))
        equivFlux = -flux(end);
    else
        equivFlux = flux(end);
    end
    link.Flux = equivFlux;

    nodes( 1 ).Neighbors(~isvalid(nodes( 1 ).Neighbors)) = [];
    nodes(end).Neighbors(~isvalid(nodes(end).Neighbors)) = [];
end
