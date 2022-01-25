function plot_simulation(SN, rounds, dims)
%PLOT_SIMULATION Summary of this function goes here
%   Detailed explanation goes here

nodes = 1:rounds;

for i = 1:length(SN.n)
    %Simple interpolation (linear) to get the position, anytime.
    %Remember that "interp1" is the matlab function to use in order to
    %get nodes' position at any continuous time.
    node(i).x = interp1(nodes,SN.n(i).Xs,5);
    node(i).y = interp1(nodes,SN.n(i).Ys,5);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);
hold on;

xlabel('X (meters)');
ylabel('Y (meters)');
title('Radom Waypoint mobility');

plot( dims('x_min'),dims('y_min'),dims('x_max'),dims('y_max') );

for i = 1:length(SN.n)
    vh_node_pos(i) = plot(SN.n(i).Xs(1), SN.n(i).Ys(1), 'o', 'MarkerFaceColor', [SN.n(i).COLs(1,:)] );
end

ht = text(min(node(1).x),max(node(1).y),cat(2,'Time (sec) = 0'));
axis([min(node(1).x) max(node(1).x) min(node(1).y) max(node(1).y)]);

hold off;
for round = 1:length(nodes);
    set(ht, 'String', cat(2,'Round = ',num2str(round)));
    for i = 1:length(SN.n)
        set(vh_node_pos(i),'XData',node(i).x(round),'YData',node(i).y(round), node(i).y(round, :));
    end
    drawnow;
end

end

