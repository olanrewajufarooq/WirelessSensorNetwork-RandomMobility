function plot_simulation(SN, rounds, dims)
%PLOT_SIMULATION Summary of this function goes here
%   Detailed explanation goes here

figure(1);
hold on;

xlabel('X (meters)');
ylabel('Y (meters)');
title('Radom Waypoint mobility');

plot( dims('x_min'),dims('y_min'),dims('x_max'),dims('y_max') );

hold off;

for i = 1:length(SN.n)
    node_plot(i) = plot(SN.n(i).Xs(1), SN.n(i).Ys(1), 'o' );
end

ht = text(dims('x_min'), dims('y_max'), cat(2,'Time (sec) = 0'));

for round = 1:rounds
    set(ht, 'String', cat(2,'Round = ', num2str(round)));
    for i = 1:length(SN.n)
        set(node_plot(i), 'XData', SN.n(i).Xs(round), 'YData', SN.n(i).Ys(round) );
    end
    drawnow;
end

end

