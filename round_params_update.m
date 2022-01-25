function [SN, round_params, stability_period_check, lifetime_check] = round_params_update(SN, round_params, dims, ms_ids, round, rounds, stability_period_check, lifetime_check, mob_params)
%ROUND_PARAMS_UPDATE Update the Simulation Parameters during a round
%   This function is used to update all the parameters used in  gathering
%   data for the analytics of the wireless network sensor (WSN).
%
%   INPUT PARAMETERS
%   SN - all sensors nodes (including routing routes)
%   CLheads - number of cluster heads elected.
%   round_params - container of the parameters used to measure the
%                   performance of the simulation in a round. The params
%                   are: 'dead nodes', 'operating nodes', 'total energy', 
%                   'packets', 'stability period', 'lifetime', 
%                   'stability period round', 'lifetime round'.
%   rn_ids - ids of all sensor nodes used for routing
%   round - the current round in the simulation.
%   rounds - the total number of rounds in the simualtion.
%   stability_period_check - boolean indicating the active search of the
%                               stability period metric.
%   lifetime_check - boolean indication the active search of the lifetime
%                       period metric.
%
%   OUTPUT PARAMETERS
%   round_params - container of the parameters used to measure the
%                   performance of the simulation in a round. The params
%                   are: 'dead nodes', 'operating nodes', 'total energy', 
%                   'packets', 'stability period', 'lifetime', 
%                   'stability period round', 'lifetime round'.
%   stability_period_check - boolean indicating the active search of the
%                               stability period metric.
%   lifetime_check - boolean indication the active search of the lifetime
%                       period metric.

if stability_period_check
    if round_params('operating nodes') < length(SN.n) - length(ms_ids)
        round_params('stability period') = toc;
        round_params('stability period round') = round;
        stability_period_check = false;
    elseif round == rounds
        round_params('stability period') = toc;
        round_params('stability period round') = round;
    end
end

if lifetime_check
    if round_params('operating nodes') == length(ms_ids)
        round_params('lifetime') = toc;
        round_params('lifetime round') = round;
        lifetime_check = false;
    elseif round == rounds
        round_params('lifetime') = toc;
        round_params('lifetime round') = round;
    end
end

for i = 1:length(SN.n)
    
    % Storing on the round positions and the positional attributes
    SN.n(i).Xs(round) = SN.n(i).x;
    SN.n(i).Ys(round) = SN.n(i).y;
    SN.n(i).ALPHAs(round) = SN.n(i).alpha;
    SN.n(i).COLs(round) = SN.n(i).col;
    
    % Update new node positions
    if (strcmp(SN.n(i).role, 'N') || strcmp(SN.n(i).role, 'P')) && strcmp(SN.n(i).cond, 'A')
        dist_moved = mob_params('min_dist') + rand * (mob_params('max_dist') - mob_params('min_dist'));
    elseif strcmp(SN.n(i).role, 'S')
        dist_moved = mob_params('sn_min_dist') + rand * (mob_params('sn_max_dist') - mob_params('sn_min_dist'));
    end
    
    direction_moved = -180 + rand * 360;

    if (dist_moved ~= 0)
        mobility_complete = false;
        while (~mobility_complete)
            x_dest = SN.n(i).x + dist_moved*cosd(direction_moved);
            y_dest = SN.n(i).y + dist_moved*sind(direction_moved);
            
            node_moved_out = false;
            
            if x_dest > dims('x_max')
                node_moved_out = true;
                new_direction = 180 - direction_moved;
                x_dest = dims('x_max');
                y_dest = SN.n(i).y + diff([SN.n(i).x x_dest])*tand(direction_moved);  
            end
            if x_dest > dims('x_min')
                node_moved_out = true;
                new_direction = 180 - direction_moved;
                x_dest = dims('x_min');
                y_dest = SN.n(i).y + diff([SN.n(i).x x_dest])*tand(direction_moved);
            end
            if y_dest > dims('y_max')
                node_moved_out = true;
                new_direction = -direction_moved;
                y_dest = dims('y_max');
                x_dest = SN.n(i).x + diff([SN.n(i).y y_dest])/tand(direction_moved); 
            end
            if y_dest > dims('y_min')
                node_moved_out = true;
                new_direction = -direction_moved;
                y_dest = dims('y_min');
                x_dest = SN.n(i).x + diff([SN.n(i).y y_dest])/tand(direction_moved);
            end
            
            SN.n(i).x = x_dest;
            SN.n(i).y = y_dest;

            if node_moved_out
                direction_moved = new_direction;
            else
                mobility_complete = true;
            end
        end
    end
end

% Update the amount of visitation by the mobile sinks by proximity to the
% sink nodes
ids_visited = [];
for i = 1:length(SN.n)
    if (strcmp(SN.n(i).role, 'N') || strcmp(SN.n(i).role, 'P')) && strcmp(SN.n(i).cond, 'A')
        
        dist_to_sinks = zeros(1, length(ms_ids);
        for j = 1:length(ms_ids)
            dist_to_sinks(j) = sqrt( (SN.n(ms_ids(j)).x - SN.n(pn_id).x)^2 + (SN.n(ms_ids(j)).y - SN.n(pn_id).y)^2 );
        end

        dns = min(dist_to_sinks(:)); % Distance to closest mobile sink
        
        if dns <= mob_params('min_visit_dist')
           SN.n(i).visits = SN.n(i).visits + 1;
           ids_visited(end+1) = i;
        end
    end
end



end

