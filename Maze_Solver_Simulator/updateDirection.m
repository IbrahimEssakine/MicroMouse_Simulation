function new_direction = updateDirection(previous_point, new_point, previous_direction)
    % Calculate the change in each dimension
    delta_x = new_point(1) - previous_point(1);
    delta_y = new_point(2) - previous_point(2);
    
    % Update the direction based on the change in each dimension
    if delta_x > 0 && delta_y > 0
        new_direction = [1, 1]; % Move right && up
    elseif delta_x < 0 && delta_y > 0
        new_direction = [-1, 1]; % Move left && up
    elseif delta_x < 0 && delta_y < 0
        new_direction = [-1, -1]; % Move left && down
    elseif delta_x > 0 && delta_y < 0
        new_direction = [1, -1]; % Move right && down
    elseif delta_x < 0
        new_direction = [1, 0]; % Move Right
    elseif delta_x < 0
        new_direction = [-1, 0]; % Move left
    elseif delta_y > 0
        new_direction = [0, 1]; % Move up
    elseif delta_y < 0
        new_direction = [0, -1]; % Move down
    else
        % No change in position
        new_direction = previous_direction;
    end
end


