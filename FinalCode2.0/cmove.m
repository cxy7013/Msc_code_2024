function [theta_x,theta_y] = cmove(I_turb,X,Y,R)
% Calculate centroid to determine beam wander
total_intensity = sum(I_turb(:));
x_centroid = sum(sum(X .* I_turb)) / total_intensity;
y_centroid = sum(sum(Y .* I_turb)) / total_intensity;

% Calculate the angular offset of the centroid position 
% (convert from meters to radians in 1 to urad)
theta_x = x_centroid / R * 1e6;  
theta_y = y_centroid / R * 1e6; 
end

