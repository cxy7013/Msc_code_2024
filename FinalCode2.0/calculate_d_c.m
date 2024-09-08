function [R86, centroid] = calculate_d_c(I, X, Y)
    % % Calculate centroid
    I_sum = sum(I(:));
    centroid_x = sum(sum(X .* I)) / I_sum;
    centroid_y = sum(sum(Y .* I)) / I_sum;

    centroid(1)=centroid_x;
    centroid(2)=centroid_y;

    % Calculate distance relative to centroid
    R = sqrt((X - centroid_x).^2 + (Y - centroid_y).^2);

    % Sort distances
    [R_sorted, idx] = sort(R(:));
    I_sorted = I(idx);

    % Calculate the cumulative power distribution
    cumulative_power = cumsum(I_sorted);

    % Find the index containing 86% of the total power
    total_power = cumulative_power(end);
    idx_86 = find(cumulative_power >= 0.86 * total_power, 1);

    %
    R86 = R_sorted(idx_86);
    D86 = 2 * R86;
end
