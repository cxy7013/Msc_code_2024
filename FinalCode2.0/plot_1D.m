function plot_1D(x,I)

mid_row = round(size(I, 1) / 2);
cross_section = I(mid_row, :); 
cs_normalized = cross_section / max(cross_section);

plot(x, cs_normalized, 'k', 'LineWidth', 1.5);
axis tight;
xlabel('x / m');
ylabel('Normalized Intensity');
title('1-D Intensity Distribution');

end