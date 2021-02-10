function plot_loss(horzs, results, add_line, plot_name, plot_legend, font_size, varargin)

    % Function for plotting results across estimation methods

    line_colors = repmat(lines(7),2,1);
    colors_indx = [3 4 1 2 5 6 7];
    line_colors = line_colors(colors_indx,:);
    line_styles = {'-', '--', '-x', '-', '-.', '-o', ':'};
    line_width  = repmat(2,14,1);
    line_width(1) = 5;
    line_width(4) = 5;
    
    if isempty(varargin)
        figure;
    end
    set(gca,'FontSize',7/8*font_size)
    set(gca,'TickLabelInterpreter','latex')
    grid on
    
    hold on;
    for i=1:size(results,2)
        plot(horzs, results(:,i)', line_styles{i}, 'Color', line_colors(i,:), 'LineWidth', line_width(i));
        xlim([min(horzs) max(horzs)])
    end
    hold off;
    
    if ~isempty(add_line)
        hold on;
        the_xlim = xlim;
        line(the_xlim, add_line*ones(1,2), 'Color', 'k', 'LineStyle', ':');
        hold off;
        xlim(the_xlim);
    end
    title(plot_name,'interpreter','latex','FontSize',9/8*font_size);
    xlabel('Horizon','interpreter','latex','FontSize',font_size);
    if isempty(varargin)
        legend(plot_legend, 'Location', 'southoutside', 'NumColumns', 3, 'interpreter', 'latex','FontSize',font_size);
    end

end