function [] = plot_surface_strain(s, lams, lamp, fig_num)
    % PLOT_SURFACE_STRAIN Plots the surface strain components.
    %
    % INPUTS:
    %   s       - Arclength coordinates.
    %   lams    - Meridional stretch ratio (\lambda_s).
    %   lamp    - Circumferential stretch ratio (\lambda_\phi).
    %   fig_num - Figure number for plotting.

    % plot the surface stresses
    figure(fig_num); hold on
    h(1) = plot(s,lams,'LineWidth',2); hold on
    h(2) = plot(s,lamp,'LineWidth',2);
    xlabel('s','FontSize',32);
    ylabel('\lambda','FontSize',32);
    legend(h,'\lambda_s','\lambda_\phi','FontSize',24,'Location','northwest');
    xlim([0,s(end)])
    ax = gca; ax.FontSize = 24;

end
