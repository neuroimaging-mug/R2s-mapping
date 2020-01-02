function [ output_args ] = savefigure( fig_name, fig,path,scale)
%% SACVEFIGURE saves figure
%@input fig_name    name of figure
%@input fig         figure handler of the current figure
%@input path        path within the figs folder is created
%%
%save figure
%get current path
    if scale == 1
    set(fig, 'Position', get(0,'ScreenSize'));
    end
    
%     set(fig, 'Position', [100, 100, 300, 300]);
    
%     path_ws = path;
%     folder_name = 'figs';
%     if (exist([path_ws,'/', folder_name]) == 0)
%         mkdir([path_ws,'/', folder_name]);
%     end
%     path_file = [path_ws, '/' folder_name '/' fig_name];

    
    path_ws = path;
    if (exist(path_ws) == 0)
        mkdir(path_ws);
    end
    path_file = [path_ws, '/' fig_name];
    
%     F_size = 18; 
% %     perform some adaptation to get a nice figure
%     set(gca,'FontSize', F_size);
%     hx=get(gca,'xlabel');
%     set(hx, 'FontSize', F_size);
%     hy=get(gca,'ylabel');
%     set(hy, 'FontSize', F_size);
      set(findall(fig, 'Type', 'Line'),'LineWidth',2);
%       set(findall(fig, 'Type', 'title'),'FontSize', 16);
%     t=get(gca,'title');    
% %     set(t, 'FontSize', 12);
    
    %to get a better area for the plot estmate the screen resolution
    %and set the figure to that size
    fig_name = [fig_name, '"'];
    %set(0, 'defaultTextInterpreter', 'latex');

    saveas(fig, path_file, 'epsc'); 
    saveas(fig, path_file, 'fig');
    saveas(fig, path_file, 'png'); 
    saveas(fig, path_file, 'svg'); 

end

