function [fit_solution, resnorm,residual, J]  = ExponentialR2starFit(te, data, fit_prop, b_plot)
%MYEXPONENTIALFIT Performs a multi-exponential fit ot gradient echo data
%   

te = te(:); 
data = data(:); 

data = abs(data); 

%perform fitting 
[fit_solution, resnorm,residual, J]...
    = fitting_routine_biexp_fit(te, data,fit_prop, b_plot);

end


function [fit_solution, resnorm,residual, J]...
    = fitting_routine_biexp_fit(t, data_points,fit_prop,b_plot)
%@param t              time vector
%@param data_points    Measured data points for each point in time

%error function - returns the residuals of the current parameters
fun =  @(x) min_fun_exp_fit(x, t, data_points);


%inital value and boundary conditions (these are optimal)
x0 = [1000,     0,      ]; % A / R2s_1 /
lb = [0,        0,      ]; 
ub = [30000,    1E4,     ]; 



%perform fitting 
opts1=  optimset('display','off','TolFun',1e-12,'TolX',1e-10);
[fit_solution, resnorm, residual,~,~,~, J] =...
    lsqnonlin(fun, x0,  lb, ub, opts1 );

%check if the values are reasonable 
if fit_solution(1) < 0.01
    fit_solution(1) = 0; fit_solution(2) = 0; 
end

if b_plot
    
    %string with the equation 
    eqn_str = [num2str(fit_solution(1),'%.0f'), ...
    ' exp(-t*(' num2str(fit_solution(2)*1E3,'%.0f'), 's^-1)', ')']; 
    
    
    %calcualte the fitted curve 
    t_fit = 0:t(end)/100:t(end); 
    f_fit = fit_solution(1)*exp(-t_fit.*fit_solution(2));
    
    if ~isfield(fit_prop, 'fig_mono_exp')
       fit_prop.fig_mono_exp = figure; 
    end

    figure(fit_prop.fig_mono_exp)
    scatter(t, data_points); 
    hold on; 
    plot(t_fit, f_fit); title('Exponetial Fit of Magnitude');  
    xlabel('time in ms '); ylabel('function value'); 
    legend('noisy data', 'fit-result'); 
%     set(fit_prop.figures.figure_biexp_fit, 'Name','BiExponentialFit','Position',...
%     [fit_prop.scrsz(3)*1/8 fit_prop.scrsz(4)*3/6 ...
%     fit_prop.scrsz(3)/4 fit_prop.scrsz(4)*2/6] );
    hold off; 
    annotation('textbox',[.5 .2 .6 .5],'String',eqn_str,'FitBoxToText','on')

end


end

function [ F ] = min_fun_exp_fit(x, t, data_points)
%Cost function for a biexponential fit. Returns the residuals for the given
%data points 
%@param x              Coefficients of the biexponetial function
%@param t              time vector
%@param data_points    Measured data points for each point in time

%parameters
A = x(1); 
R2s = x(2); 


%calculate the function value for the given parameters 
f = A*exp(-t.*R2s);

%return resdiuals 
F = f - data_points; 



end

