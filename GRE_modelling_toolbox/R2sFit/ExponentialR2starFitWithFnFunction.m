function [fit_solution, resnorm,residual, J, output]  = ExponentialR2starFitWithFnFunction(te, data, Fn, fit_prop, b_plot)
%MYEXPONENTIALFIT Fits R2* star and A from S(t) = A*exp(-t*R2s)*Fn(t).
%   @param     te       Echo times in ms
%   @param     data     Data with length(te) 
%   @param     Fn       Signal decay due to macroscopic fields (S =
%                       A*exp(-R2s)*Fn)
%   Author: Martin Soellradl
%   Department of Neurology, Medical University of Graz, Graz, Austria
%   email:martin.soellradl@medunigraz.at
%   Website: http://www.neuroimaging.com
%   November 2019; Last revision: 01-November-2019


    te = te(:); 
    data = abs(data(:)); 
    Fn = abs(Fn(:)); 

    %perform fitting 
    [fit_solution, resnorm,residual, J, output]...
        = fitting_routine_exp_fit(te, data, Fn, fit_prop, b_plot);

end


function [fit_solution, resnorm,residual, J, output]...
    = fitting_routine_exp_fit(te, data, Fn, fit_prop, b_plot)
%@param t              time vector
%@param data_points    Measured data points for each point in time

    %error function - returns the residuals of the current parameters
    fun =  @(x) min_fun_exp_fit(x, te, Fn, data);


    %inital value and boundary conditions (these are optimal)
    %    [A                R2s(1/ms)    ]
    x0 = [data(1),         40*1E-3,     ]; % T2s
    lb = [data(1).*0.5,    0,           ]; 
    ub = [data(1).*1.5,    1000         ]; 

    %perform fitting 
    opts1=  optimset('display','off','TolFun',1e-12,'TolX',1e-10, 'Jacobian','on'); %

    [fit_solution, resnorm, residual,~,output,~, J] =...
        lsqnonlin(fun, x0, lb, ub,  opts1); 


    if b_plot
        A = fit_solution(1); 
        R2s_fit = fit_solution(2); 
        %string with the equation 
        eqn_str = [num2str(A,'%.0f'), ...
        ' exp(-t(' num2str(R2s_fit*1E3,'%.0f'), 's^{-1})', ') Fn']; 

        %calcualte the fitted curve 
        t_fit = 0:te(end)/100:te(end); 
        f_fit = abs(A*exp(-t_fit*R2s_fit));

        if ~isfield(fit_prop, 'fig_mono_exp')
           fit_prop.fig_mono_exp = figure; 
        end

        figure(fit_prop.fig_mono_exp)
        scatter(te, data); 
        hold on; 
        plot(t_fit, f_fit); title('Exponetial Fit of Magnitude');  
        xlabel('time in ms '); ylabel('function value'); 
        legend('noisy data', 'fit-result'); 
        hold off; 
        annotation('textbox',[.5 .2 .6 .5],'String',eqn_str,'FitBoxToText','on')
    end

end

function [F, J ] = min_fun_exp_fit(x, t, Fn,data_points)
%Cost function for a biexponential fit. Returns the residuals for the given
%data points 
%@param x              Coefficients of the biexponetial function
%@param t              time vector
%@param data_points    Measured data points for each point in time

    %parameters
    A = x(1); 
    R2s = x(2); 

    %calculate the function value for the given parameters 
    f = abs(A*exp(-t*R2s)).*Fn; %divide by pi! sinc in Matlab
    %defined by sinc = sin(pi*t)/(pi*t)

    %return resdiuals 
    F = f - data_points; 

    %estiamte jacobian 
    if nargout > 1
        J = zeros(length(t), 2); 
        J(:,1) = exp(-x(2)*t).*Fn;
        J(:,2) = - x(1)*t.*exp(-x(2)*t).*Fn;
    end

end

