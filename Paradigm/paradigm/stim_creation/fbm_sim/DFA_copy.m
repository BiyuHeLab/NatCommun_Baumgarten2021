function [alpha, pval] = DFA_copy(data, doPlot)
        if nargin < 2, doPlot = 1; end

        %%% summation %%%
        x = cumsum(data);
        %%% DFA %%%
        l = [5 10 19 38 95];     %segment length (need to tune to data)
        n = length(x) ./ l;      %number of segments
  
        F = zeros(1,length(l));
        for i = 1:length(l)
            for k = 1:n(i)
                x0 = detrend(x((k-1)*l(i)+1:k*l(i)));
                F(i) = F(i) + std(x0,1);
            end
            F(i) = F(i) / n(i);
        end
           
        if doPlot
            figure(1); clf; loglog(l,F,'.-','MarkerSize',15); title('DFA')  
            grid on; xlabel('window length'); ylabel('F');   %slope is alpha, alpha = H for fGn
        end
        
        stat = regstats(log(F)',log(l)','linear');
        alpha = stat.beta(2);
        pval = stat.tstat.pval(2);
end

