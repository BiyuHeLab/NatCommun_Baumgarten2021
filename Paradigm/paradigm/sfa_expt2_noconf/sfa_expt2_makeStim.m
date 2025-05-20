function stim = sfa_expt2_makeStim(param, pahandle)

load series_selection_scale2Hz

% gather together the stimulus variables

ind = 0;
for i_b = 1:length(betas)
    for i_s = 1:3
        for i_tf = 1:6
            for i_rep = 1:2
                
                
                % larger sigma_e stim
                ind = ind + 1;

                stim.beta(ind)    = betas(i_b);
                stim.betaID(ind)  = i_b;
                
                stim.sigma_e(ind) = sigma_e(i_b);
                stim.sigma_e_ID(ind) = 1;
                               
                stim.logf_pen(ind)   = series_logHz{i_b}(i_s,end);
                
                stim.logf_final(ind) = series_logHz_tf_Hz{i_b}(i_s, i_tf);
                switch i_tf
                    case 1, stim.finalID(ind) = -3; % lowest value for final tone
                    case 2, stim.finalID(ind) = -2;
                    case 3, stim.finalID(ind) = -1;
                    case 4, stim.finalID(ind) =  1;
                    case 5, stim.finalID(ind) =  2;
                    case 6, stim.finalID(ind) =  3; % highest value for final tone
                end
                
                stim.logf_pred(ind)  = series_pred_logHz{i_b}(i_s);
                switch i_s
                    case 1, stim.predID(ind) = -1; % min predicted frequency
                    case 2, stim.predID(ind) = +1; % max predicted frequency
                    case 3, stim.predID(ind) =  0; % median predicted frequency
                end
                
                stim.series_f{ind}  = [series_Hz{i_b}(i_s,:) series_Hz_tf_Hz{i_b}(i_s, i_tf)];
                stim.soundwave{ind} = series2soundwave(stim.series_f{ind}, param.toneDur_inSecs, param.f_sample_inHz);
                

                % matched, smaller sigma_e stim
                ind = ind + 1;
                
                stim.beta(ind)    = betas(i_b);
                stim.betaID(ind)  = i_b;                
                
                stim.sigma_e(ind) = sigma_e_sc(i_b);
                stim.sigma_e_ID(ind) = 0;
                               
                stim.logf_pen(ind)   = series_logHz_sc{i_b}(i_s,end);
                
                stim.logf_final(ind) = series_logHz_sc_tf_Hz{i_b}(i_s, i_tf);
                switch i_tf
                    case 1, stim.finalID(ind) = -3; % lowest value for final tone
                    case 2, stim.finalID(ind) = -2;
                    case 3, stim.finalID(ind) = -1;
                    case 4, stim.finalID(ind) =  1;
                    case 5, stim.finalID(ind) =  2;
                    case 6, stim.finalID(ind) =  3; % highest value for final tone
                end
                
                stim.logf_pred(ind)  = series_pred_logHz_sc{i_b}(i_s);
                switch i_s
                    case 1, stim.predID(ind) = -1; % min predicted frequency
                    case 2, stim.predID(ind) = +1; % max predicted frequency
                    case 3, stim.predID(ind) =  0; % median predicted frequency
                end
                
                stim.series_f{ind} = [series_Hz_sc{i_b}(i_s,:) series_Hz_sc_tf_Hz{i_b}(i_s, i_tf)];                
                stim.soundwave{ind} = series2soundwave(stim.series_f{ind}, param.toneDur_inSecs, param.f_sample_inHz);

            end
        end
    end
end


% randomize the stimuli
if ~param.restarting
    rand_ind = randperm(ind);
else
    rand_ind = param.stim_ind_order;
end


stim.beta       = stim.beta(rand_ind);
stim.betaID     = stim.betaID(rand_ind);
stim.sigma_e    = stim.sigma_e(rand_ind);
stim.sigma_e_ID = stim.sigma_e_ID(rand_ind);
stim.logf_pen   = stim.logf_pen(rand_ind);
stim.logf_final = stim.logf_final(rand_ind);
stim.finalID    = stim.finalID(rand_ind);
stim.logf_pred  = stim.logf_pred(rand_ind);
stim.predID     = stim.predID(rand_ind);
stim.series_f   = stim.series_f(rand_ind);
stim.soundwave  = stim.soundwave(rand_ind);
stim.sigma_e_w  = sigma_w;
stim.nTrials    = ind;
stim.nBlocks    = stim.nTrials / param.nTrialsPerBlock;
stim.pahandle   = pahandle;
stim.ind_order  = rand_ind;

end