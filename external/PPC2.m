function ppc_out = PPC2(SpikePhases,SpikeTrialNums)
% Calculate point-field Pairwise Phase Consistency P^_2 as described by 
% Vinck et al. (2012, J Comp Neurosci) for one unit.

% angle_dot_prod = @(theta1,theta2) dot([cos(theta1),sin(theta1)], [cos(theta2),sin(theta2)]);

% convert inputs to column vectors (if they aren't already)
if size(SpikePhases,2) ~= 1
    SpikePhases = SpikePhases';
end

uniqueTrialNums = unique(SpikeTrialNums);
n_trials = length(uniqueTrialNums);

avg_dot_prods = nan(n_trials); % will contain the average of dot products across spike phases for each pair of trials

for m_ind = 1:n_trials
    m = uniqueTrialNums(m_ind);
    thetas_m = SpikePhases(SpikeTrialNums == m);
    for l_ind = 1:n_trials
        l = uniqueTrialNums(l_ind);
        if l ~= m
            thetas_l = SpikePhases(SpikeTrialNums == l);
    
            % [length(thetas_m) x length(thetas_l)] array of dot products between each pair of spike-LFP phase vectors
            % dot_prods_ml = bsxfun(@angle_dot_prod, thetas_m, thetas_l);
            % avg_dot_prods(m_ind,l_ind) = mean(dot_prods_ml,'all');
    
            % vectorized computation of dot product for each pair of spike-LFP
            % phase vectors
            % avg_dot_prods(m_ind,l_ind) = mean(cos(bsxfun(@minus,thetas_m,thetas_l')), 'omitnan');
            avg_dot_prods(m_ind,l_ind) = mean(cos(thetas_m - thetas_l'), 'all', 'omitnan');

        end
    end
end

ppc_out = mean(avg_dot_prods, 'all', 'omitnan');

end

% function dotprod = angle_dot_prod(theta1, theta2)
%     % gets the dot product of two angles' vector representations.
%     % INPUTS: theta1, theta2 are angles in radians
%     dotprod = dot([cos(theta1) sin(theta1)], [cos(theta2) sin(theta2)]);
% end