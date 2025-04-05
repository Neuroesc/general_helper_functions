function [z,p,q,obs,shuff] = stats_perm_test(v1,v2,varargin)
%stats_perm_test permutation test for two groups 
% Permutation test, with or without replacement, for two groups
% A number of different statistical measures can be selected
%
%   [z,p,obs,shuff] = stats_perm_test(v1,v2) compute permutation z-test
%   and probability of the difference between the two groups. Default
%   measure is mean difference.
%
%   [z,p,obs,shuff] = stats_perm_test(__,name,value) process with Name-Value pairs 
%
%   parameters include:
%
%   'method'        -   String, method to use when testing difference between groups.
%                       Select from: 'mean' (mean difference), 'median' (median difference),
%                       't' (t-statistic), 'f' (F-statistic), 'cohen' (Cohen's d), 'hedge' (Hedge's g),
%                       or 'p_super' (probability of superiority).
%
%                       Default value is 'mean'.
%
%   'iti'           -   Scalar, number of iterations for permutation.
%
%                       Default value is 1e3 or 1000.
%
%   'rep'           -   Scalar, permutations should be with (1) or without (0) replacement.
%
%                       Default value is 1.
%
%   outputs include:
%
%   'z'             -   Float, z-score for observed group difference relative to the shuffles
%
%   'p'             -   Float, two-tailed p-value associated with z-value.
%
%   'q'             -   [2 x 1], left and right one-tailed p-value associated with z-value 
%                       respectively
%
%   'obs'           -   Float, the observed difference between groups, value will depend on the 
%                       method chosen for calculating group differences. i.e. if 'mean' method
%                       is chosen this will be the observed differences in group means.
%
%   'shuff'         -   [iti x 1], vector of group differences as for 'obs' but for each shuffle.
%
%
%   Class Support
%   -------------
%   The input matrix in must be a real, non-sparse matrix of
%   the following classes: uint8, int8, uint16, int16, uint32, int32,
%   single or double.
%
%   Notes
%   -----
%   1. Group differences are calculated given the group order provided to the function. i.e.
%       a negative mean difference between groups will occur when v2 has a larger mean
%       value than v1.
%
%   Example
%   ---------
% 
%   See also NAME, NAME

% HISTORY:
% version 1.0.0, Release 18/11/23 Initial release
%
% Author: Roddy Grieves
% Dartmouth College, Moore Hall
% eMail: roddy.m.grieves@dartmouth.edu
% Copyright 2023 Roddy Grieves

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Heading 3
%% >>>>>>>>>>>>>>>>>>>> Heading 2
%% >>>>>>>>>> Heading 1
%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INPUT ARGUMENTS CHECK
%% Parse inputs
    p = inputParser;
    addRequired(p,'v1',@(x) ~isempty(x) && ~all(isnan(x(:))));  
    addRequired(p,'v2',@(x) ~isempty(x) && ~all(isnan(x(:))));  
    expectedmethods  = {'mean','median','t','f','cohen','hedge','cliff','p_super'};        
    addParameter(p,'method','mean',@(x) any(validatestring(x,expectedmethods)));  
    addParameter(p,'iti',1e3,@(x) ~isempty(x) && ~all(isnan(x(:))) && isscalar(x) );  
    addParameter(p,'rep',1,@(x) isscalar(x) );  
    parse(p,v1,v2,varargin{:});
    config = p.Results;
    v1 = config.v1(:); % group 1
    v2 = config.v2(:); % group 2

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FUNCTION BODY
%% >>>>>>>>>> Observed values
    obs = group_difference(v1,v2,config.method); % observed group difference value
    v3 = [v1;v2]; % combined values
    n = [numel(v1) numel(v2) numel(v3)]; % group sizes

%% >>>>>>>>>> Shuffled values
    rng(999); % for reproducibility
    shuff = NaN(config.iti,1);
    for ii = 1:config.iti
        if config.rep
            rindx = randi(n(3),n(3),1); % random order of all elements with replacement
        else
            rindx = randperm(n(3)); % random order of all elements without replacement            
        end
        r3 = v3(rindx); % shuffled data
        r1 = r3(1:n(1)); % shuffled group 1
        r2 = r3(n(1)+1:end); % shuffled group 2
        shuff(ii) = group_difference(r1,r2,config.method); % shuffled group difference value
    end

%% >>>>>>>>>> Significance
    [z,p,q] = stats_z_probability(obs,shuff);

end

%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Difference function
function d = group_difference(g1,g2,method)
    switch method
        case {'mean'} % mean difference
            d = mean(g1,'all','omitnan') - mean(g2,'all','omitnan');
    
        case {'median'} % median difference
            d = median(g1,'all','omitnan') - median(g2,'all','omitnan');
    
        case {'t'} % t statistic
            [~,~,~,s] = ttest2(g1,g2);    
            d = s.tstat;

        case {'f'} % F statistic
            [~,s,~] = anova1([g1;g2],[ones(size(g1)); ones(size(g2)).*2],'off');
            d = s{2,5};

        case {'cohen'} % Cohen's d
            e = stats_effect_size(g1,g2);
            d = e.cohen_d;

        case {'hedge'} % Hedge's g
            e = stats_effect_size(g1,g2);
            d = e.hedge_g;

        case {'cliff'}
            e = stats_effect_size(g1,g2);
            d = e.cliff_delta;

        case {'p_super'} % Probability of superiority
            e = stats_effect_size(g1,g2);
            d = e.probability_of_superiority;

        otherwise
            error('Unknown test type... exiting')
    end
end





























