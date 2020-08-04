function [ek] = GSDST_kEstimator(observation)

% Geometric Sequence Decomposition with k-simplxes Transform.
% Author: Woong-Hee Lee, Jong-Ho Lee, and Ki Won Sung
% Date: 2020-08-02

% Geometric Sequence Check

    upperk = ceil(length(observation)/2-1);
    val = zeros(1,upperk);
    
    for hatk=1:upperk
        
       % Make searchspace
       searchspace = [];
       for indexS = 1:hatk+2
          searchspace = [searchspace ; observation(indexS:indexS+hatk-1)]; 
       end
       
       % Make consecutive three simplexes and calculate volumes
       Vsimplexes = [];
       for indexS = 1:3
       Vsimplexes = [Vsimplexes det(searchspace(indexS:indexS+hatk-1,:))];
       end
       
       % Check the condition of geometric sequence

       val(hatk) = abs(Vsimplexes(2)/Vsimplexes(1) - Vsimplexes(3)/Vsimplexes(2));
       
    end
        
    [ss ek] = min(val);
    
end