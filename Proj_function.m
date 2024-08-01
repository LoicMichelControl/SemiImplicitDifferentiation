    %=================================================================
    %=================================================================

        % (c) [2024]  Nantes Université - Centrale Nantes - LS2N UMR 6004, Nantes
        % (c) [2024]  Quartz EA 7393, ENSEA, Cergy-Pontoise
        % Loïc MICHEL and Jean-Pierre BARBOT
    
    %=================================================================
    %=================================================================
    
    function  [ E, Proj_ , lq] =  Proj_function( alpha, lambda, n, e_1, mu, h )
    
    pq = n*(1 - alpha );
    
    lq = lambda * (mu * h )^n;
    
    lq_ = (lambda * (mu * h )^n)^( 1 / pq );
    
    if  ( abs( e_1 ) )^( pq ) < lq
    
        Proj_ = ( ( abs( e_1 ) )^(pq) * sign( e_1 ) )/ lq; %
    
        E = 1;
    
    else
    
        Proj_ =  sign( e_1 );
    
        E = 0;
    
    end
    
    end
