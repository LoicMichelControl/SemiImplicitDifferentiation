        
        %=================================================================
        %=================================================================
        
        % (c) [2024]  Nantes Université - Centrale Nantes - LS2N UMR 6004, Nantes
        % (c) [2024]  Quartz EA 7393, ENSEA, Cergy-Pontoise
        % Loïc MICHEL and Jean-Pierre BARBOT
        % All rights reserved under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International.
        
        % Code associated to the paper
        % "A note about high-order semi-implicit differentiation:
        % application to a numerical integration scheme with Taylor-based
        % compensated error" - preprint arXiv - Aug. 2024.
        
        %=================================================================
        %=================================================================
                
        clear all
        close all
        clc
        
        % SIMULATION SETTINGS 
        N = 2e5;
        Deltah = 1e-3;        % Same Time-step of the schemes
        

        % TO BEGIN: uncomment the Example ## ('f', 'y_exacth', 'IC_0', 'time0') you want to solve: 
        
            % Example #1 - linear ODE
               % f = @(t,y)(1 - y);
               % y_exacth = @(time)( 1 - exp( -time ) );
               % IC_0 = 1e-2;
               % time0 = 0;
        
            % Example #2 - logistic equation : taken from Griffiths, David & Higham, Desmond. (2010).
            % Numerical Methods for Ordinary Differential Equations:
            % Initial Value Problems. (p. 25)
                f = @(t,y)( 2*y*(1 - y) );
                y_exacth = @(time)( 1./(1 + 4*exp( 2*(10-time ))) );
                IC_0 = 1/5;
                time0 = 10;

            % Example #3 - taken from : Vuik, C. & Vermolen, F. & Gijzen, M.B. & Vuik, Mathea. (2023).
            % Numerical Methods for Ordinary Differential Equations. (p. 72)
                 %f = @(t,y)( -10*y.^2 + 20);
                 %y_exacth = @(time)( sqrt(2) * tanh(sqrt(300) * time ));
                 %IC_0 = 1e-2;
                 %time0 = 0;
        
        % I.C. for each solveur (Euler, RK & NSFD-SIHD3)       
        mz1_d = IC_0;
        mz2_d = 0;
        mz3_d = 0;
        mz4_d = 0;
        y_Euler = IC_0;
        y_RK = IC_0;
        y_SIHD = IC_0;
        
        % SIHD-3 parameters
        m_alpha = 1;
        m_alpha_1 = m_alpha * 0.95;
        m_alpha_2 = m_alpha * 0.95;
        m_alpha_3 = m_alpha * 0.95;
        m_alpha_4 = m_alpha * 0.95;
        
        m_lambda = 1;
        m_lambda_1 = m_lambda * 1e3;
        m_lambda_2 = m_lambda * 1e6;
        m_lambda_3 = m_lambda * 1e9;
        m_lambda_4 = m_lambda * 1e9;
        
        m_MU_ = 1;
        
        % to be investigated in the future -> should be set currently to zero
        m_theta_1 = 0;
        m_theta_2 = 0;
        m_theta_3 = 0;
        m_theta_4 = 0;
        m_lambda_p1 = abs(m_lambda_1 / (1 - m_theta_1));
        m_lambda_p2 = abs(m_lambda_2 / (1 - m_theta_2));
        m_lambda_p3 = abs(m_lambda_3 / (1 - m_theta_3));
        m_lambda_p4 = abs(m_lambda_4 / (1 - m_theta_4));
        
        % SIHD Proj. init.
        mE_1_m = 0;
        mE_2_m = 0;
        mE_2_mm = 0;
        mE_3_m = 0;
        mE_3_mm = 0;
        mE_3_mmm = 0;
        
        m_correction_factor = 1; % Enable the correction terms (enable = 1)
        for ii = 1:N
        
            % Time update
            time = time0 + Deltah * ii;
            time_vec(ii) = time;
        
            % Exact solution
            y_exact(ii) = y_exacth(time);
        
            % ======== Euler Forward
            y_Euler(ii + 1) = y_Euler(ii) + Deltah * f(time,y_Euler(ii));
        
            % ======== Runge-Kutta
            k_1 = f(0,y_RK(ii));
            k_2 = f(0+0.5*Deltah,y_RK(ii)+0.5*Deltah*k_1);
            k_3 = f((0+0.5*Deltah),(y_RK(ii)+0.5*Deltah*k_2));
            k_4 = f((0+Deltah),(y_RK(ii)+k_3*Deltah));
        
            y_RK(ii+1) = y_RK(ii) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*Deltah;
        
            % ======== SIHD-3 section ========
        
            % Compute the error e_1
            me_1 = f(time, y_SIHD(ii)) - mz1_d;
        
            % Compute the projectors
            [ mE_1, mProj_1_ ] = Proj_function( m_alpha_1, m_lambda_p1, 1, me_1, m_MU_, Deltah);
            [ mE_2, mProj_2_ ] = Proj_function( m_alpha_2, m_lambda_p2, 2, me_1, m_MU_, Deltah);
            [ mE_3, mProj_3_ ] = Proj_function( m_alpha_3, m_lambda_p3, 3, me_1, m_MU_, Deltah);
            [ mE_4, mProj_4_ ] = Proj_function( m_alpha_4, m_lambda_p4, 4, me_1, m_MU_, Deltah);
        
            % Update the differentiation
            mz4_p = mz4_d + (Deltah * mE_1_m * mE_2_m * mE_3 * m_lambda_4 * ( ( m_MU_  )^4 ) * ( abs(  me_1 ) )^(4 * m_alpha_4 - 3 ) * mProj_4_ );
        
            mz3_p = mz3_d + (Deltah * mE_1_m * mE_2_m * ( mz4_p + m_lambda_3 * ( ( m_MU_  )^3 ) * ( abs(  me_1 ) )^(3 * m_alpha_3 - 2 ) * mProj_3_ ) );
        
            mz2_p = mz2_d + Deltah * mE_1_m * (  mz3_p - ((m_correction_factor * 1/2  * Deltah * mz4_p )) + m_lambda_2 * mE_1_m * ( ( m_MU_  )^2 ) * ( abs(  me_1 ) )^(2 * m_alpha_2 - 1 ) * mProj_2_ );
        
            mz1_p = mz1_d + Deltah * ( mz2_p - (m_correction_factor * 1/2  * Deltah * mz3_p ) + ((m_correction_factor * 1/6 * Deltah^2 * mz4_p )) + m_lambda_1 * ( ( m_MU_  ) ) * ( abs( me_1 ) )^m_alpha_1 * mProj_1_  );
        
            % Proj. E_1 - update
            mE_1_m = mE_1;
        
            % Proj. E_2 - update
            mE_2_m = mE_2_mm;
            mE_2_mm = mE_2;
        
            % Proj. E_3 - update
            mE_3_m = mE_3_mm;
            mE_3_mm = mE_3_mmm;
            mE_3_mmm = mE_3;
        
            % States update
            mz1_d = mz1_p;
            mz2_d = mz2_p;
            mz3_d = mz3_p;
            mz4_d = mz4_p;
        
            % Solution update
            y_SIHD(ii + 1) = y_SIHD(ii) + Deltah * mz1_d;
            % Here Deltah can be replaced by something more complex in
            % accordance with the NSFD rules.
        
            % Vectorization
            mz1_d_vec(ii) = mz1_d;
            mz2_d_vec(ii) = mz2_d;
            mz3_d_vec(ii) = mz3_d;
            mz4_d_vec(ii) = mz4_d;
        
            x_1_vec_d(ii) = mz1_d;
            x_2_vec_d(ii) = mz2_d;
            x_3_vec_d(ii) = mz3_d;
            x_4_vec_d(ii) = mz4_d;
        
            % ======== End of the SIHD-3 section ========
        
            % Vectorization of the errors
            error_Euler(ii) = y_exact(ii) - y_Euler(ii);
            error_SIHD(ii) = y_exact(ii) - y_SIHD(ii);
            error_RK(ii) = y_exact(ii) - y_RK(ii);
        
        end
        
        for index = 1:100 % Print for example the 100 first error values
        
            fprintf('Euler = %e - SIHD = %e - RK = %e \n', error_Euler(index), error_SIHD(index), error_RK(index))
        
        end
        
        fprintf('Last errors at t = %f : \n', time)
        fprintf('Euler = %e - SIHD = %e - RK = %e \n', error_Euler(end), error_SIHD(end), error_RK(end))
        
        
        figure('name','ODE solver')
        plot(time_vec, y_Euler(1:end-1), 'b','linewidth', 3)
        hold on
        plot(time_vec, y_SIHD(1:end-1), 'r', 'linewidth', 3);
        plot(time_vec, y_RK(1:end-1), '--c','linewidth', 3);
        plot(time_vec, y_exact,'--k','linewidth', 3);
        legend('Euler', 'SIHD3', 'RK', 'Exact solution','Interpreter','latex', 'FontSize', 40);
        grid on
        xlabel('Time (s)','Interpreter','latex', 'FontSize', 40);
        ylabel('$y $','Interpreter','latex', 'FontSize', 40);
        set(gcf,'Color','w');
        set(gca,'FontSize',40);
        
        
        
        figure('name','ODE solver error')
        semilogy(time_vec, abs(error_Euler), 'b', 'linewidth', 3)
        hold on
        semilogy(time_vec, abs(error_SIHD), 'r','linewidth', 3)
        semilogy(time_vec, abs(error_RK), '--c', 'linewidth', 3)
        grid on
        set(gca, 'YScale', 'log')
        grid on
        xlim([0, 200]);
        xlabel('Time (s)','Interpreter','latex', 'FontSize', 40)
        ylabel('$|$error$|$','Interpreter','latex', 'FontSize', 40)
        legend('Euler', 'SIHD3', 'RK','Interpreter','latex', 'FontSize', 40)
        set(gcf,'Color','w');
        set(gca,'FontSize',40);
        
        
        figure('name','SIHD3 states')
        subplot(411)
        hold on
        plot( time_vec(1:end), x_1_vec_d(1:end), 'r', 'linewidth', 3);
        grid on
        ylabel('$x_1$','Interpreter','latex')
        grid on
        set(gca,'FontSize',40)
        subplot(412)
        hold on
        plot( time_vec(1:end), x_2_vec_d(1:end), 'r', 'linewidth', 3)
        grid on
        ylabel('$x_2$','Interpreter','latex')
        grid on
        set(gca,'FontSize',40)
        subplot(413)
        hold on
        plot( time_vec(1:end), x_3_vec_d(1:end), 'r', 'linewidth', 3)
        grid on
        ylabel('$x_3$','Interpreter','latex')
        set(gcf,'Color','w');
        set(gca,'FontSize',40);
        subplot(414)
        hold on
        plot( time_vec(1:end), x_4_vec_d(1:end), 'r', 'linewidth', 3)
        grid on
        xlabel('Time (s)','Interpreter','latex', 'FontSize', 40)
        ylabel('$x_4$','Interpreter','latex')
        set(gcf,'Color','w');
        set(gca,'FontSize',40);


        fprintf("\n\n\n End of the program ! \n\n\n");
        
        %%=========================================================
