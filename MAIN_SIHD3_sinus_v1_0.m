        
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
        warning off
        clc
        
        % SIMULATION SETTINGS        
        TMax = 10;      % Maximum simulation time
        Deltah = 1e-3;  % Time-step
        
        % SIHD-3 parameters
        m_alpha_1 = 0.95;
        m_alpha_2 = 0.95;
        m_alpha_3 = 0.95;
        m_alpha_4 = 0.95;
        
        m_lambda_1 = 1e3;
        m_lambda_2 = 1e6;
        m_lambda_3 = 1e9;
        m_lambda_4 = 1e9;
        
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
        
        % Noise vector creation
        NoiseAmp = 0; % Noise amplitude (0 or 1e-6 in the paper) 
        eta_noise = NoiseAmp * rand( 1, length( 0:Deltah:TMax ) );
        
        % Sine definition
        freq = 1 / (2 * pi );
        ph = 0;
        amp = 1;
        
        % I.C. of the sine function and its derivatives
        x1_d0 = amp * sin (2 * pi * freq  * (0) - ph );
        
        x2_d0 = (2 * pi * freq) * amp * cos ( 2 * pi * freq * (0) - ph);
        
        x3_d0 = -( 2 * pi * freq )^2 * amp * sin ( 2 * pi * freq  * (0) - ph);
        
        x4_d0 = -( 2 * pi * freq )^3 * amp * cos ( 2 * pi * freq  * (0) - ph);
        
        % Set the I.C. to initialize the SIHD-3 (w.r.t. the I.C. of the
        % sine function and its derivatives)
        mz1_d = x1_d0;
        mz2_d = x2_d0;
        mz3_d = x3_d0;
        mz4_d = x4_d0;
        
        fprintf('\n lambda_1 = %f \n lambda_2 = %f \n lambda_3 = %f \n lambda_4 = %f \n', m_lambda_1, m_lambda_2, m_lambda_3, m_lambda_4 );
        
        fprintf('\n theta_1 = %f \n theta_2 = %f \n theta_3 = %f \n theta_4 = %f \n', m_theta_1, m_theta_2, m_theta_3, m_theta_4 );
        
        % SIHD Proj. init.
        mE_1_m = 0;
        mE_2_m = 0;
        mE_2_mm = 0;
        mE_3_m = 0;
        mE_3_mm = 0;
        mE_3_mmm = 0;
        
        m_correction_factor = 1; % Enable the correction terms (enable = 1)
        uu = 0; % counter of the loop
        time = 0; % time
        while time <= TMax
        
            uu = uu + 1; % counter increment
        
            % Time update
            time = Deltah * uu;
            time_vec(uu) = time;
        
            % Sine update w.r.t. the time
            x1_d = amp * sin (2 * pi * freq  * time - ph );
            x2_d = (2 * pi * freq) * amp * cos ( 2 * pi * freq * time - ph);
            x3_d = -( 2 * pi * freq )^2 * amp * sin ( 2 * pi * freq  * time - ph);
            x4_d = -( 2 * pi * freq )^3 * amp * cos ( 2 * pi * freq  * time - ph);
        
            % Add noise to the output y_meas = x_1d
            y_meas = x1_d + eta_noise(uu);
        
            % ======== SIHD-3 section ========
        
            % Compute the abs. errors
            me_1_vec(uu) = abs(y_meas - mz1_d);
            me_2_vec(uu) = abs(x2_d - mz2_d);
            me_3_vec(uu) = abs(x3_d - mz3_d);
            me_4_vec(uu) = abs(x4_d - mz4_d);
        
            % Compute the error e_1
            me_1 = y_meas - mz1_d;
        
            % Compute the projectors
            [ mE_1, mProj_1_ ] = Proj_function( m_alpha_1, m_lambda_p1, 1, me_1, m_MU_, Deltah);
            [ mE_2, mProj_2_ ] = Proj_function( m_alpha_2, m_lambda_p2, 2, me_1, m_MU_, Deltah);
            [ mE_3, mProj_3_ ] = Proj_function( m_alpha_3, m_lambda_p3, 3, me_1, m_MU_, Deltah);
            [ mE_4, mProj_4_ ] = Proj_function( m_alpha_4, m_lambda_p4, 4, me_1, m_MU_, Deltah);
        
            % Update the differentiation
            mz4_p = mz4_d + (Deltah * mE_1_m * mE_2_m * mE_3 * m_lambda_4 * ( ( m_MU_  )^4 ) * ( abs(  me_1 ) )^(4 * m_alpha_4 - 3 ) * mProj_4_ );
        
            mz3_p = mz3_d + (Deltah * mE_1_m * mE_2_m * ( mz4_p + m_lambda_3 * ( ( m_MU_  )^3 ) * ( abs(  me_1 ) )^(3 * m_alpha_3 - 2 ) * mProj_3_ ) );
        
            mz2_p = mz2_d + Deltah * mE_1_m * ( (( mz3_p  )) - ((m_correction_factor * 1/2  * Deltah * mz4_p )) + m_lambda_2 * mE_1_m * ( ( m_MU_  )^2 ) * ( abs(  me_1 ) )^(2 * m_alpha_2 - 1 ) * mProj_2_ );
        
            mz1_p = mz1_d + Deltah * ( mz2_p - ((m_correction_factor * 1/2  * Deltah * (( mz3_p ))  )) + ((m_correction_factor * 1/6 * Deltah^2 * mz4_p )) + m_lambda_1 * ( ( m_MU_  ) ) * ( abs( me_1 ) )^m_alpha_1 * mProj_1_  );
        
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
        
            % Vectorization
            mz1_d_vec(uu) = mz1_d;
            mz2_d_vec(uu) = mz2_d;
            mz3_d_vec(uu) = mz3_d;
            mz4_d_vec(uu) = mz4_d;
        
            x_1_vec_t(uu) = x1_d;
            x_2_vec_t(uu) = x2_d;
            x_3_vec_t(uu) = x3_d;
            x_4_vec_t(uu) = x4_d;
        
            x_1_vec_d(uu) = mz1_d;
            x_2_vec_d(uu) = mz2_d;
            x_3_vec_d(uu) = mz3_d;
            x_4_vec_d(uu) = mz4_d;
        
            % ======== End of the SIHD-3 section ========
        
        end
        
        % Display / Plot results
        le_me_vec = length( me_1_vec );
        
        fprintf('\n mean_e1 = %e', mean( me_1_vec(floor(le_me_vec/2):end) ) );
        
        fprintf('\n mean_e2 = %e', mean( me_2_vec(floor(le_me_vec/2):end) ) );
        
        fprintf('\n mean_e3 = %e', mean( me_3_vec(floor(le_me_vec/2):end) ) );
        
        fprintf('\n mean_e4 = %e', mean( me_4_vec(floor(le_me_vec/2):end) ) );
        
        
        figure('name','SIHD3 error')
        subplot(411)
        semilogy( time_vec(1:end), me_1_vec(1:end), 'b', 'linewidth', 3)
        hold on
        grid on
        xlim([0, 10])
        ylabel('$|e_1| \quad $','Interpreter','latex', 'FontSize', 40)
        grid on
        set(gca,'FontSize',40)
        subplot(412)
        semilogy( time_vec(1:end), me_2_vec(1:end), 'b', 'linewidth', 3)
        grid on
        xlim([0, 10])
        ylabel('$|e_2| \quad $','Interpreter','latex', 'FontSize', 40)
        set(gcf,'Color','w');
        set(gca,'FontSize',40);
        subplot(413)
        semilogy( time_vec(1:end), me_3_vec(1:end), 'b', 'linewidth', 3)
        xlim([0, 10])
        grid on
        ylabel('$|e_3| \quad $','Interpreter','latex', 'FontSize', 40)
        set(gcf,'Color','w');
        set(gca,'FontSize',40);
        grid on
        subplot(414)
        semilogy( time_vec(1:end), me_4_vec(1:end), 'b', 'linewidth', 3)
        grid on
        xlim([0, 10])
        xlabel('Time (s)','Interpreter','latex')
        ylabel('$|e_4| \quad $','Interpreter','latex', 'FontSize', 40)
        set(gcf,'Color','w');
        set(gca,'FontSize',40);
        xlabel('Time (s)','Interpreter','latex', 'FontSize', 40)
        set(gcf,'Color','w');
        set(gca,'FontSize',40);
        
        
        figure('name','SIHD3 states')
        subplot(411)
        hold on
        plot( time_vec(1:end), x_1_vec_d(1:end), 'r', 'linewidth', 3);
        plot( time_vec(1:end), x_1_vec_t(1:end), '--b', 'linewidth', 3)
        legend('estimated','exact', 'FontSize', 40, 'Interpreter','latex')
        xlim([0, 10])
        ylabel('$x_1 \quad $','Interpreter','latex', 'FontSize', 40)
        grid on
        set(gca,'FontSize',40)
        subplot(412)
        hold on
        plot( time_vec(1:end), x_2_vec_d(1:end), 'r', 'linewidth', 3)
        plot( time_vec(1:end), x_2_vec_t(1:end), '--b', 'linewidth', 3)
        grid on
        legend('estimated','exact', 'FontSize', 40, 'Interpreter','latex')
        xlim([0, 10])
        ylabel('$x_2 \quad $','Interpreter','latex', 'FontSize', 40)
        set(gcf,'Color','w');
        set(gca,'FontSize',40);
        subplot(413)
        hold on
        plot( time_vec(1:end), x_3_vec_d(1:end), 'r', 'linewidth', 3)
        plot( time_vec(1:end), x_3_vec_t(1:end), '--b', 'linewidth', 3)
        legend('estimated','exact', 'FontSize', 40, 'Interpreter','latex')
        grid on
        xlim([0, 10])
        ylabel('$x_3 \quad $','Interpreter','latex', 'FontSize', 40)
        set(gcf,'Color','w');
        set(gca,'FontSize',40);
        grid on
        subplot(414)
        hold on
        plot( time_vec(1:end), x_4_vec_d(1:end), 'r', 'linewidth', 3)
        plot( time_vec(1:end), x_4_vec_t(1:end), '--b', 'linewidth', 3)
        legend('estimated','exact', 'FontSize', 40, 'Interpreter','latex')
        grid on
        xlim([0, 10])
        ylabel('$x_4 \quad $','Interpreter','latex', 'FontSize', 40)
        set(gcf,'Color','w');
        xlabel('Time (s)','Interpreter','latex', 'FontSize', 40)
        set(gcf,'Color','w');
        set(gca,'FontSize',40);
        
        
        fprintf("\n\n\n End of the program ! \n\n\n");
        
        %%=========================================================

