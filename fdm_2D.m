%% reaction-diffusion FDM solver (2D)
% Rushil Pingali
% Finite Difference Method implementation of reaction-diffusion model for
% Projection Two-Photon Lithography with operator splitting

function[X,Z,DOC,varargout] = fdm_2D(D_xz,coordx,coordz,varargin)

    %% parameters (in S.I. base units)
    
    params_default = {
    'bright_tol', 1e-5;             % time step size during exposure
    'CG_tol', 1e-3;                 % conjugate gradient tolerance
    'D_O2', single(1.2e-12);        % O2 diffusivity
    'dark_period', 0.1;             % post-exposure dark period length
    'dark_tol', 1e-5;               % time step size during dark period
    'dosage_factor', 6.5467e+10;    % dosage factor derived from quantum yield
    'h', 2e-7;                      % mesh step size     
    'kq', single(2300);             % quenching rate constant
    'kt', single(58.58);            % termination rate constant
    'make_plot',false;              % produces DOC plot if enabled
    'O2_0', 6;                      % initial O2 concentration    
    'PETA_0', single(3956);         % initial PETA concentration
    'PI_0', single(1.65);           % initial PI concentration
    'pulse_count', 20;              % number of laser pulses
    'pulse_int', 2e-4;              % time between laser pulses = 1/rep. rate
    'switch_projection_events', []; % times at which to switch projections
    'timeout', inf;                 % stop and return flag = 1 after this time
    'x_bound', [-15e-6, 15e-6];     % x geometry bounds
    'z_bound', [-5e-6, 5e-6];       % z geometry bounds
    };    
    
    parser = inputParser();
    
    for i = 1:size(params_default, 1)
        addParameter(parser, params_default{i, 1}, params_default{i, 2});
    end

    parse(parser, varargin{:});     % user supplied parameters override defaults
    
    % parameters assigned to variables
    bright_tol = parser.Results.bright_tol; 
    CG_tol = parser.Results.CG_tol;
    D_O2 = parser.Results.D_O2;
    dark_period = parser.Results.dark_period;
    dark_tol = parser.Results.dark_tol;
    dosage_factor = parser.Results.dosage_factor;
    h = parser.Results.h;
    kq = parser.Results.kq;
    kt = parser.Results.kt;
    make_plot = parser.Results.make_plot;
    O2_0 = parser.Results.O2_0;
    PETA_0 = parser.Results.PETA_0;
    PI_0 = parser.Results.PI_0;
    pulse_count = parser.Results.pulse_count;
    pulse_int = parser.Results.pulse_int;
    switch_projection_events = parser.Results.switch_projection_events;
    timeout = parser.Results.timeout;
    x_bound = parser.Results.x_bound;
    z_bound = parser.Results.z_bound;
    
    %% geometry
    [X,Z] = meshgrid(x_bound(1):h:x_bound(2),z_bound(1):h:z_bound(2)); % create meshgrid arrays
    Z = flip(Z,1);                      % 1st index in natural ordering is bottom left corner not top left
    [nz,nx] = size(X);                  % number of z and x coordinates
    N = nx*nz;                          % number of mesh nodes
    ones_N = ones(N,1);                 % preallocating empty vectors for performance
    zeros_N5 = zeros(N,5);

    bound_ind = zeros(nz,nx);           
    bound_ind(1,:) = 1;
    bound_ind(end,:) = 1;
    bound_ind(:,1) = 1;
    bound_ind(:,end) = 1;
    bound_ind = fliplr(bound_ind');     % transpose and flip to get natural ordering
    bound_ind = find(bound_ind(:));     % logical array of boundary nodes
    
    %% projections/dosages
    if ~iscell(D_xz)                    % projections supplied as cell array of dosage arrays
        D_xz = {D_xz};
    end
    dos_all = [];
    for i = 1:length(D_xz)
        dos_i = dosage_factor.*interp2(coordx,coordz,D_xz{i},X,Z);   % interpolate at nodes and scale
        dos_i(isnan(dos_i)) = 0;
        dos_i = fliplr(dos_i');         % convert from meshgrid to natural ordering
        dos_i = single(dos_i(:));
        dos_all = [dos_all, dos_i];
    end
    current_projection = 1;             % used for multiple projections
    dos = dos_all(:,current_projection);

    %% time step list
    
    t_range = [0, pulse_count*pulse_int+dark_period];                               % simulation time range
    time_steps = t_range(1):bright_tol:pulse_int*pulse_count;
    time_steps = [time_steps, time_steps(end):dark_tol:t_range(2)];
    pulse_events =1e-5:pulse_int:pulse_count*pulse_int+1e-5;
    time_steps = [time_steps, pulse_events, switch_projection_events];              % adding time steps corresponding to events
    time_steps = unique(time_steps);
    time_steps = sort(time_steps);                                                  % vector of time steps
    pulse_events_ind = ismember(time_steps,pulse_events);                           % logical masks for events
    switch_projection_events_ind = ismember(time_steps,switch_projection_events);
    %%indicators=[switch_projection_events_ind', pulse_events_ind'];
    %%save('inter_data.mat','indicators');
    dt_list = diff(time_steps);

    %% working matrices
    PI = ones_N.*PI_0;
    O2 = ones_N.*O2_0;
    kp = get_kp(0).*ones_N;
    rk1 = zeros_N5;
    rk2 = zeros_N5;
    rk3 = zeros_N5;
    rk4 = zeros_N5;
    conc = zeros_N5;                
    conc(:,2) = ones_N.*PETA_0;
    DOC = zeros(N,1);
    
    % for all N x 5 concentration matrices:
    %     conc(:,1) is [Rstar];
    %     conc(:,2) is [PETA];
    %     conc(:,3) is [Pstar];
    %     conc(:,4) is [Rx];
    %     conc(:,5) is [Px];

    %% A matrix
    a_off = ones_N.*(-1*D_O2/h^2);                                          % off-diagonal elements
    A = spdiags([a_off a_off ones_N a_off a_off], [-nx -1 0 1 nx],N,N);     % initialize sparse matrix A
    A_main_ind = find(speye(N));                                            % indices of main diagonal
    a_main = ones_N;

    %% timestepping
    flag = 0;
    dt = 2;
    DOC_prev = DOC;
    tic
    for i = 3:length(time_steps)
        
        dt_prev = dt;
        dt = single(dt_list(i-1));          % current time step
        if dt~=dt_prev
            A = A.*double((dt/dt_prev));    % scale off-diagonal terms in A
        end
        
        %% check if event
        if switch_projection_events_ind(i)
            current_projection = current_projection + 1;    % move to next projection
            dos = dos_all(:,current_projection);             
        end
        if pulse_events_ind(i)                              % laser pulse event
            conc(:,1) = conc(:,1) + 2.*PI.*dos;
            PI = PI - PI.*dos;
            %%fprintf('dosage updated\n');
        end

        %% update variable vectors

        DOC = 1 - conc(:,2)./PETA_0;
        if norm(DOC_prev-DOC)>1e-3          % update kp when DOC has changed
            kp = get_kp(DOC);
            DOC_prev = DOC;
        end

        %% update ODEs using fourth order Runge-Kutta

        conc_rk = conc;

        kpconc2 = kp .* conc_rk(:,2);
        kpconc1conc2 = kpconc2 .* conc_rk(:,1);
        kqconc1O2 = kq .* conc_rk(:,1).*O2;
        ktO2conc3 = kt .* O2 .* conc_rk(:,3);

        rk_1(:,1) = -kpconc1conc2-kqconc1O2;
        rk_1(:,2) = -kpconc1conc2-kpconc2.*conc_rk(:,3);
        rk_1(:,3) = kpconc1conc2-ktO2conc3;
        rk_1(:,4) = kqconc1O2;
        rk_1(:,5) = ktO2conc3;

        conc_rk = conc + (dt/2).*rk_1;

        kpconc2 = kp .* conc_rk(:,2);
        kpconc1conc2 = kpconc2 .* conc_rk(:,1);
        kqconc1O2 = kq .* conc_rk(:,1).*O2;
        ktO2conc3 = kt .* O2 .* conc_rk(:,3);

        rk_2(:,1) = -kpconc1conc2-kqconc1O2;
        rk_2(:,2) = -kpconc1conc2-kpconc2.*conc_rk(:,3);
        rk_2(:,3) = kpconc1conc2-ktO2conc3;
        rk_2(:,4) = kqconc1O2;
        rk_2(:,5) = ktO2conc3;

        conc_rk = conc + (dt/2).*rk_2;

        rk_3 = kp .* conc_rk(:,2);
        kpconc1conc2 = kpconc2 .* conc_rk(:,1);
        kqconc1O2 = kq .* conc_rk(:,1).*O2;
        ktO2conc3 = kt .* O2 .* conc_rk(:,3);

        rk_3(:,1) = -kpconc1conc2-kqconc1O2;
        rk_3(:,2) = -kpconc1conc2-kpconc2.*conc_rk(:,3);
        rk_3(:,3) = kpconc1conc2-ktO2conc3;
        rk_3(:,4) = kqconc1O2;
        rk_3(:,5) = ktO2conc3;

        conc_rk = conc + dt.*rk_3;

        rk_4 = kp .* conc_rk(:,2);
        kpconc1conc2 = kpconc2 .* conc_rk(:,1);
        kqconc1O2 = kq .* conc_rk(:,1).*O2;
        ktO2conc3 = kt .* O2 .* conc_rk(:,3);

        rk_4(:,1) = -kpconc1conc2-kqconc1O2;
        rk_4(:,2) = -kpconc1conc2-kpconc2.*conc_rk(:,3);
        rk_4(:,3) = kpconc1conc2-ktO2conc3;
        rk_4(:,4) = kqconc1O2;
        rk_4(:,5) = ktO2conc3;

        conc = conc + (rk_1 + 2.*rk_2 + 2.*rk_3 + rk_4).*dt./6;

        %% update PDE using backward Euler
        
        a_main_new = 1 + dt.*(kq.*conc(:,1)+kt.*conc(:,3)+(4.*D_O2)./(h.^2)); % updated main diagonal elements
        a_main_new(bound_ind) = 1;
        if norm(a_main_new-a_main)>1e-6      % if main diagonal elements have changed, update them in A
           A(A_main_ind) = a_main_new;       % ignore warning: linear indexing gives faster performance than spdiags
        end
        a_main = a_main_new;

        %% apply conjugate gradient method to update O2
        x = O2;
        b = O2;
        r_current = b - A*x;
        p_current = r_current;
        r_norm = norm(r_current);
        while r_norm > CG_tol
            if toc > timeout            % optionally can set a timeout
                flag = 1;
                break
            end      
            Ap_current = A*p_current;
            a_current = (r_norm^2)/(p_current'*Ap_current);
            x = x + a_current*p_current;
            r_next = r_current - a_current*Ap_current;
            r_next_norm = norm(r_next);
            b_current = (r_next_norm^2)/(r_norm^2);
            p_current = r_next + b_current*p_current;
            r_current = r_next;
            r_norm = r_next_norm;        
        end
        if flag == 1
            break
        end
        O2 = x;
    O2(bound_ind) = O2_0;               % reapply boundary conditions
    end

    %% output DOC
    PETA = reshape(conc(:,2),[nx,nz]);
    PETA = fliplr(PETA)';               % flip and transpose to go from natural ordering to meshgrid
    DOC = 1 - PETA./PETA_0;

    %% plot
    if make_plot
        figure
        surf(imresize(X,5),imresize(Z,5),imresize(DOC,5),'EdgeColor', 'None', 'facecolor', 'interp');
        view(2); 
        colormap jet; 
        c = colorbar;
        c.Label.String = 'Degree of Conversion';
        xlabel('X (m)');
        ylabel('Z (m)');
    end
    varargout{1} = flag;
end

%% returns kp(DOC) in m3/(s mol)
% based on literature data for k_p of PETA (see text and references of
% https://doi.org/10.1115/1.4051830)
function[kp] = get_kp(DOC)
   kp_0 = single(43);  % m3/(s mol)
   kp = (-5.95).*DOC+0.9914;
   mask_2 = DOC>=0.16&DOC<0.4;
   kp(mask_2) = (20.1).*exp(-38.8.*DOC(mask_2));
   mask_3 = DOC>=0.4;
   kp(mask_3) = (20.1).*exp(-38.8.*0.4);
   kp = kp_0.*kp;
end