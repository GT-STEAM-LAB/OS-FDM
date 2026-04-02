%% Simulation of temporal focusing system
% Sourabh Saha and Rushil Pingali

% Evaluating optical dosage and intensity in the focal volume
% f1 and f2 may be unequal for the two lenses
% Paraxial approximation
% All units in micrometers (um) unless specified

function[I_fz,D_fz,Xf,Yf,zf] = optical_dosage(H_dmd,system,varargin)

    load(system);   % system parameters

    %% simulation parameters
    params_default = {
    'N', 256;                % Number of points in the spatial field (power of 2)
    'N_omega', 256;         % Number of points in spectrum sample (note accuracy improves with higher values, best speeds if power of 2)
    'dLx', 7.56;             % Distance sampling step size (DMD pitch) (um)
    'make_plot', false;      % Produces dosage plot
    'calibration_power', calibration.power;     % Default system power
    'Zmin', -5;              % Minimum z value for dosage simulation (um)
    'Zmax', 5;              % Maximum z value for dosage simulation (um)
    'Zstep', 0.1;           % z step size for dosage simulation (um)
    };

    parser = inputParser();
    for i = 1:size(params_default, 1)
        addParameter(parser, params_default{i, 1}, params_default{i, 2});
    end
    parse(parser, varargin{:});

    N = parser.Results.N;
    N_omega = parser.Results.N_omega;
    dLx = parser.Results.dLx;
    make_plot = parser.Results.make_plot;
    calibration.power = parser.Results.calibration_power;
    Zmin=parser.Results.Zmin;
    Zmax=parser.Results.Zmax;
    Zstep=parser.Results.Zstep;
    
    %% Evaluating system parameters
    f2=f_tube/M_obj;             %Focal length of objective lens 
    Mag=f2/f1;                     %Magnification of optical system
    D_obj=2*NA*f2/n_lens;          %Pupil diameter of objective lens

    %Bandwidth parameters
    omega_o=c_light/lambda_o;    %Center frequency of light (Hz)
    omega_band=c_light*wave_band/lambda_o^2;  %Spectral bandwidth:FWHM (Hz)
    omega=omega_band/(2*sqrt(log(2)));  %Half-spectral bandwidth: 1/e2 (Hz)

    %% Simulation parameters

    %Source space definitions
    Lx=N*dLx;                    %Size of field in object plane
    x=(-N/2:N/2-1)*dLx;         %x coordinates (um)
    [Xo,Yo]=meshgrid(x);        %Mesh of x,y coordinates in object space (um)

    %Spectral sample definitions
    dOmega=c_light*0.5e-3/lambda_o^2;  %Frequency step size, corresponds to 0.5 nm bandwidth (Hz)
    f_wave=omega_o + (-N_omega/2:N_omega/2-1)*dOmega; %Spectral components
    Uw=exp(-(f_wave-omega_o).^2/omega^2);           %Spectral electric field
    lambda=c_light./f_wave;                 %Wavelength components

    %% Simulation of optical system
    %First evaluating field immediately after DMD for each wavelength (one loop over wavelength)
    %Next for each z plane evaluating sum of field for each wavelength...
    %(nested loops - z on outer loop, wavelength on inner loop)
    %Not storing field variable for each z plane; storing only intensity of
    %each z plane

    %Initialization of field 
    U=zeros(N,N,N_omega);                     %Running field variable (value updates to different planes)
    U_dmd=zeros(N,N,N_omega);                 %Field variable after DMD (holds same value across all Z steps)

    %Loop over each wavelength
    for nw=1:N_omega

        %Field before DMD
        U(:,:,nw)=Uw(nw);
        %Field after DMD (DMD pattern and linear phase from grating effect)
        U(:,:,nw)=U(:,:,nw).*H_dmd.*(exp(2i*pi.*(Xo+Yo)*m_order*(lambda(nw)-lambda_o)/(lambda_o*dLx)));

    end

    %%Field in the focal volume of objective

    U_dmd=U;                    %This value is held across all z steps

    %Propagation around the focal volume & summation over wavelengths
    zf=(Zmin:Zstep:Zmax)';                 %Z coordinates around focal volume

    I_fz=zeros(N,N,size(zf,1));     %Intensity around focal volume
    D_fz=zeros(N,N,size(zf,1));     %Dosage around focal volume
    I_fz_t=zeros(N,N,N_omega);      %Time varying intensity at each z plane
    Xf=Xo*Mag;Yf=Yo*Mag;            %X-Y coordinates in the focal volume   
    U_t=zeros(N,N,N_omega);         %Time varying electric field at each z, sum of all wavelengths
    temporal=zeros(N_omega);        %Field at any given space point at various time/frequency steps
    time_pulse=(-N_omega/2:N_omega/2-1)/(N_omega*dOmega);%time coordinates corresponding to sampled spectrum 

    lambda_3D = reshape(lambda, 1, 1, []);
    %At back focal plane of objective
    [U_back, Fx, Fy] = getOpticalFFT_nD(U_dmd, dLx);
    U_back = U_back.*(-1i./(lambda_3D.*f1));

    %Pupil function of objective
    Xbf = Fx.*lambda_3D.*f1; %Spatial coordinates in back-focal plane
    Ybf = Fy.*lambda_3D.*f1; %Spatial coordinates in back-focal plane
    Rbf_sq_3D=(Xbf.^2+Ybf.^2); %radius square in back focal plane
    P_obj_3D=((Rbf_sq_3D<(0.5*D_obj)^2)*1 + 0.5*(Rbf_sq_3D==(0.5*D_obj)^2)); 
    P_obj_3D_ind = find(P_obj_3D);              % linear indices of non-zeros in P
    [~,~,ind3] = ind2sub([N,N,nw],P_obj_3D_ind);   % row/col/page indices
    U_back_filtered = U_back(P_obj_3D_ind);
    factor_1 = (-2i.*pi.*Rbf_sq_3D(P_obj_3D_ind).*n_lens./(f2.^2.*lambda(ind3)'));
    dLb_all = (lambda'.*f1) .* reshape(abs(Fx(1,2,:) - Fx(1,1,:)),[],1); %Spatial grid spacing in back-focal plane 
    factor_3 = -1i.*n_lens./(lambda_3D.*f2);
    
    for nz=1:size(zf,1)         % for loop is necessary for memory management
        z=zf(nz);
        U = getOpticalFFT_nD(U_back_filtered.*exp(factor_1.*z), dLb_all,[N,N,nw],P_obj_3D_ind);
        factor_2 = exp(2i.*pi.*(2.*f2+z).*n_lens./lambda_3D);
        U = factor_2.*U.*factor_3;

        %Evaluating intensity at each z by summing up all wavelengths
        %Net electric field obtained by temporal Fourier transform at each space point in sample
        U_t = (dOmega)*fftshift(fft(fftshift(U,3),[],3),3);

        %Intensity at z from field at z
        I_fz_t=abs((U_t().^2));   %Time varying intensity at z
        I_fz(:,:,nz)=mean(I_fz_t,3);
        D_fz(:,:,nz)=sum(I_fz_t.^2, 3);
    end

    %Scaling intesity to TW/cm2 units
    UtoI=1e-16*sqrt(2/pi)*calibration.eta*calibration.power/(omega*calibration.Rep*calibration.Area);
    I_fz=I_fz*UtoI;
    D_fz=D_fz*UtoI^2*(time_pulse(2)-time_pulse(1));

    %% Diagram

    if make_plot
        D_xz=squeeze(D_fz(N/2+1,:,:));
        [coordx, coordz]=meshgrid(Xf(1,:,1),zf); coordx=coordx'; coordz=coordz';

        %Plotting the X-Z profile of dosage
        figure;
        surf(coordx,coordz,D_xz,'EdgeColor', 'None', 'facecolor', 'interp');
        view(2); colormap jet; 
        c = colorbar;
        xlabel('X (\mum)');
        ylabel('Z (\mum)');
        c.Label.String = 'Dosage TW^{2}s/cm^{4}';
        c.Location = 'southoutside';
        c.Limits = [0,12E-13];
        c.TickDirection = 'out';
        set(gca,'TickDir','out');
    end

    %% convert coordinates to meters before returning
    Xf = Xf.*1e-6;
    Yf = Yf.*1e-6;
    zf = zf.*1e-6;
end

%% getOpticalFFT: performs multiple 2D Optical FFTs
function [u2, Fx, Fy] = getOpticalFFT_nD(u1,dL,varargin)

    % Input parameters:
    %   u1: fields at source
    %   dL: spatial sampling periods (um)
    %   varargin{1}: dims of 3D array
    %   varargin{2}: non-zero indices of 3D array
    % Output parameters:
    %   u2: fields at observation plane
    %  Fx,Fy: spatial frequency coordinates in the observation planes 
    
    if nargin > 2       % data has been supplied in sparse col format
                        % # of rows = # of non-zero matrices in 3D array
                        % # of cols = # of 3D arrays
        dims = varargin{1};
        n_elements = dims(1)*dims(2)*dims(3);
        N = dims(1);    %Sample size
        sz = size(u1);
        n_mat = sz(2);
        u0 = (zeros(dims(1),dims(2),dims(3)*n_mat));
        for i = 1:n_mat
            inds = varargin{2}+(i-1)*n_elements;
            u0(inds) = u1(:,i);
        end
        u1 = u0;
    else                % data has been supplied in multidim array format, convert to 3D
        sz = size(u1);  
        u1 = reshape(u1, sz(1), sz(2), []);
        N=sz(1);           %Sample size
    end
    
    if nargout > 1
        % Setting up frequency coordinates
        dfs=1./(N.*dL);           %frequency step
        fy=(-N/2:N/2-1)'.*dfs';    %x freq coord
        Fy = repmat(fy, [1 1 N]);
        Fy = permute(Fy, [1 3 2]);
        Fx = permute(Fy, [2 1 3]);
    end
    
    %Evaluating 2-D FFT for each plane
    u2 = reshape((dL.^2),1,1,[]) .* fftshift(fftshift(fft(fft(fftshift(fftshift(u1,1),2),[],1),[],2),1),2);            %FFT of source field

end