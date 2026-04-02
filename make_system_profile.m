%% save projection system parameters

% set variables to correct values for your system and run this script
% all units in micrometers (um) unless specified

system_name = 'STEAM_Lab_GT';

%% System parameters

M_obj=60;                %Designed magnification of the objective lens (×)
f_tube=180e3;              %Designed value of tube lens focal length for the chosen objective (um)
f1=200e3;                  %Actual focal length of the tube lens (um)
NA=1.25;                    %Numerical aperture of objective
n_lens=1.52;               %Refractive index of immersion medium of objective lens

%Light parameters

lambda_o=804e-3;               %Center wavelength of light (um)
c_light=3e8*1e6;                  %Speed of light in vacuum (um/s)
wave_band=41e-3;                %Wavelength spectral bandwidth of light, FWHM(um)
 
%Grating parameters
m_order=4;                      %Diffraction grating order 

%Power calibration parameters

calibration.power=139.3e-9;       %(DEFAULT) Measured average power of beam at the entrance to the objective lens per pixel (W/pixel)
calibration.Rep=5e3;            %Repetition rate of laser (Hz)
calibration.Area=5.715e-11;     %Area of single pixel on the DMD (m2)
calibration.eta=0.727;            %Optical efficiency through objective lens

save([system_name '.mat']);