function [EF, varargout]= sers_enhancement_factor(waveNumber, ramanSpectrum, SERSspectrum, varargin)
% Estimates surface-enhanced Raman spectroscopy (SERS) enhancement factor
%
% SYNTAX
% [EF, paramStruct] = sers_enhancement_factor(waveNumber, ramanSpectrum, ...
%                       SERSspectrum, 'Name', Value)
%
% INPUTS
%      REQUIRED
%    waveNumber Raman shift vector (cm-1)
%  ramanSpectrum Normal Raman spectrum (a.u.)
%  SERSspectrum SERS Raman spectrum (a.u.)
% The above required input parameters, can be followed by Name/Value pairs
% to specify additional variables of the computations.
%    [OPTIONAL]
%         band: Raman band of interest (cm^-1)
%       lambda: Wavelength of excitation (nm)
%    molWeight: Molecular weight of the molecule undeer test (g/mol)
%           NA: Numerical aperture of the objective
%          rho: Density of molecule under test (g/cm^3)
%     surfArea: Surface area of molecule under test (nm^2)
%     showPlot: Describes whether to display a plot of the spectra & E.F.
% For example,
% sers_enhancement_factor(waveNumber,ramanSpectrum, SERSspectrum,'band',1262)
% will estimate E.F. at the band closest to 1262 cm^-1.
%
% OUTPUTS
%      REQUIRED
%            EF SERS enhancement factor
%    [OPTIONAL]
%   paramStruct Structure containing the following parameters:
%         band: Raman band of interest (cm^-1)
%      Inormal: Intensity of normal Raman spectrum at given band (a.u.)
%        ISERS: Intensity of SERS spectrum at given band (a.u.)
%       lambda: Wavelength of excitation (nm)
%    molWeight: Molecular weight of the molecule undeer test (g/mol)
%           NA: Numerical aperture of the objective
%      Nnormal: Number of excited molecules in normal Raman
%        NSERS: Number of excited molecules in SERS
% ramanSpectrum:Normal Raman spectrum (a.u.)
%          rho: Density of molecule under test (g/cm^3)
% SERSspectrum: SERS Raman spectrum (a.u.)
%     surfArea: Surface area of molecule under test (nm^2)
%     showPlot: Describes whether to display a plot of the spectra & E.F.
%   waveNumber: Raman shift vector (cm-1)
%            r: Radius of the laser spot (m)
%            h: Depth of focus (m)
%       volExc: Excitation volume (m^3)
%
% Example:
% load('SERS_data_example.mat');
% [EF, paramStruct] = sers_enhancement_factor(waveNumber, ramanSpectrum, SERSspectrum);
% See sers_enhancement_factor_script.m for an example of usage
%
% Reference:
% Childs, A., Vinogradova, E., Ruiz-Zepeda, F., Velazquez-Salazar, J. J., &
% Jose-Yacaman, M. (2016). Biocompatible gold/silver nanostars for
% surface-enhanced Raman scattering. Journal of Raman Spectroscopy, 47(6),
% 651–655. https://doi.org/10.1002/jrs.4888
%
%_______________________________________________________________________________
% Copyright (C) 2016 Edgar Guevara, PhD
% CONACYT-Universidad Autónoma de San Luis Potosí
% Coordinación para la Innovación y Aplicación de la Ciencia y la Tecnología
%_______________________________________________________________________________

%Make sure inputs are columns
if ~iscolumn(waveNumber)
    waveNumber = waveNumber(:);
end
if ~iscolumn(ramanSpectrum)
    ramanSpectrum = ramanSpectrum(:);
end
if ~iscolumn(SERSspectrum)
    SERSspectrum = SERSspectrum(:);
end

% Number of outputs must be >=minargs and <=maxargs.
minargs=1;  maxargs=2;
nargoutchk(minargs, maxargs);

% Check the validity of required and optional function inputs
p = inputParser;

% Default input values
[~, defaultBandIdx]= max(SERSspectrum);
defaultBand = waveNumber(defaultBandIdx);   % Default band (cm^-1)
defaultLambda = 785;        % Default excitation wavelength (nm)
defaultNA = 0.25;           % Default numerical aperture of the objective
defaultInormal = 1;         % Default intensity of the normal Raman spectrum
defaultISERS = 1;           % Default intensity of the SERS spectrum
defaultNnormal = 1;         % Default No. of molecules of the normal Raman spectrum
defaultNSERS = 1;           % Default No. of molecules of the SERS spectrum
defaultRho = 1.26;          % Default density of the molecule of interest(g/cm3)
defaultMolWeight = 479.02;  % Default weight of interest(g/mol)
defaultSurfArea = 4;        % Default surface area of the molecule (nm^2)
defaultShowPlot = true;     % Default option to show plot of the E.F. estimation
% ------------------------ Add required parameters ------------------------
addRequired(p,'waveNumber',@isnumeric);
addRequired(p,'ramanSpectrum',@isnumeric);
addRequired(p,'SERSspectrum',@isnumeric);
% ------------------------ Add optional parameters ------------------------
addOptional(p,'band',defaultBand, @isnumeric);
addOptional(p,'lambda',defaultLambda, @isnumeric);
addOptional(p,'NA',defaultNA, @isnumeric);
addOptional(p,'rho',defaultRho, @isnumeric);
addOptional(p,'molWeight',defaultMolWeight, @isnumeric);
addOptional(p,'surfArea',defaultSurfArea, @isnumeric);
addOptional(p,'showPlot',defaultShowPlot, @islogical);
% --------------------- Add parameters to be computed ---------------------
addParameter(p,'Inormal',defaultInormal);
addParameter(p,'ISERS',defaultISERS);
addParameter(p,'Nnormal',defaultNnormal);
addParameter(p,'NSERS',defaultNSERS);

% Parse input values
parse(p, waveNumber, ramanSpectrum, SERSspectrum, varargin{:});

% Compute E.F.
[EF, paramStruct] = compute_EF(p);

% Optional output argument
varargout{1} = paramStruct;
end

function [EF, paramStruct] = compute_EF(p)
% Computes enhancement factor estimation
% Copy Results to pRes structur
paramStruct = p.Results;
% Find normal intensity at the peak closest to the specified band
[paramStruct.band, closestIdx] = findClosest(paramStruct);
paramStruct.Inormal = paramStruct.ramanSpectrum(closestIdx);
% Find SERS intensity at the specified band
paramStruct.ISERS = paramStruct.SERSspectrum(closestIdx);
% lambda in m
lambda = paramStruct.lambda/1e9;
% Height of the cylinder probed by Raman laser (m)
paramStruct.h = (2*lambda) / paramStruct.NA^2;
% First zero of the Bessel function of the first kind, order 1 (divide by 2*pi)
myZero = fzero('besselj(1,x)', 3.8) / (2*pi);
% Radius of the laser spot (m)
paramStruct.r = myZero*lambda / paramStruct.NA;
% volume of the focused laser spot (optical excitation volume m^3)
paramStruct.volExc = pi*paramStruct.r^2*paramStruct.h;
% Avogadro constant (mol^-1)
AN = 6.022140857e23;
% Convert density from g/cm^3 to g/m^3
paramStruct.Nnormal = paramStruct.volExc * (paramStruct.rho / 1e-6) * AN...
    / paramStruct.molWeight;
% Approximate computation of NSERS
paramStruct.NSERS = pi*paramStruct.r^2 /...
    (paramStruct.surfArea / 1e18);  % Surface area in m^2
% Compute E.F.
EF = (paramStruct.ISERS * paramStruct.Nnormal) / ...
    (paramStruct.Inormal * paramStruct.NSERS);
% Plot spectra and display E.F. estimation
if paramStruct.showPlot
    figure; plot(paramStruct.waveNumber, paramStruct.ramanSpectrum, 'k-',...
        paramStruct.waveNumber, paramStruct.SERSspectrum, 'r-', ...
        paramStruct.waveNumber(closestIdx), paramStruct.SERSspectrum(closestIdx), 'bo')
    set(gcf,'color','w')
    text(1.025*paramStruct.band, paramStruct.ISERS, ...
        sprintf('k = %0.1f cm^{-1}\nE.F. = %0.5G',paramStruct.band, EF))
    legend({'Normal' 'SERS'}, 'Location', 'NorthWest')
    xlabel('Raman shift [cm^-1]')
    ylabel('Raman intensity [a.u.]')
end
end

function [closestValue, closestIdx] = findClosest(paramStruct)
% Finds closest Raman peak to the proposed band
% Check if signal processing toolbox is installed
if license('test', 'Signal_Toolbox')
    % Finds index of the element in paramStruct.waveNumber closest to paramStruct.band
    [~, closestWaveNumberIdx] = min(abs(paramStruct.waveNumber - paramStruct.band));
    prominenceFactor = 10/100;      % Adjust as needed
    [ramanPeaks, ramanIdx] = findpeaks(paramStruct.SERSspectrum, ...
        'MinPeakProminence', prominenceFactor*max(paramStruct.SERSspectrum));
    % Finds index of the element in ramanIdx closest to closestWaveNumberIdx
    [~, closestPeakIdx] = min(abs(ramanIdx - closestWaveNumberIdx));
    closestIdx = ramanIdx(closestPeakIdx);
    closestValue = paramStruct.waveNumber(closestIdx);
    % commented out - plot of SERS peaks
    % figure; plot(paramStruct.waveNumber, paramStruct.SERSspectrum, 'k-');
    % hold on; plot(paramStruct.waveNumber(ramanIdx), ramanPeaks, 'ro')
else
    % Simply find the wave number closest to proposed band
    [~, closestWaveNumberIdx] = min(abs(paramStruct.waveNumber - paramStruct.band));
    closestValue = paramStruct.waveNumber(closestWaveNumberIdx);
end
end

% EOF

