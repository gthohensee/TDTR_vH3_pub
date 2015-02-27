%process_yymmdd_template - Template for processing TDTR data
%This script prepares raw TDTR data for thermal analysis. Steps include:
%finding t0, setting the phase, picking roundtrip echoes or Brillouin
%oscillations, adding an offset to V(out) to correct for spurious signal 
%as determined by blocked-pump or blocked-probe configurations, examining 
%the change in V(out) across time delay to judge the extent of pump drift 
%vs. delay time, and other data verification calculations.
%
% Basic use of this template:
% First, do a Ctrl+R search-and-replace for the current date in yymmdd
% representation. Next, "Save As" this file with a new name, e.g.
% process_140322_[dataseries_ID], such as process_140322_AlDAC.
% Then, sequentially run through the code, editing it as indicated by
% the comments and common sense. This file can serve as your notebook for
% how you handled this particular data series!
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: analyze_yymmdd_template.m, TDTR_MAIN_vH.m

% Author: Gregory T. Hohensee
% University of Illinois Urbana-Champaign
% email: hohense2@illinois.edu
% Website: n/a
% April 2014; Last revision: 3-Apr-2014
%                            12-May-2014 - save results to savefolder, not
%                                          datafolder.
%                            14-July-2014 - vH2. No changes.
%                            17-Feb-2015 - vH3. No changes.
%------------- BEGIN CODE --------------

% Your name and any identifying (time)stamp as desired. %
% [Short descriptor of the sample(s) and conditions of measurement] %

%% Define the directory and filenames
directory = '/Users/gregoryhohensee/Documents/MATLAB/BeamOffset';
yymmdd = '140803'; % e.g., 140227 = February 27th, 2014.
datafolder = strcat('/',yymmdd,'/');      % folder in directory where raw data is kept
savefolder = strcat('/',yymmdd,'_edit/'); % folder in directory where analyzed data is saved
datadir = strcat(directory, datafolder);
savedir = strcat(directory, savefolder);

tagline = 'Al'; % common prefix for your series of TDTR data files.

% get filenames in alphanumerical order
fnames = dir(strcat(directory,datafolder,tagline,'*'));
fnames.name
nraw = length(fnames);
%% initialize variables
data = cell(1,nraw);
echoes = zeros(nraw,5); % assuming at most 5 echoes picked
delphases = zeros(nraw,1);
phase_shift = zeros(nraw,1);
t0 = zeros(nraw,1);
fitparam = zeros(nraw,3);
%% Process all files through SetPhase_vH.m
for k=1:nraw
  fname = fnames(k).name;
  data = dlmread(strcat(directory,datafolder,fname));
  
  % If you have measured a mean offset in V(out) for a blocked-pump
  % configuration, uncomment the next four lines to adjust for that offset.
  %Vout_offset = 0;
  %Voutmean = mean(data(:,4));
  %data(:,4) = data(:,4) - Vout_offset;
  %data(:,5) = -data(:,3)./data(:,4);
  
  % SetPhase_vH picks t0 and sets the phase. PickEchoes reads acoustics.
  t_window = [-10 10]; tmax = 80;
  [t0(k)] = set_t0(data, t_window); % Set t0
  [delphases(k),phase_shift(k),fitparam(k,:)] = SetPhase_vH3(data,t0(k),fname,savedir,t_window);
  echotime = PickEchoes(data,tmax,t0(k),t_window);
  echoes(k,1:length(echotime)) = echotime;
  clear echotime;
  
  % This conversion is (t*6.42/2)+3 nm
  % Al/SiO2 - 28.9 ps - 95.7 nm
  % Al/Si - 20.9 ps - 70 nm

end
%% Save delphases and echoes
shifts = horzcat(t0, phase_shift);
dlmwrite(strcat(directory, savefolder, yymmdd, '_delphases_', tagline, '.txt'),    {delphases});
dlmwrite(strcat(directory, savefolder, yymmdd, '_shifts_', tagline, '.txt'),       {shifts});
dlmwrite(strcat(directory, savefolder, yymmdd, '_echotimes_', tagline, '.txt'),    {echoes});
dlmwrite(strcat(directory, savefolder, yymmdd, '_Voutfitparams_', tagline, '.txt'),{fitparam});
%% Check the percent change in Vout over 4 ns; retain the sign of Vout.
for k = 1:nraw
    fname = fnames(k).name;
    data = dlmread(strcat(directory,savefolder,fname,'_shifted.txt'));
    L = length(data(:,4));
    fVout(k) = data(L,4) / abs(mean(data(1:10,4)));
end
fVout = fVout(:)
dlmwrite(strcat(directory, savefolder, yymmdd, '_fVout_', tagline, '.txt'),{fVout});
% magnitudes should be 0.6 to 0.8; any less and there's a problem with the data.

%% Check ratio, V(in), or V(out) at short and long delay times
% across your data series (for instance, for one sample as function of
% temperature or pressure).
for k = 1:nraw
    fname = fnames(k).name;
    data = dlmread(strcat(directory,savefolder,fname,'_shifted.txt'));

    % col: 3 = V(in), 4 = V(out), 5 = ratio, 6 = detector voltage
    col = 5;
    avgS = 15; % +/- local average value
    avgL = 5; % +/- local average value
    [~,index] = min(abs(data(:,2)-80));
    [~,index2] = min(abs(data(:,2)-3500));
    s80(k) =   mean(data(index-avgS:index+avgS,  col));
    s3500(k) = mean(data(index2-avgL:index2+avgL,col));

    % Get photodiode reading from raw data
    data = dlmread(strcat(directory,datafolder,fname));
    [~,index] = min(abs(data(:,2)-80));
    [~,index2] = min(abs(data(:,2)-3500));
    pd80(k) = mean(data(index-15:index+15,6));
    pd3500(k) = mean(data(index2-5:index2+1,6));
end

%% Any additional notes

% -------- END OF CODE --------