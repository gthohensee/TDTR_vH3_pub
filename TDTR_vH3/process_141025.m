
%% Define the directory and filenames
directory = '/Users/gregoryhohensee/Documents/MATLAB/TDTR_vH2_1C/TDTR_vH2';
yymmdd = '141025'; % e.g., 140227 = February 27th, 2014.
datafolder = strcat('/',yymmdd,'/');      % folder in directory where raw data is kept
savefolder = strcat('/',yymmdd,'_edit/'); % folder in directory where analyzed data is saved
datadir = strcat(directory, datafolder);
savedir = strcat(directory, savefolder);

tagline = 'acoustics'; % common prefix for your series of TDTR data files.

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
  
  % SetPhase_vH picks t0 and sets the phase. PickEchoes reads acoustics.
  tmax = 60;
  t_window = [-10 10];
  [t0(k)] = set_t0(data, t_window); % Set t0
  t_window = [-20 80];
  [delphases(k),phase_shift(k),fitparam(k,:)] = SetPhase_vH2(data,t0(k),fname,savedir,t_window);
  %echotime = PickEchoes(data,tmax,t0(k),t_window);
  %echoes(k,1:length(echotime)) = echotime;
  %clear echotime;
  % AlSiO2: 28.25 ps
  % Al_nSi: 24 to 25 ps
  % Al_SiDAC: 25 ps
  % k = 3 Afu* sample had huge change in V(out), dubious.
end
%% Save delphases and echoes
shifts = horzcat(t0, phase_shift);
dlmwrite(strcat(directory, savefolder, yymmdd, '_delphases_', tagline, '.txt'),    {delphases});
dlmwrite(strcat(directory, savefolder, yymmdd, '_shifts_', tagline, '.txt'),       {shifts});
dlmwrite(strcat(directory, savefolder, yymmdd, '_echotimes_', tagline, '.txt'),    {echoes});
dlmwrite(strcat(directory, savefolder, yymmdd, '_Voutfitparams_', tagline, '.txt'),{fitparam});
%% Check the percent change in Vout over 4 ns; retain the sign of Vout.
clear fVout;
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
    col = 3;
    avgS = 5; % +/- local average value
    avgL = 1; % +/- local average value
    [~,index] = min(abs(data(:,2)-50));
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