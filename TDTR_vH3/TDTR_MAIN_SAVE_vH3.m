%% Subfunction of TDTR_MAIN: Save and print figure from MAIN results
figsens = 202;
figfit = 203;

% Save the workspace
save(strcat(datadir,'Results_', fname(1:end-4),'.mat'))

% define the data and model arrays for the results figure
if ~doughnut
    switch sigfit
        case 1
            plot_data = Vin_data / Vin_data(Zind);
            plot_model = real(deltaR_model);
            plot_model = plot_model / (plot_model(Zind)/N); % see MANUALFIT for N
            ytext = 'normalized V(in)';
            fittext = 'FITVin_';
        case 2
            plot_data = Vout_data / Vout_data(Zind);
            plot_model = imag(deltaR_model);
            plot_model = plot_model / (plot_model(Zind)/N);
            ytext = 'normalized V(out)';
            fittext = 'FITVout_';
        otherwise
            plot_data = ratio_data;
            plot_model = ratio_model;
            ytext = '-V(in)/V(out)';
            fittext = 'FIT_';
    end
else
    [~,I] = min(abs(offset));              % get index of r = 0 data point
    Norm = mean(abs(Vout_data(I-1:I+1)));  % take 3pt average of Vout data around r = 0
    plot_data = abs(Vout_data/Norm);      % remove negative sign, rescale to match DOUGHNUT.

    Vout_model = imag(deltaR_model);        
    plot_model = abs(Vout_model) / max(abs(Vout_model));
    ytext = 'normalized V(out)';
    fittext = 'FITVout_boffset_'; 
end

% Create a figure for the last fit
figure(figfit)
clf;

if importdata
    w0 = sqrt(r_pump*r_probe);
    if ~doughnut
        hd = loglog(tdelay,plot_data,'ko',...
                   'MarkerSize',8.0,'LineWidth',1);
        hold on;
        hm = loglog(tdelay,plot_model,'r-',...
                   'LineWidth',2);
    else
        hd = plot(offset*1e-6/w0 - V,plot_data / N,'ko',...
                   'MarkerSize',8.0,'LineWidth',1);
        hold on;
        hm = plot(offset*1e-6/w0,plot_model,'r-',...
                   'LineWidth',2);
    end

else % if importdata is false, then we didn't fit anything,
     % so there's no data points to plot -- just the model.
    if ~doughnut
        hm = loglog(tdelay,plot_model,'r-',...
                   'LineWidth',2);
    else
        hm = plot(offset*1e-6/w0,plot_model,'r-',...
                   'LineWidth',2);
    end

end
fontsize = 16;

% define condtext
if P0 ~= -1
    if exist('frac','var') && length(r_pump) > 1
        condtext = sprintf('P = %0.1f GPa, dR_{pump} = %i%%, Z = %0.2d, t_{fit} = %i ps',P0,frac*100,Z,Zdelay);
    else
        condtext = sprintf('P = %0.1f GPa, Z = %0.2d, t(fit) = %i ps',P0,Z,Zdelay);
    end
elseif T0 ~= -1
    if exist('frac','var') && length(r_pump) > 1
        condtext = sprintf('T = %0.1f K, dR_{pump} = %i%%, Z = %0.2d, t_{fit} = %i ps',T0,frac*100,Z,Zdelay);
    else
        condtext = sprintf('T = %0.1f K, Z = %0.2d, t(fit) = %i ps',T0,Z,Zdelay);
    end
else
    condtext = sprintf('dR_{pump} = %i%%, Z = %0.2d, t_{fit} = %i ps',frac*100,Z,Zdelay);
end

% define axes
if sigfit ~= 2
    set(gca, 'YTick', [0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 10 20 30 50]);
    axis([tdelay_min 5e-9 min(0.1,min(plot_data)) 2*max(plot_data)])
else
    set(gca,'YScale','linear');
    set(gca, 'YTick', [0 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2]);
    if ~doughnut, axis([tdelay_min 5e-9 min(0,min(plot_data)) 1.2]);
    else axis([-4 4 -0.2 1.2]); end

end

% define labels
if doughnut, xlabel('Offset x/w0', 'FontSize',fontsize); 
else xlabel('time delay (ps)','FontSize',fontsize); end

ylabel(ytext, 'FontSize',fontsize)
title(condtext,'FontSize',fontsize)

if ~doughnut
    set(gca, 'XTick', [1e-11 2e-11 5e-11 1e-10,2e-10,5e-10,1e-9,2e-9,5e-9,1e-8]);
    set(gca, 'XTickLabel', [10 20 50 100, 200, 500, 1000, 2000, 5000, 1e4]);
    set(gca, 'XMinorTick', 'off');
    set(gca, 'YMinorTick', 'off');
end

set(gca, 'TickLength' , [.02 .02]);
set(gca,'FontSize',fontsize);
%% Provide a summary of final fit parameters in the figure.

% Warning: uses "ii" index from the analyze script's for loop!
% Compose fitstr and update LCTsol to match Xsol
XsolIJ = XguessIJ;
XsolIJ(:,1) = Xsol';
LCTEtext = genLCTEtext(LCTE, XsolIJ, stack(ii,:));

% write LCTEstr contents to a text box in the figure
pBox = annotation('textbox',[0.15,0.15,0.8,0.3]);
set(pBox,'Units','characters')
set(pBox,'HorizontalAlignment','left')
set(pBox,'FontSize',12)
set(pBox,'String',LCTEtext);
set(pBox,'LineStyle','none');

% also record the data file used
dBox = annotation('textbox',[0.15,0.85,0.8,0.07]);
set(dBox,'Units','characters')
set(dBox,'Interpreter','none')
set(dBox,'HorizontalAlignment','left')
set(dBox,'FontSize',12)
set(dBox,'String',{'Data file:'; fname(1:end)});
set(dBox,'LineStyle','none');

%% Save the figure to .fig and .eps files
saveas(figfit, strcat(datadir,fittext,fname(1:end-4),'.fig'))
print(figfit,'-depsc',strcat(datadir,fittext,fname(1:end-4),'.eps'))

if senseplot
    tag = input('Name the sensitivity plot: ','s');
    save(strcat(datadir,'/SENS_', tag, '.mat'))
    figure(figsens)
    saveas(figsens, strcat(datadir,'/SENS_',tag,'.fig'))
    print(figsens,'-depsc',strcat(datadir,'/SENS_',tag,'.eps'))
end