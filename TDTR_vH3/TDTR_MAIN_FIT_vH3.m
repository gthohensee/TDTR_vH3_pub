%% Subroutine of MAIN.
%--------------Perform Fit--------------------------
% Fitting assigns values to Xsol, Z, deltaR_model, ratio_model, LCTE.
% LCTE changes if eta changes by changing L, or if the entire LCTE
% changes from self-consistent temperature changes.
if manualfit
    if ~doughnut
        switch sigfit
            case 1, fprintf('Manual fitting to V(in)(t) normalized at %i ps...\n',tdelay(Zind)*1e12);
            case 2, fprintf('Manual fitting to V(out)(t) normalized near %i ps...\n',tdelay(Zind)*1e12);
            otherwise fprintf('Manual fitting to ratio -V(in)/V(out)...\n');
        end
    else
        fprintf('Manual fitting beam offset V(out) data (normalized to unity)...\n');
    end

    [Xsol,Z,deltaR_model,ratio_model,LCTE,T_adj,N,V] = ...
              TDTR_MANUALFIT_vH3(XguessIJ, datparams,...
            sysparams, calparams, matparams, Tparams);

    XguessIJ(:,1) = Xsol;

else
    if doughnut, warning('Sorry, autodoughnut is not coded yet.'); break; end
    Xguess = XguessIJ(:,1);

    fprintf('Please wait for automatic fitting...\n');
    tic
    Xsol = fminsearch(@(X) TDTR_FIT_vH3(X, Xij, datparams,...
                           sysparams, calparams, matparams, Tparams),...
                           Xguess,options);

    % fminsearch doesn't output anything but Xsol, so get the rest here.                
    [Z,deltaR_model,ratio_model,LCTE,T_adj]=...
        TDTR_FIT_vH3(Xsol, Xij,datparams,sysparams,calparams,matparams,Tparams);
    toc
end
fprintf('Data fit completed\n');