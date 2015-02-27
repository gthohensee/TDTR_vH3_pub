function [t_sub,f_sub] = Subtract_Background(t,f)
%Subtract_Background - subtracts spline background from data set (t,f).
%                      For instance, remove thermal signal from TDTR data.
% Part of the TDTR_vH2 package, 14-July-2014.

%Plot Data to be Background-subtracted, Zoom to Region-of-Interest
fig1=figure(44);
plot(t,f,'-o')
axis([min(t) max(t) min(f) max(f)])

question = input('Hey, does this look promising? Hit RETURN if yes, 0 if no: \n');
if ~isempty(question)
    if question == 0
        t_sub = [];
        f_sub = [];
        return; 
    end
end

fprintf('Zoom on the area you want to analyze...the press button to continue\n')
pause

keepgoing=1;
%
while keepgoing==1;
    fprintf('Click on waypoints for cubic spline, press RETURN to stop\n')
    [t_spline,f_spline]=ginput;
    
    %EXTRACT DATA
    %Choose range of time delays to fit
    t_min=min(t_spline);
    t_max=max(t_spline);
    [t_chopped,f_chopped]=extract_interior(t,f,t_min,t_max);
    
    %generate spline background predictions
    [f_BK]=spline(t_spline,f_spline,t_chopped);
    t_BK=t_chopped; %obviously the points must match up.
    
    %plot with spline
    fig2=figure(45);
    plot(t_chopped,f_chopped,'b-o',t_BK,f_BK,'g-')
    title('spline fit')
    
    t_sub=t_chopped;
    f_sub=f_chopped-f_BK;
    
    fig3=figure(46);
    plot(t_sub,f_sub,'-o')
    title('subtracted data')
    
    test = input('Again? yes = 1: ');
    if test
        keepgoing = 1;
    else
        keepgoing = 0;
    end
    
end

close(44, 45, 46)