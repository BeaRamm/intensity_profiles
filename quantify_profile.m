function [CHdata,discard]=quantify_profile(CHdata)
%%%%% function [CHdata,discard]=quantify_profile(CHdata) %%%%%
% called by analyze_bathtubs_profiles to identify the center 
% and the edges of the microcompartments
% INPUT: data structure as used in analyze_bathtubs_profiles 
% 
% For detailed information: Ramm et al, Nat. Commun. 2018
%
% Jonas Mücksch 2018


x=CHdata.reduced_normprofile_unitlength(:,1);
prof=CHdata.reduced_normprofile_unitlength(:,2);
% p=polyfit(x,prof,2);
% CHdata.x0=-p(2)/(2*p(1))

% fit a quadratic function without linear contribution to 
% profile
[p,proffit]=fit_quadratic_function_center0(x,prof);
% plot data
figure('Units','centimeters','Position',[40,25,15,15],'Color',[1,1,1]);
plot(x,prof);
hold on
plot(x,proffit);
set(gca,'XLim',[-0.5,.5]);


CHdata.x0=0; %profile center (profile was aligned to be 
			%centered at x0=0
CHdata.curvature=p(1);
CHdata.int=(0.5*(prof(1)+prof(end))+sum(prof(2:end-1)))*abs(x(1)-x(2)); 
CHdata.y0=p(2); %integrated profile

if CHdata.x0<x(1) || CHdata.x0>x(end)
    discard=true;
else
    discard=false;
end
end

function [p,yfit]=fit_quadratic_function_center0(x,y)
options=optimset('Algorithm','trust-region-reflective',...
                'Maxiter',1E5,'TolX',1E-9,'TolFun',1E-9,...
                'Display','none');
pstart=[0,min(y)];
pmin=[-10,0];
pmax=[10,max(y)];
p=lsqnonlin(@fitfunx2,pstart,pmin,pmax,options,x,y);
yfit=fitfunx2(p,x,zeros(size(y)));
end

function dy=fitfunx2(p,x,y)
dy=p(1)*x.^2+p(2)-y;
end