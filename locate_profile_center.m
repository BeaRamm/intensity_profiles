function data=locate_profile_center(data)
%%%%% function data=locate_profile_center(data) %%%%%
% called by analyze_bathtubs_profiles to identify the center 
% and the edges of the microcompartments
% INPUT: data structure as used in analyze_bathtubs_profiles 
% 
% For detailed information: Ramm et al, Nat. Commun. 2018
%
% Jonas Mücksch 2018

discard=false;
%% determine where the bathtub begins
prof=data.greenCH.rawdata(:,2); %profile to be analyzed
prof=prof/max(prof);            %normalize

diffprof=diff(prof);            %derivative to see where the signal increases
[~,outind]=findpeaks(abs(diffprof),'MinPeakdistance',round(numel(prof)*1/2),...
                    'MinPeakHeight',0.1);   %find max increase, i.e. where bathtub starts
%deal with outliers and special cases
if numel(outind)<2
    outind=[1,numel(prof)];
elseif numel(outind)>2  %too many peaks found, take extreme
    outind=[min(outind),max(outind)];
end
if (outind(1)>numel(prof)/4) && (outind(2)<numel(prof)*3/4) %peaks too close to the center
    outind=[1,numel(prof)];
elseif (outind(1)>numel(prof)/4) && (outind(2)>numel(prof)*3/4)
    outind(1)=numel(prof)-outind(2);
elseif (outind(1)<numel(prof)/4) && (outind(2)<numel(prof)*3/4)
    outind(2)=numel(prof)-outind(1);
end
prof=prof(outind(1):outind(2)); %clip the profile according to determined limits


%% determine center position of MinD profile
smprof=smooth(prof,20); %smooth profile to kill some local noise
centerpos_guess=round(numel(prof)/2);   %first guess is middle of the chamber
sel=[-1,1]*round(numel(prof)/8)+round(centerpos_guess); %cut out ROI
xsel=(sel(1):sel(2))';
prof_sel=prof(xsel);    %profile in ROI
[centerpos,~]=find_center_parabola(xsel,prof(xsel)); %fit quadratic function to profile in ROI
if centerpos>sel(2) || centerpos<sel(1)     %extreme in selected range not found
    centerpos=round(numel(prof)/2);
end

%% determine edges as maxima of MinD profile
p=polyfit(1:numel(prof),prof',4); %smooth profile by approaching with 4th order polynomial
[~,edgepos_guess(1)]=max(polyval(p,1:ceil(numel(prof)/2))); %left max is guessed as first max of the polynomial
if abs(edgepos_guess(1)-numel(prof)/4)>numel(prof)/4 %what if found position is too close to center
    edgepos_guess(1)=round(numel(prof)/4);
end
[~,edgepos_guess(2)]=max(polyval(p,(ceil(numel(prof)/2)+1):numel(prof))); %right max is guessed as other max of the polynomial
edgepos_guess(2)=edgepos_guess(2)+ceil(numel(prof)/2);
if abs(edgepos_guess(2)-numel(prof)*3/4)>numel(prof)/4 %what if found position is too close to center
    edgepos_guess(1)=round(numel(prof)*3/4);
end

for i=1:numel(edgepos_guess) %fit again quadratic function to find center of maxima
    xsel=(edgepos_guess(i)+(-20:20))';                      %ROI
    xsel(xsel<=0)=[]; xsel(xsel>numel(prof))=[];            %exclude some cases
    [val,pos_candidate]=findpeaks(smooth(prof(xsel),5),'MinPeakDistance',range(xsel)/2);
    if numel(val)>1
        [~,ind]=max(val);
        edgepos_guess(i)=pos_candidate(ind)+xsel(1)-1;
    elseif numel(val)==1
        edgepos_guess(i)=pos_candidate+xsel(1)-1;
    elseif isempty(val)
        discard=true;
    end
    
end
edgepos=edgepos_guess;

%% store all profiles to the data structure
data.discard=discard;
if ~discard
    data.greenCH.edgepos=edgepos;
    data.greenCH.centerpos=centerpos;
    data.greenCH.outind=outind;
    data.greenCH.sel=round((edgepos(1):edgepos(2)))';
    xpos=(data.greenCH.sel-centerpos);
    
    normprof=generate_normprof(xpos,prof(edgepos(1):edgepos(2)));
    
    
    
    data.greenCH.reduced_normprofile=[xpos,...
                                normprof];
    data.greenCH.reduced_normprofile_unitlength=[(linspace(-0.5,0.5,numel(data.greenCH.sel)))',...
                                normprof];

    redprof=data.redCH.rawdata(:,2);
    redprof=redprof/max(redprof); 
    redprof=redprof(outind(1):outind(2));
    normprof=generate_normprof(xpos,redprof(edgepos(1):edgepos(2)));
    data.redCH.reduced_normprofile=[xpos,...
                                normprof];
    data.redCH.reduced_normprofile_unitlength=[(linspace(-0.5,0.5,numel(data.greenCH.sel)))',...
                                normprof];

%% plot if interested 
figure('Units','centimeters','Position',[25,25,15,15],'Color',[1,1,1]);
plot(prof/max(prof));
hold on
plot(redprof/max(redprof));plot(centerpos*[1,1],[0,1],'-k')
plot(edgepos(1)*[1,1],[0,1],'-k')
plot(edgepos(2)*[1,1],[0,1],'-k')
title(data.file,'Interpreter','none');
end
    
end

function normprof=generate_normprof(x,y)
w=5;
% x1=mean(x(1:w));
% y1=mean(y(1:w));
% x2=mean(x(end-(1:w)+1));
% y2=mean(x(end-(1:w)+1));
% bl=(y2-y1)/(x2-x1)*(x-x1)+y1;

xf=[x(1:w),x(end-(1:w)+1)];
yf=[y(1:w),y(end-(1:w)+1)];
p=polyfit(xf,yf,1);

normprof=y-polyval(p,x)+max(polyval(p,x));

end

function [centerpos,p]=find_center_parabola(x,y)
p=polyfit(x,y,2);
centerpos=-p(2)/(2*p(1));
end