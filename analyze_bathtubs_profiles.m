function analyze_bathtubs_profiles(folder)
%%%%% function analyze_bathtubs_profiles(folder) %%%%%
% function analyzes time-averaged intensity profiles in
% microcompartments. 
% INPUT: folder is a string for the target folder in which the 
% source profiles are located. Source files are saved as csv files. These profiles are obtained 
% using FIJI
%
% Extracted data sets are stored in a data structure
%
% Note the convention for filenames: 
% *C2* corresponds to MinD-GFP, *C1* correspopnds to the  red % channel; for streptavidin it is the inverse
% For detailed information: Ramm et al, Nat. Commun. 2018
%
% Jonas Mücksch 2018


%% locate all files in the folder
%get all MinD-GFP profiles
files_green=dir([folder,'\*C2*.csv']);
% check for the corresponding images in the red channel
for i=1:numel(files_green)  
    %check for red channel counter file
    tempfilename=files_green(i).name;
    tempfilename=strrep(tempfilename,'C2','C1');
    files_red(i)=dir([folder,'\',tempfilename]);
    if isempty(files_red(i))
        fprintf(['file %s ', '%s\n'], files_green(i).name, 'has no counter file!');
    end
end

if ~isempty(strfind(files_red(1).name,'strep')) 
%take care of the streptavidin case:
% swap red and green, because here streptavidin was labeled in 
% the green channel
    tempfiles=files_green;
    files_green=files_red;
    files_red=tempfiles;
    clear tempfiles;
end

%% initialize figures
f=figure('Units','centimeters','Position',[2,2,16,16],'Color',[1,1,1]);
for i=1:4
    ax(i)=subplot(2,2,i); hold(ax(i),'on'); set(ax(i),'Box','on');
end

f=figure('Units','centimeters','Position',[18,2,16,16],'Color',[1,1,1]);
for i=1:4
    ax2(i)=subplot(2,2,i); hold(ax2(i),'on'); set(ax2(i),'Box','on');
end
xlabel(ax(3),'unit length');
ylabel(ax(3),'norm MinD signal [a.u.]');
xlabel(ax(4),'unit length');
ylabel(ax(3),'norm signal [a.u.]');

set(ax(3),'XLim',[-.5,.5],'YLim',[0,1]);
set(ax(4),'XLim',[-.5,.5],'YLim',[0,1]);
set(ax2(1),'XLim',[-.5,.5],'YLim',[-1.5,2]);
set(ax2(2),'XLim',[-.5,.5],'YLim',[-1.5,2]);
set(ax2(3),'XLim',[-.5,.5],'YLim',[0.5,1]);
set(ax2(4),'XLim',[-.5,.5],'YLim',[0.5,1]);

xlabel(ax2(1),'MinD min position');
xlabel(ax2(2),'pushed extreme position');
xlabel(ax2(3),'MinD min position');
xlabel(ax2(4),'pushed extreme position');
ylabel(ax2(1),'MinD profile curvature [a.u.] ???');
ylabel(ax2(2),'pushed profile curvature [a.u.] ???');
ylabel(ax2(3),'MinD AUC');
ylabel(ax2(4),'pushed AUC');

%% extract and analyze files
w=waitbar(0,'analyzing profiles');
for i=1:numel(files_green)
    fprintf('%s\n', files_green(i).name);
    %get green and red channel profiles
    	data.greenCH.rawdata=read_bathtubs_profile(folder,files_green(i).name);
    	data.redCH.rawdata=read_bathtubs_profile(folder,files_red(i).name);
    
    data.folder=folder;                 %save folder
    data.file=files_green(i).name;  	%save file name
    
    %find the center and boundaries of the profile
    data=locate_profile_center(data);   

    if data.discard
        fprintf(['MinD profile %s ', '%s\n'], files_green(i).name, 'has no proper minimum!');
        continue
    end
    [data.greenCH,data.discard]=quantify_profile(data.greenCH);
    if data.discard
        fprintf(['MinD profile %s ', '%s\n'], files_green(i).name, 'has no proper minimum!');
        continue
    end
    [data.redCH,~]=quantify_profile(data.redCH);

     
    newfilename=files_green(i).name;
    [~,filecore,~]=fileparts(newfilename);
    newfilename=[folder,'\',filecore,'_data.mat'];
    save(newfilename,'data');

    plot(ax(1),data.greenCH.rawdata(:,1),data.greenCH.rawdata(:,2));
    plot(ax(2),data.redCH.rawdata(:,1),data.redCH.rawdata(:,2));
    plot(ax(3),data.greenCH.reduced_normprofile_unitlength(:,1),...
        data.greenCH.reduced_normprofile_unitlength(:,2));
    plot(ax(4),data.redCH.reduced_normprofile_unitlength(:,1),...
        data.redCH.reduced_normprofile_unitlength(:,2));

    plot(ax2(1),data.greenCH.x0,data.greenCH.curvature,'x'); 
    plot(ax2(2),data.redCH.x0,data.redCH.curvature,'x'); 
    plot(ax2(3),data.greenCH.x0,data.greenCH.int,'x');  
    plot(ax2(4),data.redCH.x0,data.redCH.int,'x');
    
    w=waitbar(i/numel(files_green),w);
end
delete(w);

xlabel(ax(3),'unit length');
ylabel(ax(3),'norm MinD signal [a.u.]');
xlabel(ax(4),'unit length');
ylabel(ax(3),'norm signal [a.u.]');

set(ax(3),'XLim',[-.5,.5],'YLim',[0,1]);
set(ax(4),'XLim',[-.5,.5],'YLim',[0,1]);
set(ax2(1),'XLim',[-.5,.5],'YLim',[-1.5,2]);
set(ax2(2),'XLim',[-.5,.5],'YLim',[-1.5,2]);
set(ax2(3),'XLim',[-.5,.5],'YLim',[0.5,1]);
set(ax2(4),'XLim',[-.5,.5],'YLim',[0.5,1]);


xlabel(ax2(1),'MinD min position');
xlabel(ax2(2),'pushed extreme position');
xlabel(ax2(3),'MinD min position');
xlabel(ax2(4),'pushed extreme position');
ylabel(ax2(1),'MinD profile curvature [a.u.]');
ylabel(ax2(2),'pushed profile curvature [a.u.]');
ylabel(ax2(3),'MinD AUC');
ylabel(ax2(4),'pushed AUC');


end
