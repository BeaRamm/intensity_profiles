function rawdata=read_bathtubs_profile(folder,filename)

[~,datastr,~]=xlsread([folder,'\',filename]);
rawdata=zeros(numel(datastr)-1,2);
for j=2:numel(datastr)
    rawdata(j-1,:)=str2num(datastr{j});
end