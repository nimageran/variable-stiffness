function z=bucklingMoment(T)

% T=round(T);
% display(T);

global NFE;
if isempty(NFE)
    NFE=0;
end


% Updating NFE number
NFE=NFE+1;

%%%%%% Write variables into TEXT2.TXT file 
dlmwrite('proposed_T.txt',T,'\t')

%%%%%% Writing the outputs by running the python code
disp('Start python')
system(['abaqus cae nogui=pythonScript.py'])
disp('Python is finished')

newpath=['D:\optLaminatedComp\JOB\',mat2str(NFE)];

bucklingMOMENT=importdata(fullfile( newpath,'bucklingMOMENT.txt'));
z=-abs(bucklingMOMENT);
end