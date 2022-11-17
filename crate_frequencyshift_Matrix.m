function crate_frequencyshift_Matrix
try

EigenModes_Data=dir('SpinDynamicsResponse_PdB_Box*_*_gamma0*.mat');

NoOfEigenModesMat=length(EigenModes_Data);

% DimensionlessCoefficient_Matrix=[];
FrequencyShift_Matrix=[];


save('FrequencyShift_box20_0_gamma000_gamma020.mat','FrequencyShift_Matrix');

for ii=1:NoOfEigenModesMat
      ii 
      disp(EigenModes_Data(ii).name);

load(EigenModes_Data(ii).name,'A','B','EigenValue','indx','EigenVector10','gamma','gamma1');
     
      q=gamma; % step length is 0.02;

      column_number=round((q/0.02))+1

      load('FrequencyShift_box20_0_gamma000_gamma020.mat','FrequencyShift_Matrix');
     
        
      FrequencyShift_Matrix(column_number)=EigenValue(1,1);

  
     save('FrequencyShift_box20_0_gamma000_gamma020.mat','FrequencyShift_Matrix');
     disp(' successivly save ');
end

     FrequencyShift_Matrix
     
catch error

     disp(getReport(error))
     exit(1)

end
