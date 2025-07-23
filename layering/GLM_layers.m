%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Script to run GLM on layers
% This script is similar to GLM_layers_normalized.m, but no
% z-transformation (normalization using mean and standard deviation) is applied.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Locate all layers.mat fles i.e., outputs of Layering_Automatization.m
fileContent = fileread('unique_no_PA.txt');
folders = splitlines(fileContent);
folders = folders(1:end-1); % remove last empty entry
folders = regexp(folders, '^(.*?/func)/', 'tokens', 'once');
folders = vertcat(folders{:});

%% First apply layer-wise GLM on BOLD data
N_sub = 6; % number of participants
N = 30; % number of layers 
con_idx = 1; % index for memory_vs_math contrast 
N_vol = 150; % number of time-points in each run
run_indx_max = [1:150;151:300;301:450;451:600;601:750;751:900;901:1050;...
    1051:1200;1201:1350]; % maximum possible number of time-points if concatenation was done across multiple sessions. There are no more than 9 runs in total for none of the subjects.

for n = 1:N_sub
    cd(folders{n})
    load('layers_BOLD_merged.mat');    
    cd('first_level/bold')
    load('SPM.mat');
    
    % Multiply spm design matrix withlayer values

    runs = length(layers_bold_per_sub{1,1})./N_vol; % the layers_BOLD_merged.mat file is based on concatenation of all functional runs. Divide by time points to get the number of runs
   
    for run = 1:runs
        sz = size(layers_bold_per_sub{1,1}(:,run_indx_max(run,:),:));
        tmp = reshape(layers_bold_per_sub{1,1}(:,run_indx_max(run,:),:),[sz(1)*sz(3) ,sz(2)]);   
        tmp =  tmp - (SPM.xX.K(run).X0*(SPM.xX.K(run).X0.'*tmp.')).';
        layers_bold_per_sub{1,1}(:,run_indx_max(run,:),:) = reshape(tmp,sz);  
    end

   % run GLM on layers 

   e = length(con_idx);  
   W = SPM.xX.W;
   GLM_use = SPM.xX.pKX;
   con = SPM.xCon(con_idx(e)).c; 
   Subfield = size(layers_bold_per_sub{1,1},3); % 1=sub, 2=ca1, 3=ca2, 4=ca3

   T = zeros(Subfield,N); 
   Tcrit = zeros(Subfield,1);
   pmax = Tcrit;
   con_array = T;

   for f = 1: Subfield
       Subregion = layers_bold_per_sub{1,1}(:,:,f)';
       KWY = spm_filter(SPM.xX.K,W*Subregion);
       b = GLM_use*KWY;
       res = spm_sp('r',SPM.xX.xKXs,KWY);
       ResSS = sum(res.^2);                    %-Residual SSQ
       ResMS = ResSS / SPM.xX.trRV;
       Vc  = con'*SPM.xX.Bcov*con;
       SE  = sqrt(ResMS*Vc);    
       beta_index = find(abs(con) > 0);
       beta_use   = b(beta_index,:); 

       for j=1:size(beta_use,1)
           con_array(f,:) = con_array(f,:) + con(beta_index(j)) * beta_use(j,:);
       end

       T(f,:) = con_array(f,:)./SE;
       alpha = 0.05; %p-value

       %correction = {'FWE','FDR','none'}; % choose and change the next line based on choise 

       correction = 'NONE'; 

       if strcmp(correction,'FWE') == 1
          Tcrit(f) = spm_uc(alpha,[1 SPM.xX.erdf],'T',SPM.xVol.R,1,numel(b));

      elseif strcmp(correction,'FDR') == 1
           p = 2 * (1 - spm_Tcdf(T(f,:), SPM.xX.erdf));
           p = spm_P_FDR(p);
           tmp = T(f,:);
           Tcrit(f) = min(tmp(p<alpha));
           clear tmp
           if isempty(Tcrit(f))
              Tcrit(f) = nan;
           end
           pmax(f) = max(p(p<alpha));
           if isempty(pmax(f))
              pmax(f) = 1;
           end
      elseif strcmp(correction,'NONE') == 1
           Tcrit(f) = spm_u(alpha,[1 SPM.xX.erdf],'T');
           pmax(f) = alpha;
       end
   
       if  f == 4 % if subfield is Sub, CA3 -> don't include SRLM estimation 
           con_array(f,1:10) = nan;
           T(f,1:10) = nan;
      end

      results = struct('con_array',con_array(f,:),'T',T(f,:),'Tcrit',Tcrit(f),'pmax',pmax(f));
      save(strcat('results','_',string(f)),'-struct','results');
   end

cd ../../../../   
end 
clearvars -except folders N_vol N_sub con_idx N run_indx_max

%% Repeat the above for VASO

for n = 1:N_sub
    cd(folders{n})
    load('layers_VASO_merged.mat');    
    cd('first_level/vaso')
    load('SPM.mat');

    % Multiply spm design matrix withlayer values

    runs = length(layers_vaso_per_sub{1,1})./N_vol; % the layers_BOLD_merged.mat file is based on concatenation of all functional runs. Divide by time points to get the number of runs
   
    for run = 1:runs
        sz = size(layers_vaso_per_sub{1,1}(:,run_indx_max(run,:),:));
        tmp = reshape(layers_vaso_per_sub{1,1}(:,run_indx_max(run,:),:),[sz(1)*sz(3) ,sz(2)]);   
        tmp =  tmp - (SPM.xX.K(run).X0*(SPM.xX.K(run).X0.'*tmp.')).';
        layers_vaso_per_sub{1,1}(:,run_indx_max(run,:),:) = reshape(tmp,sz);  
    end

    % run GLM on layers 

    e = length(con_idx);  
    W = SPM.xX.W;
    GLM_use = SPM.xX.pKX;
    con = SPM.xCon(con_idx(e)).c; 
    Subfield = size(layers_vaso_per_sub{1,1},3); % 1=sub, 2=ca1, 3=ca2, 4=ca3

    T = zeros(Subfield,N); 
    Tcrit = zeros(Subfield,1);
    pmax = Tcrit;
    con_array = T;

    for f = 1: Subfield
       Subregion = layers_vaso_per_sub{1,1}(:,:,f)';
       KWY = spm_filter(SPM.xX.K,W*Subregion);
       b = GLM_use*KWY;
       res = spm_sp('r',SPM.xX.xKXs,KWY);
       ResSS = sum(res.^2);                    %-Residual SSQ
       ResMS = ResSS / SPM.xX.trRV;
       Vc  = con'*SPM.xX.Bcov*con;
       SE  = sqrt(ResMS*Vc);    
       beta_index = find(abs(con) > 0);
       beta_use   = b(beta_index,:); 

       for j=1:size(beta_use,1)
           con_array(f,:) = con_array(f,:) + con(beta_index(j)) * beta_use(j,:);
       end

       T(f,:) = con_array(f,:)./SE;
       alpha = 0.05; %p-value
       %correction = {'FWE','FDR','none'}; % choose and change the next line based on choise 

       correction = 'NONE';
       if strcmp(correction,'FWE') == 1
          Tcrit(f) = spm_uc(alpha,[1 SPM.xX.erdf],'T',SPM.xVol.R,1,numel(b));

      elseif strcmp(correction,'FDR') == 1
           p = 2 * (1 - spm_Tcdf(T(f,:), SPM.xX.erdf));
           p = spm_P_FDR(p);
           tmp = T(f,:);
           Tcrit(f) = min(tmp(p<alpha));
           clear tmp
           if isempty(Tcrit(f))
              Tcrit(f) = nan;
           end
           pmax(f) = max(p(p<alpha));
           if isempty(pmax(f))
              pmax(f) = 1;
           end
      elseif strcmp(correction,'NONE') == 1
           Tcrit(f) = spm_u(alpha,[1 SPM.xX.erdf],'T');
           pmax(f) = alpha;
       end
    
       if  f == 4 % if subfield is Sub, CA3 -> don't include SRLM estimation 
           con_array(f,1:10) = nan;
           T(f,1:10) = nan;
       end

      results = struct('con_array',con_array(f,:),'T',T(f,:),'Tcrit',Tcrit(f),'pmax',pmax(f));
      save(strcat('results','_',string(f)),'-struct','results');
   end

cd ../../../../   
end 