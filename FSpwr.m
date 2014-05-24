function FSpwr(cfg)


load('FSpwr.surfData.mat');

%% tests
% create outPath
oldOutPath= cfg.outputPath;
outPathExists = 1;
useNewOutPath = 0;
while outPathExists
    if exist(cfg.outputPath, 'dir') %outpath already exists try outPath+
        cfg.outputPath = [cfg.outputPath, '+'];
        useNewOutPath = 1;
    else
        outPathExists = 0;
    end
end
mkdir(cfg.outputPath);
if useNewOutPath
    disp(['ATTENTION: ', oldOutPath, ' exists.']);
end
disp(['Data will be written to ', cfg.outputPath]);

%create vars to be calculated
switch cfg.analysisType
    case 1
        cfg.N = 0;
    case 2
        cfg.power = 0;
end


%check for forbidden values
if (cfg.measure <= 0  || cfg.measure >= 5)
    error('cfg.measure must be 1, 2, 3 or 4');
end
if ~(strcmp(cfg.test,'paired')  || strcmp(cfg.test,'two-sample') )
    error('cfg.test must be paired or two-sample');
end
if ~(cfg.analysisType == 1  || cfg.analysisType == 2 || cfg.analysisType == 3)
    error('cfg.analysisType must be 1, 2 or 3');
end
if (cfg.alpha <= 0  || cfg.alpha >= 1)
    error('Alpha must be 0 < alpha < 1');
end
if ~(cfg.tail == 1 || cfg.tail == 2)
    error('tail must be 1 or 2');
end
if (cfg.power < 0 || cfg.power > 1)
    error('Power must be 0 < power < 1');
end
if (cfg.N < 0)
    error('N cannot be negative');
end
if (rem(cfg.N,1))
    error('N must be whole number');
end

%%
disp('************************');
disp('Starting Power Analysis');
disp('************************');

save(fullfile(cfg.outputPath,[brain.meas{cfg.measure},'.', cfg.test, '.configFile.mat']), 'cfg');

nHemis = 2;
nSmoothing = size(brain.smoothing,2);
if (cfg.measure == 4) %subcortical
    nHemis = 1; %set hemi to 1, because lh & rh are concatenated
    nSmoothing = 1; %set smoothing to 1, because no smoothing is applied
end

for h = 1:nHemis
    disp(brain.hemi{h});
    for s = 1:nSmoothing
        disp(brain.smoothing{s});
        switch cfg.analysisType
            case 1 % A-priori -> calc. N
                out{h,s}.reqN = aPriori(cfg.test,cfg.alpha, cfg.power, cfg.tail, cfg.difference, brain.sd{cfg.measure,h,s}, brain.rho{cfg.measure,h,s});
                out{h,s}.actualPower = postHoc(cfg.test, out{h,s}.reqN, cfg.alpha, cfg.tail, cfg.difference, brain.sd{cfg.measure,h,s}, brain.rho{cfg.measure,h,s});
                if (cfg.measure == 4) % subcortical
                    write_subcort_outputFile(cfg,brain, out);
                else
                    save_mgh(out{h,s}.reqN,fullfile(cfg.outputPath,[brain.hemi{h},'.', brain.meas{cfg.measure},'.', cfg.test, '.reqN.smoothing',brain.smoothing{s},'.mgh']),eye(4));
                    save_mgh(out{h,s}.actualPower,fullfile(cfg.outputPath,[brain.hemi{h},'.', brain.meas{cfg.measure},'.', cfg.test, '.actualPower.smoothing',brain.smoothing{s},'.mgh']),eye(4));
                end
                
            case 2 % Post-hoc
                out{h,s}.actualPower = postHoc(cfg.test, cfg.N, cfg.alpha, cfg.tail, cfg.difference, brain.sd{cfg.measure,h,s}, brain.rho{cfg.measure,h,s});
                if (cfg.measure == 4) % subcortical
                    write_subcort_outputFile(cfg,brain, out);
                else
                    save_mgh(out{h,s}.actualPower,fullfile(cfg.outputPath,[brain.hemi{h},'.', brain.meas{cfg.measure},'.', cfg.test, '.actualPower.smoothing',brain.smoothing{s},'.mgh']),eye(4));
                end
                
            case 3 % Sensitivity
                [out{h,s}.smallestDifferece, out{h,s}.smallestCohensD] = sensitivity(cfg.test,cfg.N, cfg.alpha, cfg.power, cfg.tail, brain.sd{cfg.measure,h,s}, brain.rho{cfg.measure,h,s});
                if (cfg.measure == 4) % subcortical
                    write_subcort_outputFile(cfg,brain, out);
                else
                    save_mgh(out{h,s}.smallestDifferece,fullfile(cfg.outputPath,[brain.hemi{h},'.', brain.meas{cfg.measure},'.', cfg.test, '.smallestDifference.smoothing',brain.smoothing{s},'.mgh']),eye(4));
                    save_mgh(out{h,s}.smallestCohensD,fullfile(cfg.outputPath,[brain.hemi{h},'.', brain.meas{cfg.measure},'.', cfg.test, '.smallestCohensD.smoothing',brain.smoothing{s},'.mgh']),eye(4));
                end
        end
    end
end

%write output text file
if (cfg.measure ~= 4)
    write_cortical_summary_outputFile(cfg, brain, out);
end

disp('Calculations finished');



%% functions
function reqN = aPriori(theTest,alpha, targetPower, tail, delta, sd, rho)
incN = [100,50,25,10,1];
newStartN = ones(size(sd));

disp('calculating required N');

for i=1:size(incN,2)
    powerReached = 0;
    pow = [];
    N = newStartN;
    while powerReached == 0
        switch theTest
            case 'two-sample'
                pow(end+1,:) = calcPower_2sample(N(end,:),alpha,tail,delta,sd);
            case 'paired'
                pow(end+1,:) = calcPower_paired(N(end,:),alpha,tail,delta,sd,rho);
        end
        
        lastPowCalc = pow(end,:);
        N(end+1,:) = N(end,:) + incN(i);
        % disp(N(end,1));
        if ~(all(isnan(lastPowCalc)))
            if ( all(lastPowCalc(~isnan(lastPowCalc)) > targetPower) ) %target Power reached
                powerReached = 1;
                N(end,:) = [];
            end
        end
    end
    
    newStartN = calcReqN(pow, targetPower, N);
    if (incN(i)>1)
        newStartN = newStartN-1*incN(i);
    else
        reqN = newStartN;
    end
   
end



function pow = postHoc(theTest, N, alpha, tail, delta, sd, rho)
disp('calculating actual Power');
switch theTest
    case 'two-sample'
        pow = calcPower_2sample(N,alpha,tail,delta,sd);
    case 'paired'
        pow = calcPower_paired(N,alpha,tail,delta,sd,rho);
end

function [smallestDifferece, smallestCohensD] = sensitivity(theTest,N, alpha, targetPower, tail, sd, rho)
disp('calculating smallest detectable difference');

incD = [1,.1,.01,.001, .0001];
newStartD = zeros(size(sd));

for i=1:size(incD,2)
    powerReached = 0;
    pow = [];
    D = newStartD;
    
    while powerReached == 0
        switch theTest
            case 'two-sample'
                pow(end+1,:) = calcPower_2sample(N,alpha,tail,D(end,:),sd);
            case 'paired'
                pow(end+1,:) = calcPower_paired(N,alpha,tail,D(end,:),sd,rho);
        end
        lastPowCalc = pow(end,:);
        D(end+1,:) = D(end,:) + incD(i);
        if ~(all(isnan(lastPowCalc)))
            if ( all(lastPowCalc(~isnan(lastPowCalc)) > targetPower) ) %target Power reached
                powerReached = 1;
                D(end,:) = [];
            end
        end
    end
    
    newStartD = calcReqN(pow, targetPower, D);
    if (i < length(incD) )
        newStartD = newStartD-1*incD(i);
    else
        reqD = newStartD;
    end
end

switch theTest
    case 'two-sample'
        smallestCohensD = reqD ./ sd;
    case 'paired'
        smallestCohensD = reqD ./ (sqrt(2 * sd.^2 .* (1 - rho)));
end
smallestDifferece = round(reqD .*1000) ./ 1000;
smallestCohensD = round(smallestCohensD .*1000) ./ 1000;


function reqN = calcReqN(pow, targetPower, N)
powerReached = pow > targetPower;
temp = (powerReached .* N);
temp(temp == 0) = NaN;
reqN = min(temp);


function pow = calcPower_2sample(N,alpha,tail,delta,sd)
%N per group
if (size(N,2) == 1) %if N == int make vector
    N = N * ones(size(sd));
end
dz = delta ./ sd;
d= dz .* sqrt((N.^2) ./ (2*N));
isNaN = find(isnan(d));
d(isNaN) = 0;
df = 2*N - 2;
tc=tinv(1-alpha/tail, df);
pow = 1 - nctcdf(tc, df, d );
pow(isNaN) = NaN;



function pow = calcPower_paired(N,alpha,tail,delta,sd,rho)
if (size(N,2) == 1) %if N == int make vector
    N = N * ones(size(sd));
end
dz=delta ./ sqrt(2*sd.^2 - 2*rho.*sd.^2);
d= dz .* sqrt(N);
isNaN = find(isnan(d));
d(isNaN) = 0;
df = N - 1;
tc=tinv(1-alpha/tail, df);
pow = 1 - nctcdf(tc, df, d );
pow(isNaN) = NaN;


function write_cortical_summary_outputFile(cfg, brain, out)

nSmoothing = size(brain.smoothing,2);
ck = clock;
fid = fopen(fullfile(cfg.outputPath,[ brain.meas{cfg.measure},'.',cfg.test,'.results.txt']), 'w');
fprintf(fid,'%s %u:%u\n',date,ck(4),ck(5));
fprintf(fid, 'FSpwr.m * Freesurfer power analysis * F. Liem\n');
fprintf(fid,'See %s for surface overlays in fsaverage space\n', cfg.outputPath);
fprintf(fid,'View, e.g., with "tksurfer fsaverage lh inflated -overlay lh.[...].smoothing0.mgh"\n\n');

switch cfg.analysisType
    case 1 %aPriori
        fprintf(fid,'\n**********************\n');
        fprintf(fid, 'A Priori Analysis\n');
        fprintf(fid,'**********************\n');
        fprintf(fid, '%s: %s: target Power = %f; alpha = %f; absolute difference (in mm%u) = %f; %d-tailed\n\n', brain.meas{cfg.measure},cfg.test, cfg.power, cfg.alpha,cfg.measure, cfg.difference, cfg.tail);
        
        if strcmp(cfg.test,'two-sample')
            fprintf(fid,'N is N per group!\n\n');
        end
        
        for s = 1:nSmoothing
            Ndata = [out{1,s}.reqN,out{2,s}.reqN];
            fprintf(fid, 'smoothing (mm): %s\n', brain.smoothing{s});
            fprintf(fid, 'Mean required N = %f\n', nanmean(Ndata));
            fprintf(fid, 'SD required N = %f\n', nanstd(Ndata));
            fprintf(fid, '95th percentile: required N = %f\n', prctile(Ndata,95));
            actPow = [out{1,s}.actualPower,out{2,s}.actualPower];
            fprintf(fid, 'Mean actual Power = %f\n', nanmean(actPow));
            fprintf(fid, '95th percentile: actual Power = %f\n\n', prctile(actPow,95));
        end
    case 2 %post hoc
        fprintf(fid,'\n**********************\n');
        fprintf(fid, 'Post-hoc Analysis\n');
        fprintf(fid,'**********************\n');
        fprintf(fid, '%s: %s: N = %d; alpha = %f; absolute difference (in mm%u) = %f; %d-tailed\n\n', brain.meas{cfg.measure},cfg.test, cfg.N, cfg.alpha,cfg.measure, cfg.difference, cfg.tail);
        
        if strcmp(cfg.test,'two-sample')
            fprintf(fid,'N is N per group!\n\n');
        end
       
        for s = 1:nSmoothing
            fprintf(fid, 'smoothing (mm): %s\n', brain.smoothing{s});
            actPow = [out{1,s}.actualPower,out{2,s}.actualPower];
            fprintf(fid, 'Mean actual Power = %f\n', nanmean(actPow));
            fprintf(fid, '95th percentile: actual Power = %f\n\n', prctile(actPow,95));
        end
    case 3 %sensitivity
         fprintf(fid,'\n**********************\n');
        fprintf(fid, 'Sensitivity Analysis\n');
        fprintf(fid,'**********************\n');
        fprintf(fid, '%s: %s: N = %d; alpha = %f; power = %f; %d-tailed\n\n', brain.meas{cfg.measure},cfg.test, cfg.N, cfg.alpha, cfg.power, cfg.tail);
        
        if strcmp(cfg.test,'two-sample')
            fprintf(fid,'N is N per group!\n\n');
        end
        
        for s = 1:nSmoothing
            fprintf(fid, 'smoothing (mm): %s\n', brain.smoothing{s});
            sallestDiff = [out{1,s}.smallestDifferece,out{2,s}.smallestDifferece];
            sallestCohensD = [out{1,s}.smallestCohensD,out{2,s}.smallestCohensD];
            fprintf(fid, 'Mean smallest difference (in mm%u) = %.3f\n', cfg.measure, nanmean(sallestDiff));
            fprintf(fid, 'SD smallest difference (in mm%u) = %.3f\n', cfg.measure, nanstd(sallestDiff));

            fprintf(fid, '95th percentile: smallest difference (in mm%u) = %.3f\n\n',cfg.measure, prctile(sallestDiff,95));
            fprintf(fid, 'Mean smallest Cohen''s D = %.3f\n', nanmean(sallestCohensD));
            fprintf(fid, 'SD smallest Cohen''s D = %.3f\n', nanstd(sallestCohensD));
            fprintf(fid, '95th percentile: smallest Cohen''s D  = %.3f\n\n', prctile(sallestCohensD,95));
        end
end



function write_subcort_outputFile(cfg, brain, out)
ck = clock;
fid = fopen(fullfile(cfg.outputPath,[ brain.meas{cfg.measure},'.',cfg.test,'ROIs.results.txt']), 'w');
fprintf(fid,'%s %u:%u\n',date,ck(4),ck(5));
fprintf(fid, 'FSpwr.m * Freesurfer power analysis * F. Liem\n');

switch cfg.analysisType
    case 1 %aPriori
        fprintf(fid,'\n**********************\n');
        fprintf(fid, 'A Priori Analysis\n');
        fprintf(fid,'**********************\n');
        fprintf(fid, '%s: %s: target Power = %f; alpha = %f; absolute difference (in mm3) = %f; %d-tailed\n\n', brain.meas{cfg.measure},cfg.test, cfg.power, cfg.alpha, cfg.difference, cfg.tail);
        if strcmp(cfg.test,'two-sample')
            fprintf(fid,'N is N per group!\n\n');
        end
        fprintf(fid, 'Region\trequired N\n');
        for r=1:size(out{1,1}.reqN,2)
            fprintf(fid, '%s\t%d\n', brain.asegNames{r}, out{1,1}.reqN(r));
        end
        
    case 2 %post hoc
        fprintf(fid,'\n**********************\n');
        fprintf(fid, 'Post-hoc Analysis\n');
        fprintf(fid,'**********************\n');
        fprintf(fid, '%s: %s: N = %d; alpha = %f; absolute difference (in mm3) = %f; %d-tailed\n\n', brain.meas{cfg.measure},cfg.test, cfg.N, cfg.alpha, cfg.difference, cfg.tail);
        if strcmp(cfg.test,'two-sample')
            fprintf(fid,'N is N per group!\n\n');
        end
        fprintf(fid, 'Region\tactual Power\n');
        for r=1:size(out{1,1}.actualPower,2)
            fprintf(fid, '%s\t%f\n', brain.asegNames{r}, out{1,1}.actualPower(r));
        end
        
    case 3 %sensitivity
        fprintf(fid,'\n**********************\n');
        fprintf(fid, 'Sensitivity Analysis\n');
        fprintf(fid,'**********************\n');
        fprintf(fid, '%s: %s: N = %d; alpha = %f; power = %f; %u-tailed\n\n', brain.meas{cfg.measure},cfg.test, cfg.N, cfg.alpha, cfg.power, cfg.tail);
        
        if strcmp(cfg.test,'two-sample')
            fprintf(fid,'N is N per group!\n\n');
        end
        fprintf(fid, 'Region\tsmallesd difference (mm3)\tsmallest Cohen''s D\n');
        for r=1:size(out{1,1}.smallestDifferece,2)
            fprintf(fid, '%s\t%.3f\t%.3f\n', brain.asegNames{r}, out{1,1}.smallestDifferece(r),out{1,1}.smallestCohensD(r));
        end       
end

