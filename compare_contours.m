function compare_result = compare_contours(dicom_path,manual_contour_path,auto_contour_path)
%COMPARE_CONTOURS Compare manual rawn contours with auto contours
%   COMPARE_CONTOURS(DICOM_PATH,MANUAL_CONTOUR_PATH,AUTO_CONTOUR_PATH)
%
%   Copyright: Imaging Research, Sunnybrook Health Sciences Centre, Toronto, ON, Canada
%   Author: Perry Radau, perry.radau@gmail.com
%   Version: 0.9   
%   Date: 2009/05/11

%   OUTPUT: compare_result is a struct with fields of
%-----------------------------------------------------
%   auto_number_i:  total number of auto inside contours
%   auto_number_o:  total number auto outside contours;
%   manual_number_i: total number manual inside contours
%   manual_number_o: total number manual outside contours
%   detect_percent_i: percent of detected number of auto inside contours compared to total number of  manual inside contours
%   detect_percent_o: percent of detected number of auto outside contours compared to total number of  manual outside contours
%   auto_missing_index_i: missing auto inside contours' index
%   auto_missing_index_o: missing auto outside contours' index
%   auto_bad_index_i: bad auto inside contours' index (bad contous is the contour whose average distance larger than dist_limit)
%   auto_bad_index_o: bad auto outside contours' index 
%   good_percent_i:  percent of good auto inside contours (good contous is the detected contour whose average distance smaller than dist_limit)
%   good_percent_o:percent of good auto outside contours
%   auto_ef: auto contour's ejection fraction
%   manual_ef: manual contour's ejection fraction
%   ef_diff:  auto_ef - manual_ef
%   auto_lvm: auto contour's lv mass
%   manual_lvm: manual contour's lv mass
%   lvm_diff:  auto_lvm - manual_lvm
%   avg_dist_i: average distance of inside contours
%   avg_dist_o: average distance of outside contours
%   avg_dm_i: average dice metric of inside contours
%   avg_dm_o: average dice metric of outside contours
%-----------------------------------------------------

%-initialize
compare_result = [];

%-parameters
para.name_prefix = 'IM-0001-'; %according to the IM-0001-XXXX naming standard
para.digit_length=4; %IM-0001-XXXX, XXXX is 4 digits
%contour mode
para.inside_contour_mode = 'icontour';
para.outside_contour_mode = 'ocontour';
%segmentation mode
para.auto_seg_mode='auto';
para.manual_seg_mode ='manual';
%for average distance
init_value = -999;
dist_limit = 4; %mm. If a contour's average distance larger than this number, it will be considered as bad contour

%-check path
if ~exist(dicom_path,'dir')
    disp([dicom_path ' :NOT exist!'])
    return;
end
if ~exist(manual_contour_path,'dir')
    disp([manual_contour_path ' :NOT exist!'])
    return;
end
if ~exist(auto_contour_path,'dir')
    disp([auto_contour_path ' :NOT exist!'])
    return;
end

%-check files
dicom_files = dir([dicom_path filesep '*.dcm']);
if isempty(dicom_files)
    disp([dicom_path ' :NO dicom files!'])
    return;
end
manual_contour_files = dir([manual_contour_path filesep '*.txt']);
if isempty(manual_contour_files)
    disp([manual_contour_path ' :NO manual contours!'])
    return;
end
auto_contour_files = dir([auto_contour_path filesep '*.txt']);
if isempty(auto_contour_files)
    disp([auto_contour_path ' :NO auto contours!'])
    return;
end

%-dicominfo
try
  dicom_filename = dicom_files(1).name; %use the first dicom file.
  full_dicom_filename = [dicom_path filesep dicom_filename];
  dicom_meta= dicominfo(full_dicom_filename);
  
  para.width = dicom_meta.Width;%image width
  para.height = dicom_meta.Height;%image height
  para.pixel_spacing = dicom_meta.PixelSpacing; % mm
  para.thickness = dicom_meta.SliceThickness;% mm
  para.gap = dicom_meta.SpacingBetweenSlices - para.thickness; % mm
  para.phase_number = dicom_meta.CardiacNumberOfImages; %numer of phases
  para.image_number = length(dicom_files);%number of images
catch
  s = lasterror;
  disp(s.message);
  return
end

%-calc percent of auto contours compared to manual contours
manual_number_i = 0; %number of manual inside contours
manual_number_o = 0; %number of manual outside contours
manual_index_i = []; %index of manual inside contours
manual_index_o = []; %index of manual outside contours

%count manual contours
for ix = 1:para.image_number 
    sindex = add_zero_index(ix,para.digit_length);
    %inside contour
    full_icontour_filename = get_contour_filename(manual_contour_path, para.name_prefix, sindex, para.inside_contour_mode,para.manual_seg_mode);
    if exist(full_icontour_filename,'file') 
       manual_number_i = manual_number_i + 1;
       manual_index_i = [manual_index_i; ix]; %#ok<AGROW>
    end
    %outside contour
    full_ocontour_filename = get_contour_filename(manual_contour_path,para.name_prefix, sindex, para.outside_contour_mode,para.manual_seg_mode);
    if exist(full_ocontour_filename,'file') 
       manual_number_o = manual_number_o + 1;
       manual_index_o = [manual_index_o; ix]; %#ok<AGROW>
    end
end
if (manual_number_i == 0)
   disp([manual_contour_path ' :NO inside contours!'])
   return; 
end
if (manual_number_o == 0)
   disp([manual_contour_path ' :NO outside contours!'])
   return; 
end

%count auto inside contours according to manual_index_i 
auto_number_i = 0;
auto_missing_index_i = [];
for ix = 1:length(manual_index_i)
    sindex = add_zero_index(manual_index_i(ix),para.digit_length);
    full_icontour_filename = get_contour_filename(auto_contour_path,para.name_prefix, sindex, para.inside_contour_mode, para.auto_seg_mode);
    if exist(full_icontour_filename,'file') 
       auto_number_i = auto_number_i + 1;
    else
       auto_missing_index_i = [auto_missing_index_i  manual_index_i(ix)]; %#ok<AGROW>
    end
end
if (auto_number_i == 0)
   disp([auto_contour_path  ' :NO inside contours!'])
   return; 
end

%count auto outside contours according to manual_index_o 
auto_number_o = 0;
auto_missing_index_o = [];
for ix = 1:length(manual_index_o)
    sindex = add_zero_index(manual_index_o(ix),para.digit_length);
    full_ocontour_filename = get_contour_filename(auto_contour_path,para.name_prefix, sindex, para.outside_contour_mode,para.auto_seg_mode);
    if exist(full_ocontour_filename,'file') 
       auto_number_o = auto_number_o + 1;
    else
       auto_missing_index_o = [auto_missing_index_o  manual_index_o(ix)]; %#ok<AGROW>
    end
end
if (auto_number_o == 0)
   disp([auto_contour_path ' :NO outside contours!'])
   return; 
end

%-calc perpendicular distance & Dice Metric
%inside contours
avg_dist_i=ones(length(manual_index_i),1)*init_value;
dm_i = ones(length(manual_index_i),1)*init_value;
auto_bad_index_i = [];

for ix = 1:length(manual_index_i)
    sindex = add_zero_index(manual_index_i(ix),para.digit_length);
    % manual contour
    manual_icontour_filename = get_contour_filename(manual_contour_path,para.name_prefix,sindex, para.inside_contour_mode,para.manual_seg_mode);
    % auto contour
    auto_icontour_filename = get_contour_filename(auto_contour_path,para.name_prefix,sindex, para.inside_contour_mode,para.auto_seg_mode);
    % if both exist, compare
    if (exist(manual_icontour_filename,'file') &&  exist(auto_icontour_filename,'file'))
         manual_xy =  load(manual_icontour_filename);
         auto_xy =  load(auto_icontour_filename);
         try
            %perpendicular distance
            avg_dist_i(ix) =calc_dist(auto_xy(:,1),auto_xy(:,2),manual_xy(:,1),manual_xy(:,2),para.pixel_spacing(1));
            %dice metric
            dm_i(ix) = calc_dm(auto_xy,manual_xy,para);
         catch
            s = lasterror;
            disp(s.message);
            continue;
         end
     end
end
if max(avg_dist_i) == init_value
    return;
end
auto_bad_ix_i =find(avg_dist_i >= dist_limit);
if ~isempty(auto_bad_ix_i)
    auto_bad_index_i = manual_index_i(auto_bad_ix_i);
end
auto_good_ix_i =find(avg_dist_i < dist_limit & avg_dist_i ~= init_value);


%outside contours
auto_bad_index_o = [];
avg_dist_o=ones(length(manual_index_o),1)*init_value;
dm_o = ones(length(manual_index_o),1)*init_value;

for ix = 1:length(manual_index_o)
    sindex = add_zero_index(manual_index_o(ix),para.digit_length);
    % manual contour
    manual_ocontour_filename = get_contour_filename(manual_contour_path,para.name_prefix,sindex, para.outside_contour_mode,para.manual_seg_mode);
    % auto contour
    auto_ocontour_filename = get_contour_filename(auto_contour_path,para.name_prefix,sindex, para.outside_contour_mode,para.auto_seg_mode);
    % if both exist, compare
    if (exist(manual_ocontour_filename,'file') &&  exist(auto_ocontour_filename,'file'))
         manual_xy =  load(manual_ocontour_filename);
         auto_xy =  load(auto_ocontour_filename);
         try
            %perpendicular distance
            avg_dist_o(ix)  = calc_dist(auto_xy(:,1),auto_xy(:,2),manual_xy(:,1),manual_xy(:,2),para.pixel_spacing(1));
            %dice metric
            dm_o(ix) = calc_dm(auto_xy,manual_xy,para);
         catch
            s = lasterror;
            disp(s.message);
            continue;
         end
     end
end

if max(avg_dist_o) == init_value
    return;
end
auto_bad_ix_o =find(avg_dist_o >= dist_limit);
if ~isempty(auto_bad_ix_o)
    auto_bad_index_o = manual_index_o(auto_bad_ix_o);
end
auto_good_ix_o =find(avg_dist_o < dist_limit & avg_dist_o ~= init_value);

%-calc ejection fraction and lv mass
lv_manual = calc_clinical_para(dicom_path, manual_contour_path, 'manual',para);
lv_auto= calc_clinical_para(dicom_path, auto_contour_path, 'auto', para,1,para.image_number/para.phase_number,lv_manual.es,lv_manual.ed);

%record result
compare_result.auto_number_i = auto_number_i; %total auto inside contours
compare_result.auto_number_o = auto_number_o;
compare_result.manual_number_i = manual_number_i;
compare_result.manual_number_o = manual_number_o;
compare_result.detect_percent_i = auto_number_i / manual_number_i * 100; %percent of detected auto inside contours
compare_result.detect_percent_o = auto_number_o / manual_number_o * 100;
compare_result.auto_missing_index_i = auto_missing_index_i; % missing auto inside contours
compare_result.auto_missing_index_o = auto_missing_index_o;

compare_result.auto_bad_index_i = auto_bad_index_i; %bad auto inside contours' index (average distance larger than dist_limit)
compare_result.auto_bad_index_o = auto_bad_index_o;
compare_result.good_percent_i = (auto_number_i - length(auto_bad_index_i)) /manual_number_i * 100; %percent of good auto inside contours
compare_result.good_percent_o = (auto_number_o - length(auto_bad_index_o)) /manual_number_o * 100;

compare_result.auto_ef = lv_auto.ef;%auto's ef 
compare_result.manual_ef = lv_manual.ef; %manual's ef
compare_result.ef_diff = lv_auto.ef - lv_manual.ef; %auto's ef - manual's ef

compare_result.auto_lvm = lv_auto.lvm; %auto's lv mass
compare_result.manual_lvm = lv_manual.lvm; %manual's lv mass
compare_result.lvm_diff = lv_auto.lvm - lv_manual.lvm; %auto's lvm - manual's lvm

compare_result.avg_dist_i = mean(avg_dist_i(auto_good_ix_i)); %average distance
compare_result.avg_dist_o = mean(avg_dist_o(auto_good_ix_o));
compare_result.avg_dm_i = mean(dm_i(auto_good_ix_i)); %average dice metric
compare_result.avg_dm_o = mean(dm_o(auto_good_ix_o));

end %compare contour

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sub functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sindex = add_zero_index(index,digit_length)
    %add zero before index with total length of digit_length, for instance, change 20 to 0020.
    sindex=int2str(index);
    while length(sindex) < digit_length
           sindex = ['0', sindex]; %#ok<AGROW>
    end
end 

function full_contour_filename = get_contour_filename(contour_path, name_prefix, sindex, contour_mode, seg_mode)
    %gef contour filename with full path
    full_contour_filename = [contour_path filesep name_prefix, sindex, '-',contour_mode,'-',seg_mode, '.txt'];
end

function varargout = calc_dist(xPtsAuto,yPtsAuto,xPtsExp,yPtsExp,pixel_spacing)    
    %calc average perpendicular distance

    % initialize
    mappingArray = zeros(length(xPtsAuto), 4);
    distance = zeros(length(xPtsAuto),1);

    for idx=1:length(xPtsAuto)  % ---- BEGIN FOR LOOP
        % moving to the next auto contour point to calculate
        if (idx == 1)     
            % initialize values (first point)
            currXPtAuto = xPtsAuto(1);
            currYPtAuto = yPtsAuto(1);

            prevXPtAuto = xPtsAuto(end);
            prevYPtAuto = yPtsAuto(end);

            nextXPtAuto = xPtsAuto(2);
            nextYPtAuto = yPtsAuto(2);
            
        elseif (idx == length(xPtsAuto))
            % last iteration (last point)
            prevXPtAuto = currXPtAuto;
            prevYPtAuto = currYPtAuto;

            currXPtAuto = nextXPtAuto;
            currYPtAuto = nextYPtAuto;

            nextXPtAuto = xPtsAuto(1);
            nextYPtAuto = yPtsAuto(1);
        else    
            % and all other cases
            prevXPtAuto = currXPtAuto;
            prevYPtAuto = currYPtAuto;

            currXPtAuto = nextXPtAuto;    
            currYPtAuto = nextYPtAuto;

            nextXPtAuto = xPtsAuto(idx+1);
            nextYPtAuto = yPtsAuto(idx+1);
        end        

        if (prevXPtAuto == nextXPtAuto) 
            m1 = 0;
            b1 = currYPtAuto;
            normal = 'y=currYPtAuto';
            y = currYPtAuto;
        elseif (prevYPtAuto == nextYPtAuto)
           m1 = inf;
           b1 = inf;
           normal = 'x=currXPtAuto';
           x = currXPtAuto;
        else
           p = csapi([prevXPtAuto; currXPtAuto; nextXPtAuto],[prevYPtAuto; currYPtAuto; nextYPtAuto]);
           tangent = fnder(p);
           m1 = -1/ppval(tangent,currXPtAuto);
           b1 = currYPtAuto - m1*currXPtAuto;
           normal = 'y=m1*x+b1';
        end

        capToExpPtsD = sqrt((xPtsExp - currXPtAuto).^2 + (yPtsExp - currYPtAuto).^2); 
        tolerance = max(capToExpPtsD/2);
        expPtsIdx = find(capToExpPtsD < tolerance);        
        xExpSubset = xPtsExp(expPtsIdx);
        yExpSubset = yPtsExp(expPtsIdx);

        if (strcmp(normal,'x=currXPtAuto') == true)
            expPtsToNormalD = abs(xExpSubset - currXPtAuto);
        else
            expPtsToNormalD = abs( (m1*xExpSubset - yExpSubset + b1))/sqrt(m1^2 + (-1)^2);
        end
        
        expSubset = [xExpSubset yExpSubset expPtsIdx expPtsToNormalD capToExpPtsD(expPtsIdx)];
        
        [minDToNormal] = min(expSubset(:,4));
        expSubset = expSubset(expSubset(:,4) == minDToNormal,:);

        [minDToPt, nearestExpPtIdx] = min(expSubset(:,5));
        xPtNearestExp = expSubset(nearestExpPtIdx,1);
        yPtNearestExp = expSubset(nearestExpPtIdx,2);
        nearestExpPtIdx = expSubset(nearestExpPtIdx,3);

        xPtExp = xPtNearestExp;
        yPtExp = yPtNearestExp;    

        if (nearestExpPtIdx==1)
            prevXPtExp = xPtsExp(end);
            prevYPtExp = yPtsExp(end);
            nextXPtExp = xPtsExp(2); 
            nextYPtExp = yPtsExp(2);
        elseif (nearestExpPtIdx==length(xPtsExp))
            prevXPtExp = xPtsExp(nearestExpPtIdx-1);
            prevYPtExp = yPtsExp(nearestExpPtIdx-1);
            nextXPtExp = xPtsExp(1);
            nextYPtExp = yPtsExp(1);
        else
             prevXPtExp = xPtsExp(nearestExpPtIdx-1);
             prevYPtExp = yPtsExp(nearestExpPtIdx-1);
             nextXPtExp = xPtsExp(nearestExpPtIdx+1);
             nextYPtExp = yPtsExp(nearestExpPtIdx+1);
        end
        if (prevXPtExp == nextXPtExp)
            x = xPtExp;

            if (~exist('y','var'))
                y = eval(solve(normal,'y'));    
            end
        elseif (prevYPtExp == nextYPtExp)
            y = yPtExp;
            if (~exist('x','var'))
                x = eval(solve(normal,'x'));
            end 
        else   
            LnExpParams = polyfit([prevXPtExp xPtExp nextXPtExp],[prevYPtExp yPtExp nextYPtExp],1);
            m2 = LnExpParams(1); %#ok<NASGU>
            b2 = LnExpParams(2); %#ok<NASGU>
            LnExp = 'y=m2*x+b2';
            [x,y] = solve(normal,LnExp);
            x = eval(x);
            y = eval(y);
        end

        distance(idx) = sqrt((currXPtAuto - x)^2 + (currYPtAuto - y)^2);
             
        clear x y normal m1 m2 b1 b2;

        mappingArray(idx, 1) = currXPtAuto;
        mappingArray(idx, 2) = currYPtAuto;
        mappingArray(idx, 3) = xPtExp;
        mappingArray(idx, 4) = yPtExp;    
    end   % ---- END FOR LOOP

    avgD_px = mean(distance);
    maxD_px = max(distance); 
    minD_px = min(distance);
    
    avgD_mm = avgD_px*pixel_spacing;
    maxD_mm = maxD_px*pixel_spacing;
    minD_mm = minD_px*pixel_spacing;
        
    varargout{1} = avgD_mm;
    varargout{2} = maxD_mm;
    varargout{3} = minD_mm;
    varargout{4} = distance.*pixel_spacing;
    varargout{5} = mappingArray;
end

function dm = calc_dm(autoPoints,manualPoints,para)
    %calc dice metric
    auto_mask = poly2mask (autoPoints(:,1),autoPoints(:,2),double(para.width),double(para.height));
    manual_mask = poly2mask (manualPoints(:,1),manualPoints(:,2),double(para.width),double(para.height));

    auto_size = sum(auto_mask(:)>0);
    manual_size = sum(manual_mask(:)>0);
    intersect_size = sum((auto_mask(:) + manual_mask(:))==2);
    dm = 2 * intersect_size / (auto_size + manual_size);
end

function lv = calc_clinical_para(dicomPath, contour_path, autoManual, para,varargin)
    %calc ejection fraction and lv mass
       
    [dicomCount, sliceCount] = dicom_counter(dicomPath,para.phase_number);
    contTable_i = zeros(sliceCount, para.phase_number);
    contTable_o = zeros(sliceCount, para.phase_number);
    
    if length(varargin) > 3
        startingSlice = varargin{1};
        endingSlice = varargin{2};
        systolePhase = varargin{3};
        diastolePhase = varargin{4};
        
    elseif length(varargin) > 1
        startingSlice = varargin{1};
        endingSlice = varargin{2};
        systolePhase = 0;
        diastolePhase = 0;
        
    else 
        startingSlice = 1;
        endingSlice = sliceCount;
        systolePhase = 0;
        diastolePhase = 0;    
    end    
    
    if startingSlice < 1 || startingSlice > sliceCount
        startingSlice = 1;
    end

    if endingSlice < 1 || endingSlice > sliceCount
        endingSlice = sliceCount;
    end
    
    if startingSlice > endingSlice
        endTemp = startingSlice;
        startingSlice = endingSlice;
        endingSlice = endTemp;
    end   
            
    contourList = dir(contour_path);
    for contourIdx = 1:length(contourList)
        contourFile = contourList(contourIdx).name;
        suffixTag = regexp(contourFile, ['contour-', autoManual, '.txt'], 'start');
        
        if isempty(suffixTag)
            continue
        end
        
        contourType = contourFile(suffixTag(1) - 1); % (i or o)
        if strcmpi(contourType, 'i')  %inside contour   
            match = regexp((contourFile),'-','start');
            imNum = contourFile(match(2) + 1:match(3) - 1);
            currSlice = getSlice(str2double(imNum),para.phase_number);
            currPhase = getPhase(str2double(imNum),para.phase_number);
            full_contour_filename = [contour_path filesep contourFile];
            
            try
               xy = load(full_contour_filename);
               area = polyarea(xy(:,1),xy(:,2));
               area_mm = area* para.pixel_spacing(1)^2;
               vol_cm3 = area_mm * (para.thickness + para.gap) *(1/10)^3; %change to cm3
               contTable_i(currSlice, currPhase) = vol_cm3; 
           catch
               s = lasterror;
               disp(s.message);
               return;
            end
            
        end %end if
        
        if strcmpi(contourType, 'o')  %outside contour
            match = regexp((contourFile),'-','start');
            imNum = contourFile(match(2) + 1:match(3) - 1);
            currSlice = getSlice(str2double(imNum),para.phase_number);
            currPhase = getPhase(str2double(imNum),para.phase_number);
            full_contour_filename = [contour_path filesep contourFile];
            try
               xy = load(full_contour_filename);
               area = polyarea(xy(:,1),xy(:,2));
               area_mm = area* para.pixel_spacing(1)^2;
               vol_cm3 = area_mm * (para.thickness + para.gap) *(1/10)^3;%change to cm3
               contTable_o(currSlice, currPhase) = vol_cm3; 
           catch
               s = lasterror;
               disp(s.message);
               return;
            end
        end %end if
    end     
  
    maxPhase = length(contTable_i(1,:));
    if systolePhase < 1 || systolePhase > maxPhase || diastolePhase < 1 || diastolePhase > maxPhase
        [systolePhase, diastolePhase] = esedDetermine(contTable_i);
    end
    
    [contTableChecked_i, zeroed, slicesExcluded] = efCheck(contTable_i, systolePhase, diastolePhase, startingSlice, endingSlice);

    %calc ef
    esvol = sum(contTableChecked_i(:, systolePhase));
    edvol = sum(contTableChecked_i(:, diastolePhase));
    strokeVol = edvol - esvol;
    ef = strokeVol / edvol * 100;
    
    lv.esv = esvol;
    lv.edv = edvol;
    lv.sv = strokeVol;
    lv.ef = ef;
    lv.es = systolePhase;
    lv.ed = diastolePhase;
    lv.zeroed = zeroed;
    lv.slice_excluded = slicesExcluded;
   
    %calc lv mass
    combined_table(:,1) = contTable_i(:,diastolePhase);
    combined_table(:,2) = contTable_o(:,diastolePhase);
    common_slice = ((combined_table(:,1) ~= 0) + (combined_table(:,2) ~= 0)) == 2 ;
    
    edv_i = sum(combined_table(common_slice, 1));
    edv_o = sum(combined_table(common_slice, 2));
    lv.lvm = (edv_o - edv_i) * 1.05; %1.05 (g/cm3)
end

function [systolePhase, diastolePhase] = esedDetermine(contTable)
    %determin ES and ED phase 

    contVolSum = sum(contTable);
    minVol = 9999; 
    maxVol = 0;
    
    for x = 1:length(contVolSum)
        valueUnderExamination = contVolSum(x);
        
        if valueUnderExamination > 0
            if minVol > valueUnderExamination
                minVol = contVolSum(x);
                systolePhase = x;
            end

            if maxVol < valueUnderExamination
                maxVol = contVolSum(x);
                diastolePhase = x;
            end
        end
    end
end

function sliceNum = getSlice(imageNum,phase_number)
    % calc slice number from the image number
     sliceNum = ceil(imageNum/phase_number);
end

function phaseNum = getPhase(imageNum,phase_number)
    %calc phase number from image number
    remainder = mod(imageNum,phase_number);
    if (remainder ~= 0)
        phaseNum = remainder;
    else
        phaseNum = phase_number;
    end 
end

function [volTableCropped, zeroed, slicesExcluded] = efCheck(volTable, systolePhase, diastolePhase, startSlice, endSlice)
    % modifies the contour table that is passed in by...
    % 1. removing the slices that are to be excluded using startSlice and
    % endSlice by zeroing all those values
    % 2. removing ED/ES pairs that are missing the other value.  
    % that is, if there is a zero in the systolic phase but not in the
    % diastolic phase in that slice, zero the other
    
    % 1. removing the non-included slices
    volTableCropped = zeros(length(volTable(:,1)), length(volTable(1,:)));
    volTableCropped(startSlice:endSlice, :) = volTable(startSlice:endSlice, :);
    
    % 2. finding the ES/ED slices with missing values
    sysZeros = find(volTable(:, systolePhase) == 0);
    diaZeros = find(volTable(:, diastolePhase) == 0);
    
    if ~isempty(sysZeros)
        volTableCropped(sysZeros, :) = 0;
    end
    
    if ~isempty(diaZeros)
        volTableCropped(diaZeros, :) = 0;
    end
    
    slicesExcluded = find(volTableCropped(:, systolePhase) == 0);
    zeroed = length(slicesExcluded);
end

function [dicom_count, slice_count] = dicom_counter(dicom_path, phase_number)
    %calc dicom image number and slice number
    dicom_count = 0;
    slice_count = 0;
    
    if ~exist(dicom_path,'dir')
        return
    end
    
    dicom_files = dir([dicom_path filesep '*.dcm']);
    if isempty(dicom_files)
        disp([dicom_path ' :NO DICOM files!'])
        return;
    else 
        dicom_count = length(dicom_files);
        slice_count = ceil(dicom_count / phase_number);
    end
end
