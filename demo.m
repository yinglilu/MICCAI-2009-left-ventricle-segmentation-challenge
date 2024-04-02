%demo.m
%   DEMO Demonstration of how to use compare_contours.m and write_results.m
%
%   Copyright: Imaging Research, Sunnybrook Health Sciences Centre, Toronto, ON, Canada
%   Author: Perry Radau, perry.radau@gmail.com
%   Version: 0.9   
%   Date: 2009/05/11

clear

%single mode:compare one patient's manual and auto contours
%compare_result = compare_contours('C:\segmentation_projects\all\SC-HF-I-10\DICOM','C:\segmentation_projects\all\SC-HF-I-10\contours-manual\2009-04','C:\segmentation_projects\all\SC-HF-I-10\contours-auto\yingli');

%batch mode:compare multi-patient's mauanl and auto contours
folder_prefix = 'SC*';
manual_contour_folder_name = 'contours-manual';
auto_contour_folder_name = 'contours-auto';
root_folder = 'e:\segmentation_projects\all';
dirs = dir([root_folder filesep folder_prefix]);

%check path
if isempty(dirs)
   disp( [root_folder filesep folder_prefix ' :NOT found!'])
   return;
end

%initialize
all_compare_result = cell(length(dirs),1);

for i = 1:2%length(dirs)
    display(['Processing: ' dirs(i).name]);
    %full dicom folder
    dicom_folder = [root_folder filesep dirs(i).name filesep 'DICOM'];
    %full manual contour folder
    manual_contour_folder = [root_folder filesep dirs(i).name filesep manual_contour_folder_name filesep '2009-04'];
    %full auto contour folder
    auto_contour_folder = [root_folder filesep dirs(i).name filesep auto_contour_folder_name filesep 'yingli'];
    %compare
    compare_result = compare_contours(dicom_folder, manual_contour_folder, auto_contour_folder);
    %record results
    if ~isempty(compare_result)
        all_compare_result{i}.patient_id = dirs(i).name;
        all_compare_result{i}.auto_number_i = compare_result.auto_number_i;
        all_compare_result{i}.auto_number_o= compare_result.auto_number_o;
        all_compare_result{i}.manual_number_i = compare_result.manual_number_i;
        all_compare_result{i}.manual_number_o = compare_result.manual_number_o;
        all_compare_result{i}.detect_percent_i = compare_result.detect_percent_i;
        all_compare_result{i}.detect_percent_o = compare_result.detect_percent_o;
        all_compare_result{i}.good_percent_i = compare_result.good_percent_i;
        all_compare_result{i}.good_percent_o = compare_result.good_percent_o;
        
        all_compare_result{i}.auto_ef = compare_result.auto_ef;
        all_compare_result{i}.manual_ef = compare_result.manual_ef;
        all_compare_result{i}.ef_diff= compare_result.ef_diff;
        
        all_compare_result{i}.auto_lvm= compare_result.auto_lvm;
        all_compare_result{i}.manual_lvm = compare_result.manual_lvm;
        all_compare_result{i}.lvm_diff = compare_result.lvm_diff;
        
        all_compare_result{i}.avg_dist_i = compare_result.avg_dist_i;
        all_compare_result{i}.avg_dist_o = compare_result.avg_dist_o;
        all_compare_result{i}.avg_dm_i = compare_result.avg_dm_i;
        all_compare_result{i}.avg_dm_o = compare_result.avg_dm_o;
   end
end

%write result so xls(excel) file
write_results('e:\all_compare_results', all_compare_result);

disp('Done!');
