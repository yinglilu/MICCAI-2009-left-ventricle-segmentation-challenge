function write_results(result_path, all_compare_result) %#ok<DEFNU>
%WRITE_RESULT Write compare result to xls(excel) file
%   WRITE_RESULT(RESULT_PATH,ALL_COMPARE_RESULT)
%
%   Copyright: Imaging Research, Sunnybrook Health Sciences Centre, Toronto, ON, Canada
%   Author: Perry Radau, perry.radau@gmail.com
%   Version: 0.9   
%   Date: 2009/05/11

%check path
if ~exist(result_path,'dir')
    try
       mkdir(result_path)
    catch
       s = lasterror;
       disp(s.message);
       return;
    end
end

%write file with name yyyymmddTHHMMSS.xls
full_filename = [result_path filesep datestr(now,'yyyymmddTHHMMSS') '.xls'];

try 
    result_file = fopen(full_filename,'w');
    summary.headings = {'patient_id',...
                        'auto_number_i', 'auto_number_o', ...
                        'manual_number_i','manual_number_o',...
                        'detect_percent_i(%)', 'detect_percent_o(%)',...
                        'good_percent_i(%)','good_percent_o(%)',...
                        'auto_ef(%)','manual_ef(%)','ef_diff(%)',...
                        'auto_lvm(g)','manual_lvm(g)','lvm_diff(g)',...
                        'avg_dist_i(mm)','avg_dist_o(mm)',...
                        'avg_dm_i(mm)','avg_dm_o(mm)'
                     }; 
                 
    fprintf(result_file, '%s\t', summary.headings{:});

    for i=1:length(all_compare_result)
        if ~isempty(all_compare_result{i})
            fprintf(result_file, '\n %s\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t %6.4f',...
                        all_compare_result{i}.patient_id,...
                        all_compare_result{i}.auto_number_i, all_compare_result{i}.auto_number_o, ...
                        all_compare_result{i}.manual_number_i,all_compare_result{i}.manual_number_o,...
                        all_compare_result{i}.detect_percent_i, all_compare_result{i}.detect_percent_o,...
                        all_compare_result{i}.good_percent_i,all_compare_result{i}.good_percent_o,...
                        all_compare_result{i}.auto_ef,all_compare_result{i}.manual_ef,all_compare_result{i}.ef_diff,...
                        all_compare_result{i}.auto_lvm,all_compare_result{i}.manual_lvm,all_compare_result{i}.lvm_diff,...
                        all_compare_result{i}.avg_dist_i,all_compare_result{i}.avg_dist_o,...
                        all_compare_result{i}.avg_dm_i,all_compare_result{i}.avg_dm_o);
        end
    end
    fclose(result_file);
catch
    s = lasterror;
    disp(s.message);
    return;
end

