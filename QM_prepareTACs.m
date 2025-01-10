function [Cpet, roiTacs, num_rois] = QM_prepareTACs(handles,waitbarhandle,num_files)
% QM_prepareTACs function
% [Cpet roiTacs num_rois] = QM_prepareTACs(handles,waitbarhandle,num_files)
%
% Function which calcules a concentration vector of a PET study, TACs for
% the ROIs and the number of ROIs for the current PET study and mask. It is
% executed from QModeling GUI
%---------------------------------------------------------------------------------
% QM_prepareTACs is part of QModeling.
%
% QModeling is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% QModeling is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with QModeling.  If not, see <http://www.gnu.org/licenses/>.
%---------------------------------------------------------------------------------
% Copyright (C) 2013, 2020 Francisco Javier Lopez Gonzalez, Jose Paredes Pacheco, 
% Karl-Khader Thurnhofer Hemsi

mask_paths = get(handles.mask_Textbox,'UserData');
masks=0;

if strcmp(mask_paths{2},'ERROR') %file format is not correct
    ME = MException('PrepareTACs:fileFormatIncorrect', ...
                               'TACs file format not accepted');
    throw(ME);
elseif strcmp(mask_paths{2},'.txt')
    try
        TACs=importdata(mask_paths{1});
        if (length(TACs.data)~=num_files)
            ME = MException('PrepareTACs:dimensionsDismatch', ...
                               'Dimensions of the PET and mask dismatch');
            throw(ME);
        else
           num_rois=length(TACs.textdata);
           roiTacs = cell(2,num_rois);
           for i=1:num_rois
               roiTacs{1,i}=TACs.textdata{i};
               roiTacs{2,i}=TACs.data(:,i);
           end
               
            
        end
    catch ME
%         uiwait(errordlg(char({'An error ocurred reading TACs file';...
%                 'Please, read the manual to ensure your data is well formatted.'}),'Reading error','modal'));
        ME = MException('PrepareTACs:badFormattedData', ...
                               'Something wrong in data text');
        throw(ME);
    end
        
elseif (strcmp(mask_paths{2},'.xls') || strcmp(mask_paths{2},'.xlsx'))
    try
        [TACs, TACs_names, ~]=xlsread(mask_paths{1});
        if (length(TACs)~=num_files)
            ME = MException('PrepareTACs:dimensionsDismatch', ...
                               'Dimensions of the PET and mask dismatch');
            throw(ME);
        else
           num_rois=length(TACs_names);
           roiTacs = cell(2,num_rois);
           for i=1:num_rois
               roiTacs{1,i}=TACs_names{i};
               roiTacs{2,i}=TACs(:,i);
           end             
            
        end
    catch ME
%         uiwait(errordlg(char({'An error ocurred reading TACs file';...
%                 'Please, read the manual to ensure your data is well formatted.'}),'Reading error','modal'));
        ME = MException('PrepareTACs:badFormattedData', ...
                               'Something wrong in data text');
        throw(ME);
    end
    
elseif strcmp(mask_paths{2},'.mat')
    try
        aux=load(mask_paths{1});
        if (length(aux.TACs)~=num_files)
            ME = MException('PrepareTACs:dimensionsDismatch', ...
                               'Dimensions of the PET and mask dismatch');
            throw(ME);
        else
           num_rois=length(aux.TACs_names);
           roiTacs = cell(2,num_rois);
           for i=1:num_rois
               roiTacs{1,i}=aux.TACs_names{i};
               roiTacs{2,i}=aux.TACs(:,i);
           end             
            
        end
    catch ME
%         uiwait(errordlg(char({'An error ocurred reading TACs file';...
%                 'Please, read the manual to ensure your data is well formatted.'}),'Reading error','modal'));
        ME = MException('PrepareTACs:badFormattedData', ...
                               'Something wrong in data text');
        throw(ME);
    end
    
else
    masks=1;
    mask_hdr = spm_vol(mask_paths{1});
    mask = spm_read_vols(mask_hdr);

    id_rois = unique(mask);
    id_rois = id_rois(2:end);
    num_rois = length(id_rois);

    %Create a cell with mask path and results
    roiTacs = cell(2,num_rois);

    %Load ROI's names from txt
    exist_codes = false;
    if ~isempty(mask_paths{2})
        try
            fid = fopen(mask_paths{2});
            C = textscan(fid, '%s%d');
            fclose(fid);
            exist_codes = true;
        catch ME
            uiwait(errordlg(char({'An error ocurred reading ROIs name file';...
                'Default names for the ROIs have been loaded.'}),'Reading error','modal'));
    %         errordlg(char({'An error ocurred reading ROIs name file';...
    %             'Default names for the ROIs have been loaded.'}),'Reading error','modal');
        end
    else
        uiwait(warndlg('Default names for the ROIs have been loaded.','Reading error','modal'));
    %     warndlg('Default names for the ROIs have been loaded.','Reading error','modal')
    end


    %Inicialize the cell of results
    for i = 1:num_rois
        if exist_codes == true
            if isempty(C{1}(C{2}==id_rois(i)))
                roiTacs{1,i} = strcat('TAC ',num2str(i));
            else
                roiTacs{1,i} = char(C{1}(C{2}==id_rois(i)));
            end
        else
            roiTacs{1,i} = strcat('TAC ',num2str(i));
        end
        roiTacs{2,i} = zeros(num_files,1);      %Results of roi
    end
end

paths = get(handles.study_Textbox,'UserData');
num_files = size(paths,1);

if num_files == 1 

        %Open the multiframe PET study
    %     multifile = paths;
    multifile = paths(1,1:max(strfind(paths,','))-1);      % To raise ',1'
    multifile_hdr = spm_vol(multifile);
    num_files = length(multifile_hdr);

        %Are dimensions match?
    if (masks==1 && ~isequal(multifile_hdr(1).dim, mask_hdr.dim))
        ME = MException('PrepareTACs:dimensionsDismatch', ...
                               'Dimensions of the PET and mask dismatch');
        throw(ME);
    end

        %Variable to store all data 
    Cpet = zeros(num_files,prod(multifile_hdr(1).dim));

    for frame = 1:num_files

        current_im = spm_read_vols(multifile_hdr(frame));

            %Almaceno las tac por voxel
        Cpet(frame,:) = current_im(:);

        if masks==1     %Make the average for the voxel concentrations
            for i = 1:num_rois
                %Eliminate posible NaNs from mean calculation
                current_region=current_im(mask == id_rois(i));
                roiTacs{2,i}(frame) = mean(current_region(~isnan(current_region)));
                %roiTacs{2,i}(frame) = mean(current_im(mask == id_rois(i)));
            end
        end
        waitbar(frame/num_files,waitbarhandle)
    end

else

        %Are dimensions match?
    file = paths(1,:);
    file_hdr = spm_vol(file);
    if (masks==1 && ~isequal(file_hdr(1).dim, mask_hdr.dim))
        ME = MException('PrepareTACs:dimensionsDismatch', ...
                               'Dimensions of the PET and mask dismatch');
        throw(ME);
    end

        %Creo una variable donde almacenar todos los datos
    Cpet = zeros(num_files,prod(file_hdr(1).dim));

    for frame = 1:num_files

            %Open the current frame for the PET study
        file = paths(frame,:);
        file_hdr = spm_vol(file);
        current_im = spm_read_vols(file_hdr);   

            %Almaceno las tac por voxel
        Cpet(frame,:) = current_im(:);

        if masks==1    %Make the average for the voxel concentrations
            for i = 1:num_rois
                %Eliminate posible NaNs from mean calculation
                current_region=current_im(mask == id_rois(i));
                roiTacs{2,i}(frame) = mean(current_region(~isnan(current_region)));
                %roiTacs{2,i}(frame) = mean(current_im(mask == id_rois(i)));
            end
        end
        waitbar(frame/num_files,waitbarhandle)

    end
end


