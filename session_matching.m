%% match sessions between this dataset and v240130
datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsesv240130 = {nwbdir.name};
nwbsesv240130 = nwbsesv240130(~contains(nwbsesv240130, 'Placeholder') & ...
    ( contains(nwbsesv240130, 'sub-') | contains(nwbsesv240130, 'sub_') ));
neulocaggv240130 = cell(size(nwbsesv240130));
for ises = 1:numel(nwbsesv240130)
    clearvars neuloc unit_peakch electrode_location electrode_id
    pathpp = ['S:\OpenScopeData\00248_v240130\postprocessed\' nwbsesv240130{ises} '\'];
    load([pathpp 'info_electrodes.mat']) %'electrode_probeid', 'electrode_localid', 'electrode_id', 'electrode_location', '-v7.3')
    load([pathpp 'info_units.mat']) %'unit_ids', 'unit_peakch', 'unit_times_idx', 'unit_wfdur') %'unit_times_data',
    elecid = electrode_id+1;
    revmapelecid = NaN(max(elecid),1);
    revmapelecid(elecid) = 1:numel(elecid);
    neuloc = electrode_location(revmapelecid(unit_peakch+1));
    neulocaggv240130{ises} = neuloc;
end

datadir = 'S:\OpenScopeData\000248\';
nwbdir = dir(datadir);
nwbsesv0 = {nwbdir.name};
nwbsesv0 = nwbsesv0(~contains(nwbsesv0, 'Placeholder') & ...
    ( contains(nwbsesv0, 'sub-') | contains(nwbsesv0, 'sub_') ));
neulocaggv0 = cell(size(nwbsesv0));
for ises = 1:numel(nwbsesv0)
    clearvars neuloc unit_peakch electrode_location electrode_id
    pathpp = ['S:\OpenScopeData\000248\postprocessed\' nwbsesv0{ises} '\'];
    load([pathpp 'info_electrodes.mat']) %'electrode_probeid', 'electrode_localid', 'electrode_id', 'electrode_location', '-v7.3')
    load([pathpp 'info_units.mat']) %'unit_ids', 'unit_peakch', 'unit_times_idx', 'unit_wfdur') %'unit_times_data',
    elecid = electrode_id+1;
    revmapelecid = NaN(max(elecid),1);
    revmapelecid(elecid) = 1:numel(elecid);
    electrode_location = cellstr(electrode_location);
    neuloc = electrode_location(revmapelecid(unit_peakch+1));
    neulocaggv0{ises} = neuloc;
end

sesindsv240130 = zeros(size(nwbsesv240130));
for ises = 1:numel(nwbsesv240130)
    c = cellfun(@isequal, neulocaggv0, repmat(neulocaggv240130(ises),size(neulocaggv0)), 'UniformOutput',false);
    sesv0ind = find(cat(1,c{:}));
    if numel(sesv0ind)==1
        sesindsv240130(ises) = sesv0ind;
    else
        warning('%d %s does not have a direct match with the old version', ises, nwbsesv240130{ises})
    end
end
% Warning: 1 sub-619293 does not have a direct match with the old version 
% Warning: 6 sub-625554 does not have a direct match with the old version 

isequal(neulocaggv240130{6}, neulocaggv0{4})

Nneuv0 = cellfun(@numel, neulocaggv0);
Nneuv240130 = cellfun(@numel, neulocaggv240130);
[r,c]= find(Nneuv0 == Nneuv240130');

isequal(neulocaggv240130{6}, neulocaggv0{1})
