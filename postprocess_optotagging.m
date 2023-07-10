nwb.processing.keys

%ot = nwb.processing.get('optotagging').nwbdatainterface.get('optotagging').data.load();

optocond = nwb.processing.get('optotagging').dynamictable.get('optogenetic_stimulation').vectordata.get('condition').data.load();
optodur = nwb.processing.get('optotagging').dynamictable.get('optogenetic_stimulation').vectordata.get('duration').data.load();
optolevel = nwb.processing.get('optotagging').dynamictable.get('optogenetic_stimulation').vectordata.get('level').data.load();
optostim = nwb.processing.get('optotagging').dynamictable.get('optogenetic_stimulation').vectordata.get('stimulus_name').data.load();

o = cellstr([optocond num2str(optodur)]);
v = unique(o);
c = zeros(size(v));
for ii = 1:numel(v)
c(ii) = nnz(strcmp(o,v(ii)));
fprintf('%s %d\n', v{ii},c(ii))
end

optostim = cellstr(optostim);
v = unique(optostim);
c = zeros(size(v));
for ii = 1:numel(v)
c(ii) = nnz(strcmp(optostim,v(ii)));
fprintf('%s %d\n', v{ii},c(ii))
end


optocond = cellstr(optocond);
v = unique(optocond);
c = zeros(size(v));
for ii = 1:numel(v)
c(ii) = nnz(strcmp(optocond,v(ii)));
fprintf('%s %d\n', v{ii},c(ii))
end


optostarttime = nwb.processing.get('optotagging').dynamictable.get('optogenetic_stimulation').start_time.data.load();
optostoptime = nwb.processing.get('optotagging').dynamictable.get('optogenetic_stimulation').stop_time.data.load();

optotimeseries = nwb.processing.get('optotagging').dynamictable.get('optogenetic_stimulation').timeseries.data.load();
optostoptime

