function saveAsMALK(locs_ROI,xCol, yCol, framesCol, photonsCol, filename);

locs_subset = [];
locs_subset (:,1) = locs_ROI(:,xCol);
locs_subset (:,2) = locs_ROI(:,yCol);
locs_subset (:,3) = locs_ROI(:,framesCol);
locs_subset (:,4) = locs_ROI(:,photonsCol);

name_for_LAMA = [filename, '_MALK.txt'];

fid = fopen(name_for_LAMA,'wt');
fprintf(fid, '# localization file (Malk format) \n# x[nm] \t	y[nm]	\t [frame] \t	I[a.u.]');

dlmwrite(name_for_LAMA, locs_subset,'delimiter', '\t', '-append')
fclose(fid);

fprintf('\n -- Data saved in Malk format --\n')

end