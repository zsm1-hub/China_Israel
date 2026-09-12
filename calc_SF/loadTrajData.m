function traj = loadTrajData(fname, ii)
    chunk_fname = [fname(1:end-3), 'chunk', num2str(ii, '%02d'), 'traj.mat'];
    data = load(chunk_fname);
    traj = data.traj;
end
