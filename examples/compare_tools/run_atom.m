% Types one structure (IN.txt from dump_clay.py) with clayff_atom of the atom Toolbox (Holmboe, Clays Clay Miner. 67,
% 419 (2019); https://github.com/mholmboe/atom), with 'clayff' and SPC water (Cygan's ClayFF; atom's default is SPC/E), no healing (default),
% metal-O cutoff 2.25 A.  Writes OUT.txt (index, type, charge per atom) and OUT.neigh (atom's neighbour strings).
% usage: octave-cli run_atom.m ATOM_DIR IN.txt OUT.txt      (MATLAB: set args = {ATOM_DIR, IN, OUT} and run)
if exist('argv', 'builtin') || exist('argv', 'file'), args = argv(); end
atomdir = args{1}; infile = args{2}; outfile = args{3};
addpath(genpath(atomdir));
cd(atomdir);  % bond_atom loads general_functions/Revised_Shannon_radii.mat by a relative path
fid = fopen(infile, 'r');
Box_dim = str2double(strsplit(strtrim(fgetl(fid))));
C = textscan(fid, '%f %s %s %f %f %f');
fclose(fid);
atom = struct('molid', {}, 'resname', {}, 'type', {}, 'fftype', {}, 'index', [], ...
    'neigh', struct('type', {{}}, 'index', zeros(6, 1), 'dist', zeros(6, 1)), ...
    'bond', struct('type', zeros(6, 1), 'index', zeros(6, 1)), ...
    'angle', struct('type', zeros(6, 1), 'index', zeros(6, 1)), ...
    'x', [], 'y', [], 'z', [], 'vx', [], 'vy', [], 'vz', []);
for i = 1:numel(C{1})  % the fields as import_atom_gro sets them
    atom(i).molid = C{1}(i);
    atom(i).resname = C{2}(i);
    atom(i).type = C{3}(i);
    atom(i).fftype = C{3}(i);
    atom(i).index = i;
    atom(i).neigh.index = zeros(6, 1); atom(i).neigh.dist = zeros(6, 1);
    atom(i).bond.type = zeros(6, 1); atom(i).bond.index = zeros(6, 1);
    atom(i).angle.type = zeros(6, 1); atom(i).angle.index = zeros(6, 1);
    atom(i).x = C{4}(i); atom(i).y = C{5}(i); atom(i).z = C{6}(i);
    atom(i).vx = NaN; atom(i).vy = NaN; atom(i).vz = NaN;
end
t0 = tic;
atom = clayff_atom(atom, Box_dim, 'clayff', 'spc');
dt = toc(t0);
fid = fopen(outfile, 'w');
fprintf(fid, '# seconds %.3f natoms %d\n', dt, numel(atom));
for i = 1:numel(atom)
    fprintf(fid, '%d %s %.6f\n', atom(i).index, char(atom(i).type), atom(i).charge);
end
fclose(fid);
fid = fopen([outfile(1:end-4) '.neigh'], 'w');
for k = 1:size(All_Neighbours, 1)
    r = All_Neighbours(k, :);
    fprintf(fid, '%s %s %s %d\n', char(r{2}), char(r{3}), char(r{4}), r{5});
end
fclose(fid);
