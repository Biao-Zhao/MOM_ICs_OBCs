function write_uv_to_newini(old_ini, new_ini, u_new, v_new, it)
% ------------------------------------------------------------
% write_uv_to_newini
%
% Copy an existing NetCDF initialization file and overwrite
% the variables 'u' and 'v' with newly computed fields.
%
% Parameters:
%   old_ini : string   Path to the original NetCDF file
%   new_ini : string   Path for the new NetCDF file (output)
%   u_new   : array    New u-field  [z, y_u, x_u]
%   v_new   : array    New v-field  [z, y_v, x_v]
%   it      : integer  Time index to write (1-based)
%
% Notes:
%   - NetCDF uses zero-based indexing, so actual write index = it - 1
%   - u_new and v_new are automatically converted to single precision
%   - All dimensions, attributes, and metadata from old_ini are preserved
% ------------------------------------------------------------

    arguments
        old_ini (1,1) string
        new_ini (1,1) string
        u_new (:,:,:) {mustBeNumeric}
        v_new (:,:,:) {mustBeNumeric}
        it (1,1) double {mustBeInteger, mustBePositive}
    end

    % === 1. Check and copy file ===
    if ~isfile(old_ini)
        error('L File does not exist: %s', old_ini);
    end

    copyfile(old_ini, new_ini);
    fprintf('Copied %s %s\n', old_ini, new_ini);

    % === 2. Open the new file for writing ===
    ncid = netcdf.open(new_ini, 'NC_WRITE');

    % Get variable IDs for 'u' and 'v'
    try
       varid_u = netcdf.inqVarID(ncid, 'u');
       [name_u, ~, dimids_u, ~] = netcdf.inqVar(ncid, varid_u);
       for k = 1:numel(dimids_u)
           [dname, dlen] = netcdf.inqDim(ncid, dimids_u(k));
           fprintf('u dimension %d: %s (length=%d)\n', k, dname, dlen);
       end

       varid_v = netcdf.inqVarID(ncid, 'v');
       [name_v, ~, dimids_v, ~] = netcdf.inqVar(ncid, varid_v);
       for k = 1:numel(dimids_v)
           [dname, dlen] = netcdf.inqDim(ncid, dimids_v(k));
           fprintf('v dimension %d: %s (length=%d)\n', k, dname, dlen);
       end
    catch
        netcdf.close(ncid);
        error('L Variables "u" or "v" not found in file.');
    end
    
    % === 3. Write data for time index 'it' ===
    start = [0 0 0 it-1];  % [time, z, y, x] start indices
    count_u = [size(u_new,1) size(u_new,2) size(u_new,3) 1];
    count_v = [size(v_new,1) size(v_new,2) size(v_new,3) 1];

    % Write directly (no redefine mode needed for existing vars)
    fprintf('Writing time slice %d ...\n', it);
    netcdf.putVar(ncid, varid_u, start, count_u, double(u_new));
    netcdf.putVar(ncid, varid_v, start, count_v, double(v_new));

    % === 4. Close file ===
    netcdf.close(ncid);
    fprintf('Successfully wrote u/v variables to %s\n', new_ini);
end
