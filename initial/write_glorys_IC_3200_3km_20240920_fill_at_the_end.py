#!/usr/bin/env python3
"""
script for preparing model IC (ssh,T,S,u,v) from Glorys data
How to use
./write_glorys_initial.py --config_file glorys_ic.yaml
./write_glorys_IC_3200_3km_20240920_fill_at_the_end.py  --config_file  glorys_ic_20240920_3200_3km_fill_at_the_end.yaml
"""

# author: 'Jing Chen'
# description: 'Initial conditions for MOM6, generated from GLORYS PHY fields'
# created: '2025-08-05'

import sys
import os
import argparse
import yaml

import numpy as np
import xarray
import xesmf

from HCtFlood import kara as flood

# Get the directory of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))

#
#sys.path.append(os.path.join(script_dir, './depths'))
from depths import vgrid_to_interfaces, vgrid_to_layers

#
sys.path.append(os.path.join(script_dir, '../boundary'))
from boundary import rotate_uv



def fill_from_deepest_valid(col):
    """
    Fill missing (NaN) values below the deepest valid point in each vertical column.
    col: 1D array (depth profile along 'zl')
    """
    if np.all(np.isnan(col)):
        return col  # no valid data at all
    last_valid_idx = np.where(~np.isnan(col))[0][-1]
    filled = col.copy()
    filled[last_valid_idx:] = col[last_valid_idx]
    return filled





def write_initial(config):
    # 1) Extract file paths from the new top-level YAML keys
    temp_file = config["glorys_temperature"]
    sal_file  = config["glorys_salinity"]
    ssh_file  = config["glorys_sea_surface_height"]
    u_file    = config["glorys_zonal_velocity"]
    v_file    = config["glorys_meridional_velocity"]


    vgrid_file = config['vgrid_file']
    grid_file = config['grid_file']
    output_file = config['output_file']
    reuse_weights = config.get('reuse_weights', False)

    # 2) Print them for debugging
    print("Reading from the following GLORYS files:")
    print(f"  Temperature: {temp_file}")
    print(f"  Salinity: {sal_file}")
    print(f"  SSH:       {ssh_file}")
    print(f"  U (zonal): {u_file}")
    print(f"  V (merid.):{v_file}")

    # 3) Retrieve variable names from the unchanged 'variable_names' dict
    variable_names = config["variable_names"]
    temp_var = variable_names["temperature"]           # 'thetao'
    sal_var  = variable_names["salinity"]             # 'so'
    ssh_var  = variable_names["sea_surface_height"]    # 'zos'
    u_var    = variable_names["zonal_velocity"]        # 'uo'
    v_var    = variable_names["meridional_velocity"]   # 'vo'

#    with xarray.open_dataset(ssh_file) as ds:
#        print("Original SSH data:")
#        print(ds[ssh_var])
#    ds = xarray.open_dataset(ssh_file)
#    ssh = ds["total_sea_level"]
#    print("NaN count:", ssh.isnull().sum().values)

#    ssh_roi = ssh.isel(time=0, depth=0).sel(lat=slice(0,60), lon=slice(-100,-34))
#    print("NaN count in region of interest:", ssh_roi.isnull().sum().values)

#    print("Min:", ssh_roi.min().values)
#    print("Max:", ssh_roi.max().values)




#    input("Press Enter to continue...")  # Pauses execution

     # 4) Open each NetCDF file and select the original variable name
    #    We do NOT rename to 'temp','sal','ssh','u','v' here.
    
    # Define the longitude and latitude range
    lon_min, lon_max = -101, -30
    lat_min, lat_max = 15, 52





    ds_temp = (
        xarray.open_dataset(temp_file)[temp_var]
 #       .isel(latitude=slice(None, None, 6), longitude=slice(None, None, 12))  # Downsample
        .sel(longitude=slice(lon_min, lon_max), latitude=slice(lat_min, lat_max))  # Subset region 
       .rename({"longitude": "lon", "latitude": "lat"})  # only rename coords
    )








    ds_sal = (
        xarray.open_dataset(sal_file)[sal_var]
        .sel(longitude=slice(lon_min, lon_max), latitude=slice(lat_min, lat_max))  # Subset region
#        .isel(latitude=slice(None, None, 6), longitude=slice(None, None, 12))  # Downsample
        .rename({"longitude": "lon", "latitude": "lat"})
    )

    ds_ssh = (
        xarray.open_dataset(ssh_file)[ssh_var]
        .sel(longitude=slice(lon_min, lon_max), latitude=slice(lat_min, lat_max))  # Subset region
#        .isel(latitude=slice(None, None, 6), longitude=slice(None, None, 12))  # Downsample
        .isel(time=0)
        .isel(depth=0, drop=True)
#        .sel(depth=0, method="nearest", drop=True)  # Select the surface layer and drop depth dimension
        .rename({"longitude": "lon", "latitude": "lat"})
    )

#    print("NaN count in region of interest:", ds_ssh.isnull().sum().values)
#    print("Min:", ds_ssh.min().values)
#    print("Max:", ds_ssh.max().values)

#    input("Press Enter to continue...")  # Pauses execution

    ds_u = (
        xarray.open_dataset(u_file)[u_var]
        .sel(longitude=slice(lon_min, lon_max), latitude=slice(lat_min, lat_max))  # Subset region
#        .isel(latitude=slice(None, None, 6), longitude=slice(None, None, 12))  # Downsample
        .rename({"longitude": "lon", "latitude": "lat"})
    )

    ds_v = (
        xarray.open_dataset(v_file)[v_var]
        .sel(longitude=slice(lon_min, lon_max), latitude=slice(lat_min, lat_max))  # Subset region
#        .isel(latitude=slice(None, None, 6), longitude=slice(None, None, 12))  # Downsample
        .rename({"longitude": "lon", "latitude": "lat"})
    )


    vgrid = xarray.open_dataarray(vgrid_file)
    z = vgrid_to_layers(vgrid)
    ztarget = xarray.DataArray(
        z,
        name='zl',
        dims=['zl'], 
        coords={'zl': z}, 
    )

    
    # 5) Merge into a single dataset
    glorys = xarray.merge([ds_temp, ds_sal, ds_ssh, ds_u, ds_v])
    print("Before Subsample", glorys.dims)

    # 6) Subsample: take every other point in x and y
    glorys = glorys.isel(lon=slice(0, None, 2), lat=slice(0, None, 2))
    print("After Subsample", glorys.dims)

    # Round time down to midnight
    glorys['time'] = (('time', ), ds_temp['time'].dt.floor('1d').data)


    # Interpolate GLORYS vertically onto target grid.
    # Depths below bottom of GLORYS are filled by extrapolating the deepest available value.
#    revert = glorys.interp(depth=ztarget, kwargs={'fill_value': 'extrapolate'}).ffill('zl', limit=None)
    #Interpolates only within the valid range;Leaves anything outside (e.g., deeper than GLORYS) as NaN;
    revert = glorys.interp(depth=ztarget) 
    
    # Flood temperature and salinity over land. 
    flooded = xarray.merge((
        flood.flood_kara(revert[v], zdim='zl') for v in [temp_var, sal_var, u_var, v_var]
    ))

    # flood zos separately to avoid the extra z=0 added by flood_kara.
#    flooded[ssh_var] = flood.flood_kara(revert[ssh_var]).isel(z=0).drop('z')
   
#    print("Before flooding:", revert[ssh_var])
   # Print number of NaNs before flooding
    before_nan_count = revert[ssh_var].isnull().sum().values
    print("NaN count before flooding:", before_nan_count)

    # Flood the data
    flooded_ssh = flood.flood_kara(revert[ssh_var])
#    print("After flooding:", flooded_ssh)

    # Print number of NaNs after flooding
    after_nan_count = flooded_ssh.isnull().sum().values
    print("NaN count after flooding:", after_nan_count)
    print("Min after flooding:", flooded_ssh.min().values)
    print("Max after flooding:", flooded_ssh.max().values)

    surface_ssh=flood.flood_kara(revert[ssh_var]).isel(z=0).drop_vars('z')
    surface_ssh['time'] = flooded.time
    print("surface_ssh dims:", surface_ssh.dims)
    print("surface_ssh shape:", surface_ssh.shape)
 
    print("==== flooded.time ====")
    print(flooded.time)
    print(flooded.time.values)
    print(flooded.time.dtype)

    print("==== surface_ssh.time ====")
    print(surface_ssh.time)
    print(surface_ssh.time.values)
    print(surface_ssh.time.dtype)

#    input("Press Enter to continue...")  # Pauses execution 



 
    surface_ssh_da = surface_ssh.to_dataset(name=ssh_var)

    # Merge in a way that ignores dimension mismatches:
    flooded = xarray.merge(
        [flooded, surface_ssh_da],
        compat='override'
    )


    print("NaN count (after surface selection):", surface_ssh.isnull().sum().values)
    print("Surafce ssh:", surface_ssh)
    print("Min surface_ssh:", surface_ssh.min().values)
    print("Max surface_ssh:", surface_ssh.max().values)

#    input("Press Enter to continue...")  # Pauses execution
# flooded[ssh_var] = flood.flood_kara(revert[ssh_var]).isel(z=0).drop_vars('z')    

    # Horizontally interpolate the vertically interpolated and flooded data onto the MOM grid. 
    target_grid = xarray.open_dataset(grid_file)
    
    # Adjust GLORYS longitudes to match the ocean_hgrid.nc range
    target_max_lon = target_grid['x'].max().item()
    print("Max longitude in ocean_hgrid.nc (target_max_lon):", target_max_lon)

    glorys['lon'] = xarray.where(glorys['lon'] > target_max_lon, glorys['lon'] - 360, glorys['lon'])

    
    target_t = (
        target_grid
        [['x', 'y']]
        .isel(nxp=slice(1, None, 2), nyp=slice(1, None, 2))
        .rename({'y': 'lat', 'x': 'lon', 'nxp': 'xh', 'nyp': 'yh'})
    )
    # Interpolate u and v onto supergrid to make rotation possible
    target_uv = (
        target_grid
        [['x', 'y']]
        .rename({'y': 'lat', 'x': 'lon'})
    )
    
    print("target_t", target_t)
    print("glorys", glorys)
    print("target_uv", target_uv)



    regrid_kws = dict(method='nearest_s2d', reuse_weights=reuse_weights, periodic=False)

    print("t")


    glorys_to_t = xesmf.Regridder(glorys, target_t, filename='regrid_glorys_tracers.nc', **regrid_kws)
    print("uv:")
    glorys_to_uv = xesmf.Regridder(glorys, target_uv, filename='regrid_glorys_uv.nc', **regrid_kws)

    print("GLORYS lon:", glorys["lon"].min().values, glorys["lon"].max().values)
    print("GLORYS lat:", glorys["lat"].min().values, glorys["lat"].max().values)

    print("Target_t lon:", target_t["lon"].min().values, target_t["lon"].max().values)
    print("Target_t lat:", target_t["lat"].min().values, target_t["lat"].max().values)

    print("flooded lon:", flooded["lon"].min().values, flooded["lon"].max().values)
    print("flooded lat:", flooded["lat"].min().values, flooded["lat"].max().values)

    ssh_pre = flooded[ssh_var]
    print("SSH pre-regrid min/max:", ssh_pre.min().values, ssh_pre.max().values)

#    input("Press Enter to continue...")  # Pauses execution


    interped_t = glorys_to_t(flooded[[temp_var, sal_var, ssh_var]])

    # Interpolate u and v, rotate, then extract individual u and v points
    interped_uv = glorys_to_uv(flooded[[u_var, v_var]])
    urot, vrot = rotate_uv(interped_uv[u_var], interped_uv[v_var], target_grid['angle_dx'])
    uo = urot.isel(nxp=slice(0, None, 2), nyp=slice(1, None, 2)).rename({'nxp': 'xq', 'nyp': 'yh'})
    uo.name = 'uo'
    vo = vrot.isel(nxp=slice(1, None, 2), nyp=slice(0, None, 2)).rename({'nxp': 'xh', 'nyp': 'yq'})
    vo.name = 'vo'
    
    interped = (
        xarray.merge((interped_t, uo, vo))
        .transpose('time', 'zl', 'yh', 'yq', 'xh', 'xq')
    )


    # === Merge interpolated results ===
    interped = (
        xarray.merge((interped_t, uo, vo))
        .transpose('time', 'zl', 'yh', 'yq', 'xh', 'xq')
    )

    
    # Rename to match MOM expectations.
    interped = interped.rename({
        temp_var: 'temp',
        sal_var: 'salt',
        ssh_var: 'ssh',
        u_var: 'u',
        v_var: 'v'
    })

    print("\nChunks for 'temp':")
    print(interped['temp'].chunks)

    #input("\nPress Enter to continue...")
    

    print("\n--- SHAPES BEFORE DEEP FILL ---")
    for var in ['temp']:
        if var in interped:
            da = interped[var]
            print(f"{var}: dims={da.dims}, shape={da.shape}")

    # Capture original dims for each var
    orig_dims = {
        var: interped[var].dims
        for var in ['temp','salt','u','v']
        if var in interped
    }


    # === Apply deep-ocean fill ===
    print("\nApplying deep-ocean fill...")
    for var in ['temp', 'salt', 'u', 'v']:
        if var in interped:
            print(f"Filling deep NaNs for {var}...")
            interped[var] = xarray.apply_ufunc(
                fill_from_deepest_valid,
                interped[var],
                input_core_dims=[['zl']],
                output_core_dims=[['zl']],
                vectorize=True,
                dask='parallelized',
                output_dtypes=[interped[var].dtype]
            )

             # Restore original dim‐order
            interped[var] = interped[var].transpose(*orig_dims[var])

            print("\n--- SHAPES AFTER DEEP FILL ---")
            da = interped[var]
            print(f"{var}: dims={da.dims}, shape={da.shape}")


# === Pause ===
    #input("Press Enter to continue...")










    # === Pause ===
    #input("\nDeep-ocean filling completed. Press Enter to continue...")



    # Fix output metadata, including removing all _FillValues.
    all_vars = list(interped.data_vars.keys()) + list(interped.coords.keys())
    encodings = {v: {'_FillValue': None} for v in all_vars}
    encodings['time'].update({'dtype':'float64', 'calendar': 'gregorian'})
    interped['zl'].attrs = {
        'long_name': 'Layer pseudo-depth, -z*',
         'units': 'meter',
         'cartesian_axis': 'Z',
         'positive': 'down'
    }

    # Extract the directory from the output_file pat
    output_folder = os.path.dirname(output_file)
    
    # Check if the output folder exists, and if not, create it
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)   

    print("Variables in final dataset:", list(interped.data_vars))
    print("SSH stats:", interped["ssh"].min().values, interped["ssh"].max().values)
    #input("Press Enter to continue...")  # Pauses execution


    # output results
    interped.to_netcdf(
        output_file,
        format='NETCDF3_64BIT',
        engine='netcdf4',
        encoding=encodings,
        unlimited_dims='time'
    )


def main():

    parser = argparse.ArgumentParser(description='Generate ICs from Glorys.')
    parser.add_argument('--config_file', type=str, default='glorys_ic.yaml' , help='Path to the YAML config file')
    args = parser.parse_args()

    if not args.config_file:
        parser.error('Please provide the path to the YAML config file.')

    with open(args.config_file, 'r') as yaml_file:
        config = yaml.safe_load(yaml_file)

    if not all(key in config for key in [
        'glorys_temperature',
        'glorys_salinity',
        'glorys_sea_surface_height',
        'glorys_zonal_velocity',
        'glorys_meridional_velocity',
        'vgrid_file',
        'grid_file',
        'output_file'
    ]):
        parser.error('Please provide all required parameters in the YAML config file.')


 #   if not all(key in config for key in ['glorys_file', 'vgrid_file', 'grid_file', 'output_file']):
 #       parser.error('Please provide all required parameters in the YAML config file.')

    write_initial(config)

if __name__ == '__main__':
    main()
