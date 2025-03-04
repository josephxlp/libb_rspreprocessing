import os 
import yaml
from .rgrid import format_tile_fpath, gdal_regrid, get_raster_info
from .rtransforms import scale_raster,raster_calc
from .rfilter import (classify_lwm_CopWBM, classify_lwm_CopWBM, classify_lwm_TanDEMX_LCM,
                     filter_tandemx_noise,filter_water,combine_water_masks,classify_lwm_ESAWC)


# vars = (tdem_dem_fpath,tdem_hem_fpath,tdem_wam_fpath,cdem_wbm_fpath,cdem_wbm_fpath,
#             dtm_fpath,pdem_fpath,edem_egm_fpath,edem_wgs_fpath,edem_lcm_fpath,
#             esawc_fpath,etchm_fpath,etchm_fpath,egm08_fpath,s1_fpath,s2_fpath)

def retile_datasets(
    ds_tiles_dpath, tilename, xmin, ymin, xmax, ymax, xres, yres,
    tdem_dem_fpath, tdem_hem_fpath, tdem_wam_fpath, cdem_wbm_fpath,
    dtm_fpath, pdem_fpath, edem_egm_fpath, edem_wgs_fpath,
    edem_lcm_fpath, esawc_fpath, etchm_fpath, egm08_fpath,
    s1_fpath,s2_fpath,tdem_com_fpath,fbchm_fpath
):
    tilename_dpath = os.path.join(ds_tiles_dpath, tilename)
    os.makedirs(tilename_dpath,exist_ok=True)
    nmode = "num" 
    cmode = "cat"
    ds = {}
    s1_tile = format_tile_fpath(tilename_dpath, tilename, s1_fpath)
    gdal_regrid(s1_fpath, s1_tile, xmin, ymin, xmax, ymax, xres, yres, mode=nmode)
    ds["s1"] = s1_tile
    
    print(s2_fpath)
    s2_tile = format_tile_fpath(tilename_dpath, tilename, s2_fpath)
    gdal_regrid(s2_fpath, s2_tile, xmin, ymin, xmax, ymax, xres, yres, mode=nmode)
    ds["s2"] = s2_tile

    tdem_dem_tile = format_tile_fpath(tilename_dpath, tilename, tdem_dem_fpath)
    gdal_regrid(tdem_dem_fpath, tdem_dem_tile, xmin, ymin, xmax, ymax, xres, yres, mode='num')
    ds['tdem_dem'] = tdem_dem_tile

    tdem_hem_tile = format_tile_fpath(tilename_dpath, tilename, tdem_hem_fpath)
    gdal_regrid(tdem_hem_fpath, tdem_hem_tile, xmin, ymin, xmax, ymax, xres, yres, mode='num')
    ds['tdem_hem'] = tdem_hem_tile

    tdem_wam_tile = format_tile_fpath(tilename_dpath, tilename, tdem_wam_fpath)
    gdal_regrid(tdem_wam_fpath, tdem_wam_tile, xmin, ymin, xmax, ymax, xres, yres, mode='cat')
    ds['tdem_wam'] = tdem_wam_tile

    tdem_com_tile = format_tile_fpath(tilename_dpath, tilename, tdem_com_fpath)
    gdal_regrid(tdem_com_fpath, tdem_com_tile, xmin, ymin, xmax, ymax, xres, yres, mode='cat')
    ds['tdem_com'] = tdem_com_tile

    cdem_wbm_tile = format_tile_fpath(tilename_dpath, tilename, cdem_wbm_fpath)
    gdal_regrid(cdem_wbm_fpath, cdem_wbm_tile, xmin, ymin, xmax, ymax, xres, yres, mode='cat')
    ds['cdem_wbm'] = cdem_wbm_tile

    edem_lcm_tile = format_tile_fpath(tilename_dpath, tilename, edem_lcm_fpath)
    gdal_regrid(edem_lcm_fpath, edem_lcm_tile, xmin, ymin, xmax, ymax, xres, yres, mode='cat')
    ds['edem_lcm'] = edem_lcm_tile

    esawc_tile = format_tile_fpath(tilename_dpath, tilename, esawc_fpath)
    gdal_regrid(esawc_fpath, esawc_tile, xmin, ymin, xmax, ymax, xres, yres, mode='cat')
    ds['esawc'] = esawc_tile

    etchm_tile = format_tile_fpath(tilename_dpath, tilename, etchm_fpath)
    gdal_regrid(etchm_fpath, etchm_tile, xmin, ymin, xmax, ymax, xres, yres, mode='cat')
    ds['etchm'] = etchm_tile

    fbchm_tile = format_tile_fpath(tilename_dpath, tilename, fbchm_fpath)
    gdal_regrid(fbchm_fpath, fbchm_tile, xmin, ymin, xmax, ymax, xres, yres, mode='cat')
    ds['fbchm'] = fbchm_tile
    



    ####################################################################
    ####################################################################
    # add if already exist 
    wbm_lwm_fn = cdem_wbm_tile.replace('cdem_WBM.tif', 'cdem_WBM_LWM.tif')
    if not os.path.isfile(wbm_lwm_fn):
        classify_lwm_CopWBM(cdem_wbm_tile, wbm_lwm_fn)
    ds['cdem_wbm'] = wbm_lwm_fn

    esa_lwm_fn = esawc_tile.replace('_multi_ESAWC.tif', '_ESAWC_LWM.tif')
    if not os.path.isfile(esa_lwm_fn):
        classify_lwm_ESAWC(esawc_tile, esa_lwm_fn,water_code=80)
    ds['esawc_lwm'] = esa_lwm_fn
 
    lwm_a_fn = edem_lcm_tile.replace('.tif', '_LWM_A.tif')
    lwm_b_fn = edem_lcm_tile.replace('.tif', '_LWM_B.tif')
    if not os.path.isfile(lwm_a_fn):
        classify_lwm_TanDEMX_LCM(edem_lcm_tile,lwm_a_fn,lwm_b_fn)
        #elcm_fna, elcm_fnb = classify_lwm_TanDEMX_LCM(edem_lcm_tile,lwm_a_fn,lwm_b_fn)

    ds['edem_lcm_lwma'] = lwm_a_fn
    ds['edem_lcm_lwmb'] = lwm_b_fn
    lcm_lwm_fn = lwm_a_fn

    fdem_file = tdem_dem_tile.replace('.tif', '_F.tif')
    mask_file = tdem_dem_tile.replace('.tif', '_M.tif')

    dem_fw = fdem_file.replace('F.tif', '_Fw.tif')
    dem_mw = mask_file.replace('M.tif', '_Mw.tif')

    if not os.path.isfile(dem_fw):
        fdem_file, mask_file = filter_tandemx_noise(tdem_dem_tile, 
                                                    tdem_hem_tile, 
                                                    tdem_com_tile, 
                                                    n_iter=1)  
    if not os.path.isfile(dem_fw):
        filter_water(fdem_file, mask_file, lcm_lwm_fn, esa_lwm_fn, wbm_lwm_fn)

    ds['tdem_dem_fw'] = dem_fw
    ds['tdem_dem_mw'] = dem_mw

    tname = os.path.basename(lcm_lwm_fn).split('_')[0]
    combined_mask_file = os.path.join(os.path.dirname(lcm_lwm_fn), f'{tname}_LWM.tif')
    if not os.path.isfile(combined_mask_file):
        combined_mask_file = combine_water_masks(lcm_lwm_fn, esa_lwm_fn, wbm_lwm_fn)
    ds['lwm'] = combined_mask_file

    ####################################################################
    ####################################################################

    ldar_tile = format_tile_fpath(tilename_dpath, tilename, dtm_fpath)
    gdal_regrid(dtm_fpath, ldar_tile, xmin, ymin, xmax, ymax, xres, yres, mode='num')
    ds['ldtm'] = ldar_tile

    pdem_tile = format_tile_fpath(tilename_dpath, tilename, pdem_fpath)
    gdal_regrid(pdem_fpath, pdem_tile, xmin, ymin, xmax, ymax, xres, yres, mode='num')
    ds['pdem'] = pdem_tile

    # cdem_dem_tile = format_tile_fpath(tilename_dpath, tilename, cdem_dem_fpath)
    # gdal_regrid(cdem_dem_fpath, cdem_dem_tile, xmin, ymin, xmax, ymax, xres, yres, mode='num')
    # ds['cdem_dem'] = cdem_dem_tile

    edem_egm_tile = format_tile_fpath(tilename_dpath, tilename, edem_egm_fpath)
    gdal_regrid(edem_egm_fpath, edem_egm_tile, xmin, ymin, xmax, ymax, xres, yres, mode='num')
    ds['edem_egm'] = edem_egm_tile

    edem_wgs_tile = format_tile_fpath(tilename_dpath, tilename, edem_wgs_fpath)
    gdal_regrid(edem_wgs_fpath, edem_wgs_tile, xmin, ymin, xmax, ymax, xres, yres, mode='num')
    ds['edem_wgs'] = edem_wgs_tile

    egm08_tile = format_tile_fpath(tilename_dpath, tilename, egm08_fpath)
    gdal_regrid(egm08_fpath, egm08_tile, xmin, ymin, xmax, ymax, xres, yres, mode='num')
    ds['egm08'] = egm08_tile

    ####################################################################
    ######### EXTRACT AND SHIFT GEOID EGM08
    ####################################################################
    edem_geoid_tile = f'{tilename_dpath}/{tilename}_EDEM_GEOID.tif'
    raster_calc(edem_egm_tile, edem_wgs_tile, "subtract", edem_geoid_tile)
    ds['edem_goid'] = edem_geoid_tile

    egm08p_tile = f'{tilename_dpath}/{tilename}_EGM08p.tif'
    ds['egm08p'] = egm08p_tile

    yaml_tile = os.path.join(tilename_dpath, f'{tilename}_ds.yaml')
    write_yaml(ds, yaml_tile)
    print('yaml_tile:', yaml_tile)

    
def write_yaml(yaml_data, yaml_file_path):
    with open(yaml_file_path, 'w') as yaml_file:
        yaml.dump(yaml_data, yaml_file)
    print('write_yaml')
    
def read_yaml(filename):
    with open(filename, 'r') as yaml_file:
        yaml_data = yaml.safe_load(yaml_file)
        print('read_yaml')
        return yaml_data
    


def get_tilename_from_tdem_basename(fpath):
    return os.path.basename(fpath).split('_')[4]

def process_tile(
    basefile, ds_tiles_dpath,
    tdem_dem_fpath, tdem_hem_fpath, tdem_wam_fpath, cdem_wbm_fpath,
    dtm_fpath, pdem_fpath, edem_egm_fpath, edem_wgs_fpath,
    edem_lcm_fpath, esawc_fpath, etchm_fpath, egm08_fpath,
    s1_fpath,s2_fpath,tdem_com_fpath,fbchm_fpath
):
    tilename = get_tilename_from_tdem_basename(basefile)
    print(os.path.basename(basefile))
    print(tilename)
    tilename_dpath = os.path.join(ds_tiles_dpath, tilename)
    os.makedirs(tilename_dpath,exist_ok=True)
    tile_fpath = format_tile_fpath(tilename_dpath, tilename, tdem_dem_fpath) 
    proj, xres, yres, xmin, xmax, ymin, ymax, w, h = get_raster_info(basefile)
    #print('dst size:', w, h)
    retile_datasets(
        ds_tiles_dpath, tilename, xmin, ymin, xmax, ymax, xres, yres,
        tdem_dem_fpath, tdem_hem_fpath, tdem_wam_fpath, cdem_wbm_fpath,
        dtm_fpath, pdem_fpath, edem_egm_fpath, edem_wgs_fpath,
        edem_lcm_fpath, esawc_fpath, etchm_fpath, egm08_fpath,
        s1_fpath,s2_fpath,tdem_com_fpath,fbchm_fpath
    )

    
    