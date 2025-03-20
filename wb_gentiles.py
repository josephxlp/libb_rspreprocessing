import os 
import time 
import multiprocessing
from pprint import pprint
from config.uvars import brchieve_dpath12
from src.ub_tiles import process_tile

# f1path = "/media/ljp238/12TBWolf/ARCHIEVE/ARCHIVE_VRT/format_variable_and_files.yaml"
# f2path = "/media/ljp238/12TBWolf/ARCHIEVE/ARCHIVE_VRT/vars_vrts.yaml"

f1path = "/media/ljp238/12TBWolf/ARCHIEVE/ARCHIVE_VRTV2/loadfiles_byvariable.yaml"
f2path = "/media/ljp238/12TBWolf/ARCHIEVE/ARCHIVE_VRTV2/vars_txts.yaml"
#=============================utils
import yaml

def write_yaml(yaml_data, yaml_file_path):
    with open(yaml_file_path, 'w') as yaml_file:
        yaml.dump(yaml_data, yaml_file)
    print('write_yaml')
    
def read_yaml(filename):
    with open(filename, 'r') as yaml_file:
        yaml_data = yaml.safe_load(yaml_file)
        print('read_yaml')
        return yaml_data

#=============================utils

if __name__ == "__main__":
    ti = time.perf_counter()
    bpaths = read_yaml(f1path) 
    gpaths = read_yaml(f2path)
    os.makedirs(brchieve_dpath12,exist_ok=True)
    os.chdir(brchieve_dpath12)

    basefiles = bpaths["tdem_dem"]
    print(f"basefiles: {len(basefiles)}")

    tdem_dem_fpath = gpaths["tdem_dem"]
    tdem_hem_fpath = gpaths["tdem_hem"]
    tdem_wam_fpath = gpaths["tdem_wam"]  
    tdem_com_fpath = gpaths["tdem_com"]
    

    dtm_fpath = gpaths["ldem"] # fix lidar in VRT
    #dsm_fpath = gpaths['multi_DSM_LiDAR']
    pdem_fpath = gpaths["pdem"]
    #cdem_dem_fpath = gpaths['cdem_DEM'] to be downloaded by CDEMx
    edem_egm_fpath = gpaths["edem_egm"]
    edem_wgs_fpath = gpaths["edem_wgs"]
    edem_lcm_fpath = gpaths["edem_lcm"]
    
    esawc_fpath = gpaths["esawc"]
    etchm_fpath = gpaths["etchm"]
    
    egm08_fpath = gpaths["egm08"]
    

    cdem_wbm_fpath = gpaths["cdem_wbm"] 
    s1_fpath = gpaths["s1"]
    s2_fpath = gpaths["s2"]
    # why is gdtm and gdsm here 
    fbchm_fpath = gpaths["fbchm"]
    fbcha_fpath = gpaths["fbcha"]
    lgeoid_fpath = gpaths["lgeoid"]
    glow_fpath = gpaths["gedi_dtm"]
    gtop_fpath = gpaths["gedi_dsm"]
    wsfbh_fpath = gpaths["wsfbh"]

    num_processes = int(multiprocessing.cpu_count() * 0.75)
    pool = multiprocessing.Pool(processes=num_processes)

    vars = (tdem_dem_fpath,tdem_hem_fpath,tdem_wam_fpath,cdem_wbm_fpath,cdem_wbm_fpath,
            dtm_fpath,pdem_fpath,edem_egm_fpath,edem_wgs_fpath,edem_lcm_fpath,
            esawc_fpath,etchm_fpath,etchm_fpath,egm08_fpath,s1_fpath,s2_fpath,
            tdem_com_fpath,fbchm_fpath,fbcha_fpath,glow_fpath,gtop_fpath,lgeoid_fpath,
            wsfbh_fpath)

    for i, basefile in enumerate(basefiles):
        print(f"{i}/{len(basefiles)} @{basefile}")
        pool.apply_async(
            process_tile, (basefile, brchieve_dpath12,
                           tdem_dem_fpath, tdem_hem_fpath, tdem_wam_fpath, cdem_wbm_fpath,
                           dtm_fpath, pdem_fpath, edem_egm_fpath, edem_wgs_fpath,
                           edem_lcm_fpath, esawc_fpath, etchm_fpath, egm08_fpath,
                           s1_fpath,s2_fpath,tdem_com_fpath,fbchm_fpath,
                           fbcha_fpath,glow_fpath,gtop_fpath,lgeoid_fpath,wsfbh_fpath)
        )
    pool.close()
    pool.join()

    # for i, basefile in enumerate(basefiles):
    #     print(f"{i}/{len(basefiles)} @{basefile}")
    #     if i > 0 : break
    
    #     process_tile(basefile, brchieve_dpath12,
    #                        tdem_dem_fpath, tdem_hem_fpath, tdem_wam_fpath, cdem_wbm_fpath,
    #                        dtm_fpath, pdem_fpath, edem_egm_fpath, edem_wgs_fpath,
    #                        edem_lcm_fpath, esawc_fpath, etchm_fpath, egm08_fpath,
    #                        s1_fpath,s2_fpath,tdem_com_fpath)

    print("All tasks completed")
    tf = time.perf_counter() - ti
    print(f'run.time: {tf/60} min(s)')
    

 
    
    



    