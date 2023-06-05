import glob
import shutil
base_path = '/home/brendan/research/FINAL_LIGANDS'
# for each file in ligand subdirectory
for file in glob.glob(f'{base_path}/**/**/*.mol'):
    # print(file)
    # get file name, left strip U, right strip *.xyz, readd .xyz 
    _,_,_,_,_,denticity,connection_type,core_file_name = file.split(sep='/')
    # print(core_file_name)
    new_file_name = core_file_name[1:]
    # print(new_file_name)
    try:
        csd_id,_,_,_,_ = new_file_name.split(sep='_')
        if len(csd_id) != 6:
            print('shortened id!')
        new_file_name = new_file_name[:-5] + '.mol'
        new_file_path = f'{base_path}/{denticity}/{connection_type}/{new_file_name}'
        # print(new_file_path)
        shutil.move(file, new_file_path)
    except ValueError:
        # print('Value error: ', file)
        pass

    
    # rename file