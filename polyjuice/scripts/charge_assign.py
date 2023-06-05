import pandas as pd 
import glob
import os 

def update_header(file_name, charge, confidence):
    with open(file_name, 'r') as f:
        charge_data = 'q = ' + str(charge) + ', confidence = ' + str(confidence)
        file_data = f.readlines()
        file_data[1] = charge_data + '\n'

    with open(file_name, 'w') as file:
        file.writelines(file_data)

df = pd.read_csv('./analysis/predictions.csv')

# make column for new names
df['new_names'] = df['Name'].str[1:-5] + '.mol'

ligands_not_found = []
base_path = '/home/brendan/research/FINAL_LIGANDS'
count = 0
for file in glob.glob('/home/brendan/research/FINAL_LIGANDS/**/**/*.mol'):
    try: 
        # print(file)
        _,_,_,_,_,denticity,connection_type,core_file_name = file.split(sep='/')
        # print(core_file_name)
        info = df.loc[df['new_names'] == core_file_name]
        # print(info)
        charge = float(info['charge'].iloc[0])
        # print(charge)
        confidence = float(info['confidence'].iloc[0])
        # print(confidence)
        update_header(file, charge, confidence)

        confidence_percentage = str(int(round(confidence, ndigits=2) * 100))
        # print(confidence_percentage)

        # rename file with attached confidence level
        file_base = core_file_name[:-4] 
        # print(file_base)
        new_file_name = f'{base_path}/{denticity}/{connection_type}/{file_base}_{confidence_percentage}.mol'
        # print(new_file_name)

        os.rename(file, new_file_name)
    except IndexError:
        count += 1
        print('Ligand not found in charge sheet')
        ligands_not_found.append(file)
    except ValueError:
        print("value error")
    
print(f'{count} ligands from directories not found in charge predicitions sheet')