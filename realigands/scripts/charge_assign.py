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
base_path = ''
count = 0
for file in glob.glob('/**/**/*.mol'):
    try: 
        _,_,_,_,_,denticity,connection_type,core_file_name = file.split(sep='/')
        info = df.loc[df['new_names'] == core_file_name]
        charge = float(info['charge'].iloc[0])
        confidence = float(info['confidence'].iloc[0])
        update_header(file, charge, confidence)

        confidence_percentage = str(int(round(confidence, ndigits=2) * 100))

        # rename file with attached confidence level
        file_base = core_file_name[:-4] 
        new_file_name = f'{base_path}/{denticity}/{connection_type}/{file_base}_{confidence_percentage}.mol'

        os.rename(file, new_file_name)
    except IndexError:
        count += 1
        print('Ligand not found in charge sheet')
        ligands_not_found.append(file)
    except ValueError:
        print("value error")
    
print(f'{count} ligands from directories not found in charge predicitions sheet')