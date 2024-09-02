import yaml
import pandas as pd
import json
import shutil
import subprocess

output_path = 'resources/SRA/'


# Load the YAML file
with open('config/ecoli.yaml', 'r') as file:
    yaml_data = yaml.safe_load(file)

def read_json_file(file_path):
    with open(file_path, 'r') as json_file:
        data = json.load(json_file)
    return data

# Use the function
used_samples = read_json_file('used_samples.json')

# Convert the YAML data to a pandas DataFrame
# Assuming the YAML structure is a list of dictionaries
# for sample in yaml_data["sample"]:


# Convert the YAML data (nested dictionaries) to a list of dictionaries
data = []
for key, value in yaml_data['meta'].items():
    if key in yaml_data['sample'] and value['title'] in used_samples:
        entry = {'id': key}  # Add the key as an 'id' column
        entry.update(value)
        data.append(entry)

# Convert the list of dictionaries to a pandas DataFrame
df = pd.DataFrame(data)

df = df.drop(['medium', 'title', 'damageSite'], axis=1)
# df['read1'] = df['id'].apply(lambda x: f'{x}.fastq')
df['read1'] = df['id'].apply(lambda x: '{}.fastq'.format(x))
df['read2'] = 'NA'
df['technology'] = 'XR-seq'
df['species'] = 'E. coli'

df_repair = df

# for index, row in df_repair.iterrows():
#     subprocess.run(f'gzip -c resources/samples/{row["read1"]} > {output_path}{row["read1"].split("/")[-1]}.gz', shell=True)

# Load the YAML file
with open('config/ecoli_rna_custom.yaml', 'r') as file:
    yaml_data = yaml.safe_load(file)

# Use the function
used_samples = read_json_file('used_samples_rna.json')

# Convert the YAML data (nested dictionaries) to a list of dictionaries
data = []
for key, value in yaml_data['meta'].items():
    if key in yaml_data['sample'] and value['title'] in used_samples:
        entry = {'id': key}  # Add the key as an 'id' column
        entry.update(value)
        data.append(entry)

# Convert the list of dictionaries to a pandas DataFrame
df = pd.DataFrame(data)

df = df.drop(['treatment', 'title'], axis=1)
df['read1'] = df['id'].apply(lambda x: '{}_1.fq.gz'.format(x))
df['read2'] = df['id'].apply(lambda x: '{}_2.fq.gz'.format(x))
# df['read1'] = df['id'].apply(lambda x: f'{x}_1.fq.gz')
# df['read2'] = df['id'].apply(lambda x: f'{x}_2.fq.gz')
df['technology'] = 'RNA-seq'
df['species'] = 'E. coli'

df_transcriptome = df

# for index, row in df_transcriptome.iterrows():
#     shutil.copy('resources/RNA-seq/ecoli/' + row['id'] + '/' + row['read1'], output_path)
#     shutil.copy('resources/RNA-seq/ecoli/' + row['id'] + '/' + row['read2'], output_path)


df = pd.concat([df_repair, df_transcriptome])

print(df)

# Save the DataFrame to a CSV file
df.to_csv('samples.csv', index=False)
