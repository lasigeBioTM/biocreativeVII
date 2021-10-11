import ontology_preprocessing

from fuzzywuzzy import process
from fuzzywuzzy import fuzz

is_a_graph, name_to_id, synonym_to_id, id_to_name, id_to_index = ontology_preprocessing.load_go()

all_data = ontology_preprocessing.load_data()

entities_info = {}
with open('corpora/training/drugprot_training_entities.tsv', encoding='utf-8') as ef:
    for line in ef:
        line = line.strip()
        line = tuple(line.split('\t'))
        if line[5] in entities_info:
            pass
        else:
            entities_info[line[5]] = line[2]

for k in entities_info:
    if 'GENE' in entities_info[k]:
        if k in name_to_id:
            all_data[k] = name_to_id[k]
        elif k in synonym_to_id:
            all_data[k] = synonym_to_id[k]
        else: #add syn
            print('NO ID')
            e = process.extractOne(k, name_to_id.keys(), scorer = fuzz.token_sort_ratio)
            if e[0] in name_to_id:
                entity_id = name_to_id[e[0]]
            elif e[0] in synonym_to_id:
                entity_id = synonym_to_id[e[0]][0]
            all_data[k] = entity_id
        if type(all_data[k]) is list:
            all_data[k] = all_data[k][0]

with open('corpora/training/out2.txt', mode='w', encoding='utf-8') as out_file:
    for key, value in all_data.items(): 
        out_file.write('%s %s\n' % (key, value))
    out_file.close()

