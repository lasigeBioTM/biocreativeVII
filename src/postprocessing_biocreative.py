from itertools import product, combinations
from parse_biocreative_vii import get_new_offsets_sentences


def lookup_dict(entities_file: str):
    """

    :param entities_file:
    :return:
    """

    file_entities = open(entities_file, 'r', encoding='utf-8')
    entities = file_entities.readlines()
    file_entities.close()

    output_pmids = open('pmids.txt', 'w', encoding='utf-8')

    abstracts = []

    dict_entities = {}

    for entity_elements in entities:
        abstract = entity_elements.split('\t')[0]
        entity_id = entity_elements.split('\t')[1]
        entity_type = entity_elements.split('\t')[2]
        entity_name = entity_elements.split('\t')[-1][:-1]

        abstracts.append(abstract)
        unique_id = abstract + entity_id

        dict_entities[unique_id] = [entity_type, entity_name]

    for abstract in set(abstracts):
        output_pmids.write(abstract + '\n')

    output_pmids.close()

    return dict_entities, set(abstracts)


def rules_fn(entities_file: str):

    dict_entities, abstracts = lookup_dict(entities_file)

    save_entities = {}

    for entity, entity_elements in dict_entities.items():
        abstract = entity.split('T')[0]

        if abstract not in save_entities:
            save_entities[abstract] = [[], []]

            if entity_elements[0] == 'CHEMICAL':
                save_entities[abstract][0].append(entity_elements[1])
            elif entity_elements[0] == 'GENE':
                save_entities[abstract][1].append(entity_elements[1])

        else:
            if entity_elements[0] == 'CHEMICAL':
                save_entities[abstract][0].append(entity_elements[1])
            elif entity_elements[0] == 'GENE':
                save_entities[abstract][1].append(entity_elements[1])

    save_pairs = {}  # 17512723: [(retinol, retinol dehydrogenase), ...]

    for abstract, entities in save_entities.items():

        save_pairs[abstract] = []
        chemical_list = entities[0]
        gene_list = entities[1]

        combinations = product(chemical_list, gene_list)

        for combination in combinations:

            for chemical_in_list in chemical_list:
                if combination[0] in chemical_in_list: #and combination[0] != chemical_in_list:
                    if (combination[0], combination[1]) not in save_pairs[abstract]:
                        save_pairs[abstract].append((combination[0], combination[1]))

            for gene_in_list in gene_list:
                if combination[0] in gene_in_list: #and combination[0] != gene_in_list:
                    if (combination[0], combination[1]) not in save_pairs[abstract]:
                        save_pairs[abstract].append((combination[0], combination[1]))

    return save_pairs


#print(rules_fn(entities_file='../drugprot-training-development-test-background/drugprot-gs-training-development/development/drugprot_development_entities.tsv'))


def txt_to_tsv(abstracts_file: str, entities_file: str, input_file: str, output_file: str):
    """

    :param abstracts_file:
    :param entities_file:
    :param input_file:
    :param output_file:
    :return:
    """

    dict_entities, abstracts = lookup_dict(entities_file)
    rules = rules_fn(entities_file)

    rules_relations = rules

    save_written = []

    input_txt = open(input_file, 'r', encoding='utf-8')
    input_txt.readline()
    input_txt_lines = input_txt.readlines()
    input_txt.close()

    output_tsv = open(output_file, 'w', encoding='utf-8')

    for line in input_txt_lines:
        e1 = line.split('\t')[0]
        e2 = line.split('\t')[1]
        relation = line.split('\t')[2][:-1]

        abstract = e1.split('.')[0][1:]
        entity_1 = e1.split('u')[-1].split('.')[0]
        entity_2 = e2.split('u')[-1].split('.')[0]

        unique_id_1 = abstract + entity_1
        unique_id_2 = abstract + entity_2

        if dict_entities[unique_id_1][0] == dict_entities[unique_id_2][0]:
            continue

        else:
            if relation != 'NO_RELATION':
                if dict_entities[unique_id_1][0] == 'GENE':
                    output_tsv.write(abstract + '\t' + relation + '\t' + 'Arg1:' + entity_2 + '\t' + 'Arg2:' + entity_1 + '\n')
                    save_written.append([abstract, relation, entity_2, entity_1])

                    if (dict_entities[unique_id_2][1], dict_entities[unique_id_1][1]) in rules[abstract]:
                        rules_relations[abstract].remove((dict_entities[unique_id_2][1], dict_entities[unique_id_1][1]))
                        rules_relations[abstract].append((dict_entities[unique_id_2][1], dict_entities[unique_id_1][1], relation))

                else:
                    output_tsv.write(abstract + '\t' + relation + '\t' + 'Arg1:' + entity_1 + '\t' + 'Arg2:' + entity_2 + '\n')
                    save_written.append([abstract, relation, entity_1, entity_2])

                    if (dict_entities[unique_id_1][1], dict_entities[unique_id_2][1]) in rules[abstract]:
                        rules_relations[abstract].remove((dict_entities[unique_id_1][1], dict_entities[unique_id_2][1]))
                        rules_relations[abstract].append((dict_entities[unique_id_1][1], dict_entities[unique_id_2][1], relation))

            # if abstract == '23123662':
            #     print(rules_relations['23123662'])
            #     print(len(rules_relations['23123662']))

    # WITH FN RULES FOR NO RELATION

    for line in input_txt_lines:
        e1 = line.split('\t')[0]
        e2 = line.split('\t')[1]
        relation = line.split('\t')[2][:-1]

        if relation != 'NO_RELATION':
            continue

        abstract = e1.split('.')[0][1:]
        entity_1 = e1.split('u')[-1].split('.')[0]
        entity_2 = e2.split('u')[-1].split('.')[0]

        unique_id_1 = abstract + entity_1
        unique_id_2 = abstract + entity_2

        if dict_entities[unique_id_1][0] == dict_entities[unique_id_2][0]:
            continue

        else:
            if dict_entities[unique_id_1][0] == 'GENE':
                for pair in rules_relations[abstract]:
                    if dict_entities[unique_id_2][1] == pair[0] and dict_entities[unique_id_1][1] == pair[1]:
                        if len(pair) == 3:
                            relation = pair[2]
                            output_tsv.write(abstract + '\t' + relation + '\t' + 'Arg1:' + entity_2 + '\t' + 'Arg2:' + entity_1 + '\n')
                            save_written.append([abstract, relation, entity_2, entity_1])

                        continue

            else:
                for pair in rules_relations[abstract]:
                    if dict_entities[unique_id_1][1] == pair[0] and dict_entities[unique_id_2][1] == pair[1]:
                        if len(pair) == 3:
                            relation = pair[2]
                            output_tsv.write(abstract + '\t' + relation + '\t' + 'Arg1:' + entity_1 + '\t' + 'Arg2:' + entity_2 + '\n')
                            save_written.append([abstract, relation, entity_1, entity_2])

                        continue

    # WITH SENTENCE_ENTITIES FOR NOT IN OUTPUT FILE

    abstracts_file = open(abstracts_file, 'r', encoding='utf-8')
    abstracts_lines = abstracts_file.readlines()
    abstracts_file.close()

    for abstract in abstracts_lines:
        abstract_id = abstract.split('\t')[0]
        entities_per_sentence = get_new_offsets_sentences(abstract_id)

        for sentence, entities in entities_per_sentence.items():
            abstract = sentence[0].split('.')[0][1:]

            entities_combinations = combinations([(e[3], e[0]) for e in entities], 2)

            pairs = rules_relations[abstract]
            for entities_combination in entities_combinations:
                for pair in pairs:

                    if entities_combination[0][0] == pair[0] and entities_combination[1][0] == pair[1]:
                        if len(pair) == 3:
                            if [abstract, pair[2], entities_combination[0][1], entities_combination[1][1]] not in save_written:
                                output_tsv.write(abstract + '\t' + pair[2] + '\t' + 'Arg1:' + entities_combination[0][1] + '\t' + 'Arg2:' + entities_combination[1][1] + '\n')

                    elif entities_combination[0][0] == pair[1] and entities_combination[1][0] == pair[0]:
                        if len(pair) == 3:
                            if [abstract, pair[2], entities_combination[1][1], entities_combination[0][1]] not in save_written:
                                output_tsv.write(abstract + '\t' + pair[2] + '\t' + 'Arg1:' + entities_combination[1][1] + '\t' + 'Arg2:' + entities_combination[0][1] + '\n')

    output_tsv.close()

    return

# DEV

# txt_to_tsv('../drugprot-training-development-test-background/drugprot-gs-training-development/test-background/test_background_abstracts.tsv',
#             '../drugprot-training-development-test-background/drugprot-gs-training-development/test-background/test_background_entities.tsv',
#            'run_23_09_test/model_class_weights_low_drug_gene_results.txt',
#            'run_23_09_test/1-lasigebiotm_deep_learning_ontologies_model_class_weights_lower.tsv')

# txt_to_tsv('../drugprot-training-development-test-background/drugprot-gs-training-development/test-background/test_background_abstracts.tsv',
#             '../drugprot-training-development-test-background/drugprot-gs-training-development/test-background/test_background_entities.tsv',
#            'run_23_09_test/model_no_ancestors_3_drug_gene_results.txt',
#            'run_23_09_test/2-lasigebiotm_deep_learning_ontologies_model_no_ancestors_precision.tsv')

# txt_to_tsv('../drugprot-training-development-test-background/drugprot-gs-training-development/test-background/test_background_abstracts.tsv',
#             '../drugprot-training-development-test-background/drugprot-gs-training-development/test-background/test_background_entities.tsv',
#            'run_23_09_test/model_no_ancestors_class_weights_low_drug_gene_results.txt',
#            'run_23_09_test/3-lasigebiotm_deep_learning_ontologies_model_no_ancestors_class_weights_lower.tsv')

# txt_to_tsv('../drugprot-training-development-test-background/drugprot-gs-training-development/test-background/test_background_abstracts.tsv',
#             '../drugprot-training-development-test-background/drugprot-gs-training-development/test-background/test_background_entities.tsv',
#            'run_23_09_test/model_no_ancestors_1_drug_gene_results.txt',
#            'run_23_09_test/4-lasigebiotm_deep_learning_ontologies_model_no_ancestors_score.tsv')

txt_to_tsv('../drugprot-training-development-test-background/drugprot-gs-training-development/test-background/test_background_abstracts.tsv',
            '../drugprot-training-development-test-background/drugprot-gs-training-development/test-background/test_background_entities.tsv',
           'run_23_09_test/model_no_ancestors_1_class_weights_drug_gene_results.txt',
           'run_23_09_test/5-lasigebiotm_deep_learning_ontologies_model_no_ancestors_class_weights.tsv')


# txt_to_tsv('../drugprot-training-development-test-background/drugprot-gs-training-development/development/drugprot_development_abstracs.tsv',
#            '../drugprot-training-development-test-background/drugprot-gs-training-development/development/drugprot_development_entities.tsv',
#            'run_22_09_dev_class_weights/model_no_ancestors_1_class_weights_drug_gene_results.txt',
#            'lasigebiotm_deep_learning_ontologies_model_no_ancestors_1_class_weights_dev_rules.tsv')

# txt_to_tsv('../drugprot-training-development-test-background/drugprot-gs-training-development/development/drugprot_development_abstracs.tsv',
#            '../drugprot-training-development-test-background/drugprot-gs-training-development/development/drugprot_development_entities.tsv',
#            'run_23_09_dev_low/model_no_ancestors_class_weights_low_drug_gene_results.txt',
#            'lasigebiotm_deep_learning_ontologies_model_no_ancestors_class_weights_low_dev_rules.tsv')