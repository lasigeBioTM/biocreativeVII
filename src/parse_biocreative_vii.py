from __future__ import unicode_literals, print_function

import multiprocessing
import os
import sys
import logging
import random

from spacy.lang.en import English
from itertools import product
from subprocess import PIPE, Popen
from spacy.tokens import Span
import en_core_web_sm
import networkx as nx
import string
import time
from multiprocessing import Process
from fuzzywuzzy import process
from fuzzywuzzy import fuzz
from itertools import combinations
import numpy as np

from ssmpy import ssm
import ontology_preprocessing
import ontology_preprocessing_test
import ontology_preprocessing_dev

sys.path.append('bin/DiShIn/')  # access bin/DiShIn/

# Input Parameters
sst_light_directory = 'bin/sst-light-0.4/'
temporary_directory = 'temp/'  # temp_sample/
biomedical_entity_to_ontology = {'gene': 'G', 'phenotype': 'H', 'drug': 'C', 'disease': 'D'}

# To exclude (Needs Update)
# neg_gv_list = {'cerubidine', 'trial', '5-fu', 'total', 'multivitamins', 'elemental', 'nitroprusside', 'chlortetracycline', 'transdermal', 'altered', 'promethazine', 'ml', 'fluoroquinolones', 'cephalothin_sodium', 'amiloride', 'tambocor', 'blocking_agents', 'immunosuppressives', 'weight', 'than', 'nabumetone', 'entacapone', 'fexofenadine', 'cytosine_arabinoside', 'drug', 'metaclopramide', 'divalproex_sodium', 'desloratadine', 'database', 'hydantoins', 'benazepril', 'amoxicillin', 'restricted', 'tendency', 'iron_supplements', 'azathioprine', 'exist', 'imidazole', 'half', 'anxiolytics', 'regimen', 'angiotensin-_converting_enzyme_(ace)_inhibitors', 'uroxatral', 'cefoperazone', 'other', 'wshow', 'andusing', '12', 'dobutamine', 'addiction', '500', 'potential', 'lead', 'eliminated', 'transferase', 'leflunomide', 'digitalis_preparations', 'stadol_ns', 'desbutyl_levobupivacaine', 'glibenclamide', 'vinblastine', 'aripiprazole', 'appear', 'oxidase', 'blunt', 'seriraline', 'bedtime', 'arimidex', 'dextromethorphan', 'lanoxin', 'cabergoline', 'oxacillin', 'naprosyn', 'users', 'iloprost', 'local', 'trifluoperazine', 'cefmenoxime', 'plaquenil', 'excess', 'chlorpromazine', 'misused', 'antibiotic', 'involving', 'stanozolol', 'antimycobacterial', 'zdv', 'antidiabetic_products', 'chlorothiazide', 'orlistat', 'bleomycin', 'latamoxef', 'somatostatin_analog', 'slows', 'alternatives', 'make', 'atenolol', 'corresponding', 'seen', 'l50', 'ribavirin', 'dynacirc', 'coumarin_derivatives', 'glyceryl_trinitrate', 'propofol', 'tacrine', 'mepron', 'excreted', 'examining', 'triflupromazine', 'iron_supplement', 'deltasone', 'amlodipine', 'nandrolone', 'antidiabetic_agents', 'antipsychotic_drugs', 'pa', 'containing_compounds', 'er', 'trimethoprim', 'glycoprotein', 'calcitriol', 'multiple', 'angiotensin_ii_receptor_antagonists', 'coa_reductase_inhibitor', 'nonsteroidal_antiinflammatories', 'infused', 'fluvastatin', 'reversible', 'mycophenolate', 'fell', 'vitamin_b3', 'maoi_antidepressants', '-treated', 'aeds', 'induction', 'hypoglycemic_agents', 'antifungal', 'salicylic_acid', 'gabapentin', 'fibrates', 'carvedilol', 'neuromuscular_blocking_agents', 'mesoridazine', 'require', 'fibrinogen', 'predispose', 'anakinra', 'somatostatin_analogs', 'magnesium_hydroxide_antacids', 'pregnancy', ';', 'therefore', 'antiarrhythmic_agents', 'surgery', 'conversion', 'monoamine_oxidase_inhibitor', 'serum', 'cardiac_glycosides', 'fosphenytoin', 'adrenergic_receptor_blockers', 'detected', 'grepafloxacin', 'systemic_corticosteroids', 'nucleoside_reverse_transcriptase_inhibitors', 'divalproex', 'thiocyanate', 'metrizamide', 'included', 'immunosuppressants', 'terbutaline', 'mycophenolate_mofetil', 'modify', 'blocker', 'valsartan', 'sulfoxone', 'distribution', 'famciclovir', 'minutes', 'chelating', 'immunosuppressive_drugs', 'accelerate', 'thrombolytic_agents', 'twice', 'promazine', 'bactrim', 'psychotropic_drugs', 'borne', 'novoseven', 'hivid', 'cromolyn_sodium', 'converting_enzyme_(ace)_inhibitors', 'cleared', 'transport', 'oruvail', 'experience', 'depletion', 'synkayvite', 'chlorthalidone', 'cyp1a2', 'produces', 'hypoglycemia', 'pegasys', 'diagnostic', 'mixing', 'oxc', 'hydroxyurea', 'and/or', 'requiring', 'mtx', 'lithium_carbonate', 'fibric_acid_derivatives', 'rifapentine', 'furafylline', 'dihydropyridine_calcium_antagonists', 'intensified', 'withdrawal', 'ameliorate', 'levonorgestrol', 'rofecoxib', 'ganglionic', 'anaprox', 'hiv_protease_inhibitors', 'studied', 'phenobarbitol', 'threohydrobupropion', 'antithyroid_drugs', 'alg', 'intoxication', 'anagrelide', 'assessed', 'nothiazines', 'terminated', 'coa_reductase_inhibitors', 'ticlopidine', 'cefazolin', 'cyp3a4', 'oxcarbazepine', 'hypokalemia', 'yielded', 'descarboethoxyloratadine', 'oxandrolone', 'leads', 'tranexamic_acid', 'dexmedetomidine', 'pancuronium', 'antacid', 'resorcinol', 'going', 'lenalidomide', 'influence', 'modified', 'pyrantel', 'droperidol', 'replacement', 'benzylpenicillin', 'acting_beta2-agonists', 'n=29', 'sequence', 'utilize', 'gram', 'interferences', 'nicotinic_acid', 'influenced', 'examples', 'min', 'salicylate', 'sulfur', 'keppra', 'iodoquinol', 'hours', 'trimeprazine', 'vitamin_d2', 'tolerated', 'procarbazine', 'volunteers', 'anions', 'increasing', 'etretinate', 'p450', 'nafcillin', 'cyp2c9', 'considered', 'prednisone', 'zofran', 'drawn', 'isradipine', 'lodosyn', 'substrates', 'orencia', 'debrisoquin', 'indicate', 'peginterferon', 'fortified', 'sulfisoxazole', 'tranylcypromine', 'antacid_products', 'antipsychotic_agents', 'antidiabetic_drugs', 'sucralfate', 'hemostasis', 'medrol', 'aminoglutethimide', 'clotrimazole', 'propanolol', 'monotherapy', 'irinotecan', 'identified', '/', 'somatrem', 'acetophenazine', 'gold', 'dirithromycin', 'sympathomimetics', 'erbitux', 'catalyzed', 'indanavir', 'ergonovine', 'lowered', 'infusion', 'combination', 'linezolid', 'substrate', 'differences', 'lowers', 'concomitant', 'nondepolarizing', 'meq', 'sparfloxacin', 'parameters', 'r', 'adjustments', 'prednisolone_sodium_succinate', 'nimodipine', 'tolerance', 'motrin', 'pill', 'sulfadoxine', 'mayuse', 'occurred', 'ci', 'flucoxacillin', 'metoclopramide', 'rifamipicin', 'responsive', 'cycles', 'trials', 'loop_diuretics', 'exhibits', 'folic_acid', 'ceftazidime', 'h2_antagonists', 'lansoprazole', 'escitalopram', 'methylprednisolone', 'antidepressant', 'accounts', 'vitamin_d3', 'gestodene', 'blocking_drugs', 'contribution', 'substances', 'tranylcypromine_sulfate', 'ritanovir', 'nizatidine', 'ingesting', 'buride', 'wthionamide', 'pravastatin', 'gleevec', 'index', 'tikosyn', 'cefotetan', 'antipsychotic_medications', 'aralen', 'performed', 'phenelzine', 'plicamycin', 'possibility', 'betablockers', 'isoenzymes', 'diphenylhydantoin', 'propatenone', 'eproxindine', 'alone', 'determined', 'evaluated', 'profiles', 'bioavailabillty', 'protamine', 'hyperreflexia', 'vitamin_a', 'vitamin_k_antagonists', 'medicine', 'cytokines', 'hydrocodone', 'vs.', 'methylthiotetrazole', 'tested', 'insert', 'antiacid', 'an', 'differ', 'invalidate', 'antiemetics', 'mellaril', 'dosed', 'range', 'bepridil', 'activated_prothrombin_complex_concentrates', 'inactivate', 'exercised', 'etomidate', 'vecuronium', 'coronary_vasodilators', 'dependent', 'anticholinesterases', 'prochlorperazine', 'r-', 'oxymetholone', 'aprepitant', 'ics', 'iressa', 'mephenytoin', 'ramipril', 'novum', 'medication', 'contains', 'diminished', 'activate', 'lam', 'sterilize', 'methandrostenolone', 'antipyrine', 'hydralazine', 'celecoxib', 'hydramine', 'exists', 'antipyretics', 'adenocard', 'besides', 'alpha-', 'cinacalcet', 'demonstrate', 'lomefloxacin', 'cephalothin', 'prolixin', 'concentrates', 'tests', 'analyses', 'proton_pump_inhibitors', 'mean', 'maintained', 'interferon', 'anticholinergic_agents', 'phenformin', 'failed', 'utilization', 'codeine', 'pediapred', 'isosorbide_dinitrate', 'oxaprozin', 'calcium_channel_antagonists', 'magnesium_sulfate', 'nonsteroidal_antiinflammatory_drugs', 'albuterol', 'prazosin', 'replacing', 'expanders', 'showed', 'hypercalcemia', 'benzothiadiazines', 'aza', 'humira', 'aminopyrine', 'cefamandole_naftate', '1/35', 'tolazoline', 'channel_blockers', 'thyroid_hormones', 'orudis', 'selegiline', 'analgesics', 'antagonists', 'ganglionic_blocking_agents', 'antagonism', 'pseudoephedrine', 'calcium_channel_blocking_drugs', 'oxide', 'chemotherapeutic_agents', 'cations', 'tend', 'undergo', 'includes', 'butazone', 'peak', 'sulfonamide', 'enzymes', '%', 'gabitril', 'acarbose', 'simvastatin', 'mixed', 'ethionamide', 'a', 'cyp2d6', 'ergot', 'metabolites', 'interrupted', 'carmustine', 'antianxiety_drugs', 'about', 'decarboxylase_inhibitor', 'hctz', 'advil', 'isosorbide_mononitrate', 'naltrexone', 'experienced', 'niacin', 'potassium_chloride', 'andtolbutamide', 'established', 'streptomycin', 'circulating', 'components', 'induces', 'dihydropyridine_derivative', 'caution', 'clonidine', 'piroxicam', 'phenylpropanolamine', 'label', 'indicated', 'pharmacokinetics', 'im', 'potassium_sparing_diuretics', 'adrenocorticoids', 'ocs', 'penicillin', 'conducted', 'desethylzaleplon', 'felbatol', 'nitrates', 'reviewed', 'smx', 'disease', 'cream', 'control', 'adefovir_dipivoxil', 'ethotoin', 'corticosteroid', 'voltaren', 'antivirals', 'protease_inhibitor', 'furazolidone', 'estrogen', 'investigated', 'mix', 'dapsone_hydroxylamine', 'cefamandole', 'mitotane', 'poisoning', 'metoprolol', 'dopa_decarboxylase_inhibitor', 'incombination', 'nisoldipine', 'diltiazem_hydrochloride', 'adjustment', 'tnf_blocking_agents', 'etodolac', 'phenelzine_sulfate', 'minus', 'formed', 'lower', 'show', 'cardiovasculars', 'sympathomimetic_bronchodilators', 'nitrofurantoin', 'calcium_channel_blocking_agents', 'oxymetazoline', 'neuroleptic', 'tetracyclic_antidepressants', 'steroid_medicine', 'arb', 'phenytoin_sodium', '5-dfur', 'bronchodilators', 'confirmed', 'among', 'sulphenazole', 'antiretroviral_nucleoside_analogues', 'binding', 'imatinib', 'cylates', 'plasmaconcentrations', 'acetohydroxamic_acid', 'inducing'}
neg_gv_list = {}

label_to_pair_type = {'NO_RELATION': 0, 'INHIBITOR': 1, 'PART-OF': 2, 'SUBSTRATE': 3, 'ACTIVATOR': 4,
                      'INDIRECT-DOWNREGULATOR': 5,
                      'ANTAGONIST': 6, 'INDIRECT-UPREGULATOR': 7, 'AGONIST': 8, 'DIRECT-REGULATOR': 9, 'PRODUCT-OF': 10,
                      'AGONIST-ACTIVATOR': 11, 'AGONIST-INHIBITOR': 12, 'SUBSTRATE_PRODUCT-OF': 13}
pair_type_to_label = {v: k for k, v in label_to_pair_type.items()}


# --------------------------------------------------------------
#                 DIVIDED BY SENTENCES ABSTRACTS
# --------------------------------------------------------------


def divided_by_sentences(abstract):
    """Divides abstracts by sentences

    :param abstract:load
    :return: list of sentences of the abstract
    """

    nlp_l = English()
    nlp_l.add_pipe(nlp_l.create_pipe('sentencizer'))
    doc = nlp_l(abstract)
    sentences = [sent.string.strip() for sent in doc.sents]

    return sentences


def get_abstracts_info(abstract_file):
    abstracts_info_dict = {}
    with open(abstract_file, encoding='utf-8') as af:
        for line in af:
            line = line.strip()
            line = line.split('\t')
            # k: abstract_id v: text
            abstracts_info_dict[line[0]] = [line[1]] + divided_by_sentences(line[2])
        print('Abstract file processed')

    return abstracts_info_dict


def get_entities_info(entities_file):
    entities_info_dict = {}
    with open(entities_file, encoding='utf-8') as ef:
        for line in ef:
            line = line.strip()
            line = tuple(line.split('\t'))
            if line[0] in entities_info_dict:
                entities_info_dict[line[0]].append(line[1:])
            else:
                entities_info_dict[line[0]] = [line[1:]]
    print('Entities file processed')

    return entities_info_dict


def get_relations_info(relations_file):
    relations_info_dict = {}
    if relations_file:
        with open(relations_file, encoding='utf-8') as rf:
            for line in rf:
                line = line.strip()
                line = line.split('\t')
                if line[0] in relations_info_dict:
                    relations_info_dict[line[0]].append((line[1], line[2].split(':')[1], line[3].split(':')[1]))
                else:
                    relations_info_dict[line[0]] = [(line[1], line[2].split(':')[1], line[3].split(':')[1])]
            print('Relations file processed')

    return relations_info_dict


if temporary_directory == 'temp/':
    abstracts_info = get_abstracts_info('corpora/training/drugprot_training_abstracs.tsv')
    entities_info = get_entities_info('corpora/training/drugprot_training_entities.tsv')
    relations_info = get_relations_info('corpora/training/drugprot_training_relations.tsv')

elif temporary_directory == 'temp_test/':
    abstracts_info = get_abstracts_info('corpora/test-background/test_background_abstracts.tsv')
    entities_info = get_entities_info('corpora/test-background/test_background_entities.tsv')
    relations_info = get_relations_info('')

elif temporary_directory == 'temp_sample/':
    abstracts_info = get_abstracts_info('container/biocreative_vii/sample_corpora/sample_training_abstracts.tsv')
    entities_info = get_entities_info('container/biocreative_vii/sample_corpora/sample_training_entities.tsv')
    relations_info = get_relations_info('container/biocreative_vii/sample_corpora/sample_training_relations.tsv')

elif temporary_directory == 'temp_dev/':
    abstracts_info = get_abstracts_info('corpora/development/drugprot_development_abstracs.tsv')
    entities_info = get_entities_info('corpora/development/drugprot_development_entities.tsv')
    relations_info = get_relations_info('corpora/development/drugprot_development_relations.tsv')



# --------------------------------------------------------------
#                    UPDATE SENTENCES OFFSETS
# --------------------------------------------------------------

def get_new_offsets_sentences(abstract_id):
    entities_per_sentence = {}
    abstract = abstracts_info[abstract_id]
    annotation_lines = entities_info[abstract_id]

    limit_1 = 0
    limit_2 = 0
    sentence_id = 0

    for sentence in abstract:

        limit_2 += len(sentence)
        entities_per_sentence[('a' + abstract_id + '.s' + str(sentence_id), sentence)] = []

        for annotation in annotation_lines:

            offset_1 = int(annotation[2]) + 1
            offset_2 = int(annotation[3]) + 1

            if limit_1 <= int(offset_1 - 1) <= limit_2 and limit_1 <= int(offset_2 - 1) <= limit_2:
                updated_offset_1 = int(offset_1) - limit_1
                updated_offset_2 = int(offset_2) - limit_1

                entities_per_sentence[('a' + abstract_id + '.s' +
                                       str(sentence_id), sentence)].append(
                    (annotation[0], updated_offset_1, updated_offset_2, annotation[4], annotation[1]))

        sentence_id += 1
        limit_1 += len(sentence) + 1

    return entities_per_sentence


# --------------------------------------------------------------
#                     GET SENTENCE ENTITIES
# --------------------------------------------------------------

def get_sentence_entities():
    """

    :return:
    """
    entities = {}
    positive_entities = set()

    if temporary_directory == 'temp_test/':
        all_data = ontology_preprocessing_test.load_data()
    elif temporary_directory == 'temp_dev/':
        all_data = ontology_preprocessing_dev.load_data()
    else:
        all_data = ontology_preprocessing.load_data()

    for k in entities_info.keys():
        entities_per_sentence = get_new_offsets_sentences(k)

        for sentence, entities_sentence in entities_per_sentence.items():

            entities_sentence = sorted(entities_sentence, key=lambda x: x[1])
            sentence_entities = {}

            entity_number = 1
            for entity in entities_sentence:

                if k in relations_info:
                    for pair in relations_info[k]:

                        if entity[0] == pair[1]:
                            positive_entities.add(sentence[0] + '.u' + entity[0] + '.e' + str(entity_number))

                        elif entity[0] == pair[2]:
                            positive_entities.add(sentence[0] + '.u' + entity[0] + '.e' + str(entity_number))

                sentence_entities[sentence[0] + '.u' + entity[0] + '.e' + str(entity_number)] = (
                eval('[' + str(entity[1]) + ', ' + str(entity[2]) + ']'), entity[3], all_data[entity[3]])
                entity_number += 1

            entities[sentence[0]] = sentence_entities

    return entities, positive_entities


# --------------------------------------------------------------
#                     GET ENTITIES ANCESTORS
# --------------------------------------------------------------

def get_common_ancestors(id1, id2):
    """

    :param id1:
    :param id2:
    :return:
    """

    if id1.startswith('CHEBI'):
        ssm.semantic_base('bin/DiShIn/chebi.db')
    if id1.startswith('GO'):
        ssm.semantic_base('bin/DiShIn/go.db')

    e1 = ssm.get_id(id1.replace(':', '_'))

    if id2.startswith('CHEBI'):
        ssm.semantic_base('bin/DiShIn/chebi.db')
    if id2.startswith('GO'):
        ssm.semantic_base('bin/DiShIn/go.db')

    e2 = ssm.get_id(id2.replace(':', '_'))

    a = ssm.common_ancestors(e1, e2)
    a = [ssm.get_name(x) for x in a]

    return a


def get_path_to_root(entity_id):
    """

    :param entity_id:
    :return:
    """

    if entity_id.startswith('CHEBI'):
        ssm.semantic_base('bin/DiShIn/chebi.db')
    if entity_id.startswith('GO'):
        ssm.semantic_base('bin/DiShIn/go.db')

    e1 = ssm.get_id(entity_id.replace(':', '_'))

    a = ssm.common_ancestors(e1, e1)
    a = [ssm.get_name(x) for x in a]

    return a


def get_ancestors(sentence_labels, sentence_entities):
    """Obtain the path to lowest common ancestor of each entity of each pair and path from LCA to root

    :param sentence_labels: list of (e1, e2)
    :param sentence_entities: dictionary mapping entity ID to ((e_start, e_end), text, paths_to_root)
    :return: left and right paths to LCA
    """

    right_paths = []
    left_paths = []
    common_ancestors = []

    for p in sentence_labels:
        instance_ancestors = get_common_ancestors(sentence_entities[p[0]][2], sentence_entities[p[1]][2])
        left_path = get_path_to_root(sentence_entities[p[0]][2])
        right_path = get_path_to_root(sentence_entities[p[1]][2])

        # print('Common ancestors:', sentence_entities[p[0]][1:], sentence_entities[p[1]][1:], instance_ancestors)

        instance_ancestors = [i for i in instance_ancestors if
                              i.startswith('CHEBI') or i.startswith('HP') or i.startswith('GO') or i.startswith('DOID')]
        left_path = [i for i in left_path if
                     i.startswith('CHEBI') or i.startswith('HP') or i.startswith('GO') or i.startswith('DOID')]
        right_path = [i for i in right_path if
                      i.startswith('CHEBI') or i.startswith('HP') or i.startswith('GO') or i.startswith('DOID')]

        common_ancestors.append(instance_ancestors)
        left_paths.append(left_path)
        right_paths.append(right_path)

    return common_ancestors, (left_paths, right_paths)


# --------------------------------------------------------------
#               PARSE CORPUS SENTENCES USING SPACY
# --------------------------------------------------------------

def prevent_sentence_segmentation(doc):
    """

    :param doc:
    :return:
    """

    for token in doc:
        # This will entirely disable spaCy's sentence detection
        token.is_sent_start = False

    return doc


nlp = en_core_web_sm.load(disable=['ner'])
nlp.add_pipe(prevent_sentence_segmentation, name='prevent-sbd', before='parser')


def parse_sentence_spacy(sentence_text, sentence_entities):
    """

    :param sentence_text:
    :param sentence_entities:
    :return:
    """

    # Use spacy to parse a sentence
    for e in sentence_entities:
        idx = sentence_entities[e][0]
        sentence_text = sentence_text[:idx[0] - 1] + sentence_text[idx[0] - 1:idx[1] - 1].replace(' ', '_').replace('-', '_')\
            .replace(':', '_') + sentence_text[idx[1] - 1:]

    # Clean text to make tokenization easier
    sentence_text = sentence_text.replace(';', ',')
    sentence_text = sentence_text.replace('*', ' ')
    sentence_text = sentence_text.replace(':', ',')
    sentence_text = sentence_text.replace('-', ' ')
    sentence_text = sentence_text.replace(']', ' ')
    sentence_text = sentence_text.replace('[', ' ').replace('\'', ' ').replace('/', ' ')

    parsed = nlp(sentence_text)

    return parsed


def run_sst(base_dir, token_seq):
    """

    :param base_dir:
    :param token_seq:
    :return:
    """

    chunk_size = 500
    wordnet_tags = {}
    sent_ids = list(token_seq.keys())

    chunks = [sent_ids[i:i + chunk_size] for i in range(0, len(sent_ids), chunk_size)]

    for i, chunk in enumerate(chunks):
        sentence_file = open('{}/sentences_{}.txt'.format(temporary_directory + base_dir.split('/')[1], i), 'w',
                             encoding='utf-8')

        for sent in chunk:
            sentence_file.write("{}\t{}\t.\n".format(sent, '\t'.join(token_seq[sent])))

        sentence_file.close()
        sst_args = [sst_light_directory + 'sst', 'bitag',
                    '{}/MODELS/WSJPOSc_base_20'.format(sst_light_directory),
                    '{}/DATA/WSJPOSc.TAGSET'.format(sst_light_directory),
                    '{}/MODELS/SEM07_base_12'.format(sst_light_directory),
                    '{}/DATA/WNSS_07.TAGSET'.format(sst_light_directory),
                    '{}/sentences_{}.txt'.format(temporary_directory + base_dir.split('/')[1], i), '0', '0']

        p = Popen(sst_args, stdout=PIPE)
        p.communicate()

        with open('{}/sentences_{}.txt.tags'.format(temporary_directory + base_dir.split('/')[1], i),
                  encoding='utf-8') as f:
            output = f.read()

        sstoutput = parse_sst_results(output)
        wordnet_tags.update(sstoutput)

    return wordnet_tags


def parse_sst_results(results):
    """

    :param results:
    :return:
    """

    sentences = {}
    lines = results.strip().split('\n')

    for l in lines:
        values = l.split('\t')
        wntags = [x.split(' ')[-1].split('-')[-1] for x in values[1:]]
        sentences[values[0]] = wntags

    return sentences


def parse_sentences_spacy(base_dir, entities):
    """

    :param base_dir:
    :param entities:
    :return:
    """

    parsed_sentences = {}

    # First iterate all documents, and preprocess all sentences
    token_seq = {}

    for k in entities_info.keys():

        entities_per_sentence = get_new_offsets_sentences(k)

        for sentence, entities_sentence in entities_per_sentence.items():

            parsed_sentence = parse_sentence_spacy(sentence[1], entities[sentence[0]])
            parsed_sentences[sentence[0]] = parsed_sentence

            tokens = []

            for t in parsed_sentence:
                tokens.append(t.text.replace(' ', '_').replace('\t', '_').replace('\n', '_'))

            token_seq[sentence[0]] = tokens

    wordnet_tags = run_sst(base_dir, token_seq)

    return parsed_sentences, wordnet_tags


def get_network_graph_spacy(document):
    """Convert the dependencies of the spacy document object to a networkX graph

    :param document: spacy parsed document object
    :return: networkX graph object and nodes list
    """

    edges = []
    nodes = []

    # Ensure that every token is connected
    for s in document.sents:
        edges.append(('ROOT', '{0}-{1}'.format(s.root.lower_, s.root.i)))

    for token in document:
        nodes.append('{0}-{1}'.format(token.lower_, token.i))

        for child in token.children:
            edges.append(('{0}-{1}'.format(token.lower_, token.i), '{0}-{1}'.format(child.lower_, child.i)))

    return nx.Graph(edges), nodes


def get_head_tokens_spacy(entities, sentence, positive_entities):
    """

    :param entities: dictionary mapping entity IDs to (offset, text)
    :param sentence: sentence parsed by spacy
    :param positive_entities:
    :return: dictionary mapping head tokens word-idx to entity IDs
    """

    sentence_head_tokens_type_1 = {}
    sentence_head_tokens_type_2 = {}
    pos_gv = set()
    neg_gv = set()

    for eid in entities:

        offset = (entities[eid][0][0], entities[eid][0][-1])
        entity_tokens = sentence.char_span(offset[0] - 1, offset[1] - 1)

        # i = 1
        # while not entity_tokens and i + offset[1] < len(sentence.text) + 1:
        #     entity_tokens = sentence.char_span(offset[0], offset[1] + i)
        #     i += 1
        #
        # i = 0
        # while not entity_tokens and offset[0] - i > 0:
        #     entity_tokens = sentence.char_span(offset[0] - i, offset[1])
        #     i += 1

        # print('-------------------------')
        # if not entity_tokens:
        #     sentence_text = sentence.text.replace('_', ' ')
        #     #print(entities[eid][1])
        #     if entities[eid][1] == sentence_text[offset[0] - 1:offset[1] - 1]:
        #         try:
        #             #print(sentence_text[offset[0] - 1:offset[1] - 1])
        #             print(sentence.text)
        #             print(offset[0] - 1)
        #             print(offset[1] - 1)
        #             entity_tokens = Span(sentence, offset[0] - 1, offset[1] - 1)
        #             print(entity_tokens)
        #         except IndexError:
        #             #print(sentence_text)
        #             #print(sentence.text)
        #             print()

        if not entity_tokens:
            logging.warning(('No tokens found:', entities[eid], sentence.text, '|'.join([t.text for t in sentence])))

        else:

            head_token = '{0}-{1}'.format(entity_tokens.root.lower_.replace(' ', '_'), entity_tokens.root.i)

            if eid in positive_entities:
                pos_gv.add(entity_tokens.root.head.lower_)
            else:
                neg_gv.add(entity_tokens.root.head.lower_)

            if head_token in sentence_head_tokens_type_1:
                logging.warning(('Head token conflict:', sentence_head_tokens_type_1[head_token], entities[eid]))
            elif head_token in sentence_head_tokens_type_2:
                logging.warning(('Head token conflict:', sentence_head_tokens_type_2[head_token], entities[eid]))

            if biomedical_entity_to_ontology['drug'] == entities[eid][2][0]:
                sentence_head_tokens_type_1[head_token] = eid
            elif biomedical_entity_to_ontology['gene'] == entities[eid][2][0]:
                sentence_head_tokens_type_2[head_token] = eid

    return sentence_head_tokens_type_1, sentence_head_tokens_type_2, pos_gv, neg_gv


def process_sentence_spacy(base_dir, sentence, sentence_entities, sentence_pairs, positive_entities, wordnet_tags=None,
                           mask_entities=True, min_sdp_len=0, max_sdp_len=15):
    """Process sentence to obtain labels, instances and classes for a ML classifier

    :param base_dir:
    :param sentence: sentence processed by spacy
    :param sentence_entities: dictionary mapping entity ID to ((e_start, e_end), text, paths_to_root)
    :param sentence_pairs: dictionary mapping pairs of known entities in this sentence to pair types
    :param positive_entities:
    :param wordnet_tags:
    :param mask_entities:
    :param min_sdp_len:
    :param max_sdp_len:
    :return: labels of each pair (according to sentence_entities, word vectors and classes (pair types according to sentence_pairs)
    """

    left_word_vectors = []
    right_word_vectors = []
    left_wordnets = []
    right_wordnets = []
    classes = []
    labels = []

    graph, nodes_list = get_network_graph_spacy(sentence)
    sentence_head_tokens_type_1, sentence_head_tokens_type_2, pos_gv, neg_gv = get_head_tokens_spacy(sentence_entities,
                                                                                                     sentence,
                                                                                                     positive_entities)
    entity_offsets = [sentence_entities[x][0][0] for x in sentence_entities]

    for (e1, e2) in product(sentence_head_tokens_type_1, sentence_head_tokens_type_2):

        if sentence_head_tokens_type_1.get(e1):
            if int(sentence_head_tokens_type_1[e1].split('e')[-1]) > int(
                    sentence_head_tokens_type_2[e2].split('e')[-1]):
                e1, e2 = e2, e1
        else:
            if int(sentence_head_tokens_type_2[e1].split('e')[-1]) > int(
                    sentence_head_tokens_type_1[e2].split('e')[-1]):
                e2, e1 = e1, e2

        try:
            if sentence_head_tokens_type_1.get(e1):
                e1_text = sentence_entities[sentence_head_tokens_type_1[e1]]
                e2_text = sentence_entities[sentence_head_tokens_type_2[e2]]

            else:
                e1_text = sentence_entities[sentence_head_tokens_type_1[e2]]
                e2_text = sentence_entities[sentence_head_tokens_type_2[e1]]

        except KeyError:
            print(e1, e2)
            print(sentence)
            print(sentence_head_tokens_type_1)
            print(sentence_head_tokens_type_2)

        if e1_text[1].lower() == e2_text[1].lower():
            continue

        middle_text = sentence.text[e1_text[0][-1]:e2_text[0][0]]

        if middle_text.strip() in string.punctuation:
            continue

        # if 'train' in base_dir:
        #     middle_text = sentence.text[e1_text[0][-1]:e2_text[0][0]]
        #
        #     if middle_text.strip() in string.punctuation:
        #         continue

        try:
            sdp = nx.shortest_path(graph, source=e1, target=e2)

            if len(sdp) < min_sdp_len or len(sdp) > max_sdp_len:
                continue

            neg = False
            is_neg_gv = False
            for i, element in enumerate(sdp):
                token_text = element.split('-')[0]
                if (i == 1 or i == len(sdp) - 2) and token_text in neg_gv_list:
                    pass
                    print('Skipped gv {} {}:'.format(token_text, str(sdp)))

            if neg or is_neg_gv:
                continue

            vector = []
            wordnet_vector = []
            negations = 0
            head_token_position = None

            for i, element in enumerate(sdp):
                if element != 'ROOT':
                    token_idx = int(element.split('-')[-1])  # get the index of the token
                    sdp_token = sentence[token_idx]  # get the token object

                    if mask_entities and sdp_token.idx in entity_offsets:
                        vector.append('entity')
                    else:
                        vector.append(sdp_token.text)
                    if wordnet_tags:
                        wordnet_vector.append(wordnet_tags[token_idx])

                    head_token = '{}-{}'.format(sdp_token.head.lower_, sdp_token.head.i)  # get the key of head token

                    # Head token must not have its head in the path, otherwise that would be the head token
                    # In some cases the token is its own head
                    if head_token not in sdp or head_token == element:
                        head_token_position = i + negations

            if head_token_position is None:
                print('Head token not found:', e1_text, e2_text, sdp)
                sys.exit()
            else:
                left_vector = vector[:head_token_position + 1]
                right_vector = vector[head_token_position:]
                left_wordnet = wordnet_vector[:head_token_position + 1]
                right_wordnet = wordnet_vector[head_token_position:]

            left_word_vectors.append(left_vector)
            right_word_vectors.append(right_vector)
            left_wordnets.append(left_wordnet)
            right_wordnets.append(right_wordnet)

        except nx.exception.NetworkXNoPath:
            logging.warning('No path:', e1_text, e2_text, graph.nodes())
            left_word_vectors.append([])
            right_word_vectors.append([])
            left_wordnets.append([])
            right_wordnets.append([])

        except nx.NodeNotFound:
            logging.warning(('Node not found:', e1_text, e2_text, e1, e2, list(sentence), graph.nodes()))
            left_word_vectors.append([])
            right_word_vectors.append([])
            left_wordnets.append([])
            right_wordnets.append([])

        if sentence_head_tokens_type_1.get(e1):
            labels.append((sentence_head_tokens_type_1[e1], sentence_head_tokens_type_2[e2]))
            if (sentence_head_tokens_type_1[e1], sentence_head_tokens_type_2[e2]) in sentence_pairs:
                classes.append(sentence_pairs[(sentence_head_tokens_type_1[e1], sentence_head_tokens_type_2[e2])])
            else:
                classes.append(0)
        else:
            labels.append((sentence_head_tokens_type_1[e2], sentence_head_tokens_type_2[e1]))
            if (sentence_head_tokens_type_1[e2], sentence_head_tokens_type_2[e1]) in sentence_pairs:
                classes.append(sentence_pairs[(sentence_head_tokens_type_1[e2], sentence_head_tokens_type_2[e1])])
            else:
                classes.append(0)

    return labels, (left_word_vectors, right_word_vectors), (left_wordnets, right_wordnets), classes, pos_gv, neg_gv


# --------------------------------------------------------------
#    PARSE CORPUS WITH SDP VECTORS FOR EACH RELATION INSTANCE
# --------------------------------------------------------------

# def get_sdp_instances(base_dir, parser='spacy'):
#     """Parse corpus, return vectors of SDP of each relation instance
#
#     :param base_dir: directory containing the sentences
#     :param parser:
#     :return: labels (eid1, eid2), instances (vectors), classes (0/1), common ancestors, l/r ancestors, l/r wordnet
#     """
#
#     s = time.time()
#     entities, positive_entities = get_sentence_entities()
#
#     if parser == 'spacy':
#         parsed_sentences, wordnet_sentences = parse_sentences_spacy(base_dir, entities)
#
#     else:
#         parsed_sentences = wordnet_sentences = None
#
#     left_instances = []
#     right_instances = []
#     left_ancestors = []
#     right_ancestors = []
#     common_ancestors = []
#     left_wordnet = []
#     right_wordnet = []
#     classes = []
#     labels = []
#     all_pos_gv = set()
#     all_neg_gv = set()
#
#     for k in entities_info.keys():
#             entities_per_sentence = get_new_offsets_sentences(k)
#
#             for sentence, entities_sentence in entities_per_sentence.items():
#
#                 sentence_pairs = {}
#                 sentence_entities = entities_sentence
#                 #entity id; offset1; offset2; entity text; entity type.
#
#                 all_pairs = {}
#                 if k in relations_info:
#                     for pair in relations_info[k]:
#                         if pair[1] in all_pairs:
#                             all_pairs[pair[1]].append((pair[0], pair[2]))
#                         else:
#                             all_pairs[pair[1]] = []
#                             all_pairs[pair[1]].append((pair[0], pair[2]))
#
#                 comb = combinations([e[0] for e in sentence_entities], 2)
#
#                 for pair in comb:
#                     sentence_pairs[(sentence[0] + '.e' + str(pair[0].replace('T', '')), sentence[0] + '.e' + str(pair[1].replace('T', '')))] = 0
#
#                 for entity in sentence_entities:
#                     if entity[0] in all_pairs:
#                         for matching_entities in sentence_entities:
#                             for pair in all_pairs[entity[0]]:
#                                 if pair[1] == matching_entities[0]:
#                                     sentence_pairs[(sentence[0] + '.e' + str(entity[0].replace('T', '')), sentence[0] + '.e' + str(matching_entities[0].replace('T', '')))] = label_to_pair_type[pair[0]]
#
#                 if len(sentence_pairs) > 0:  # skip sentences without pairs
#                     sentence_entities = entities[sentence[0]]
#                     parsed_sentence = parsed_sentences[sentence[0]]
#                     wordnet_sentence = wordnet_sentences[sentence[0]]
#
#                     if parser == 'spacy':
#
#                         sentence_labels, sentence_we_instances, sentence_wn_instances, sentence_classes, pos_gv, neg_gv = \
#                             process_sentence_spacy(base_dir, parsed_sentence, sentence_entities, sentence_pairs, positive_entities, wordnet_sentence)
#
#                     else:
#                         sentence_labels = sentence_we_instances = sentence_classes = sentence_wn_instances = pos_gv = neg_gv = None
#
#                     sentence_ancestors, sentence_subpaths = get_ancestors(sentence_labels, sentence_entities)
#
#                     labels += sentence_labels
#                     left_instances += sentence_we_instances[0]
#                     right_instances += sentence_we_instances[1]
#                     classes += sentence_classes
#                     common_ancestors += sentence_ancestors
#                     left_ancestors += sentence_subpaths[0]
#                     right_ancestors += sentence_subpaths[1]
#
#                     left_wordnet += sentence_wn_instances[0]
#                     right_wordnet += sentence_wn_instances[1]
#
#                     all_pos_gv.update(pos_gv)
#                     all_neg_gv.update(neg_gv)
#     e = time.time()
#     print(e-s)
#     return labels, (left_instances, right_instances), classes, common_ancestors, (left_ancestors, right_ancestors), (left_wordnet, right_wordnet), all_neg_gv, all_pos_gv
#


#### WORK ZONE ##########################

def get_sdp_instances(base_dir, parser='spacy'):
    """Parse corpus, return vectors of SDP of each relation instance

    :param base_dir: directory containing the sentences
    :param parser:
    :return: labels (eid1, eid2), instances (vectors), classes (0/1), common ancestors, l/r ancestors, l/r wordnet
    """

    s = time.time()
    entities, positive_entities = get_sentence_entities()

    if parser == 'spacy':
        parsed_sentences, wordnet_sentences = parse_sentences_spacy(base_dir, entities)

    else:
        parsed_sentences = wordnet_sentences = None

    left_instances = []
    right_instances = []
    left_ancestors = []
    right_ancestors = []
    common_ancestors = []
    left_wordnet = []
    right_wordnet = []
    classes = []
    labels = []
    all_pos_gv = set()
    all_neg_gv = set()

    s = time.time()
    counter = 0

    print(len(entities_info.keys()))
    for k in entities_info.keys():  # 3500 train 10750 test
        entities_per_sentence = get_new_offsets_sentences(k)

        q = multiprocessing.Queue()
        processes = []

        for sentence, entities_sentence in entities_per_sentence.items():
            p = Process(target=add_helper, args=(
            q, base_dir, parser, entities, positive_entities, parsed_sentences, wordnet_sentences, k, sentence, entities_sentence,))
            processes.append(p)
            p.start()

        for proc in processes:
            try:
                out = q.get(True, 5)

                if out[0]:
                    labels += out[0]
                else:
                    labels += []

                if out[1]:
                    left_instances += out[1]
                else:
                    left_instances += []

                if out[2]:
                    right_instances += out[2]
                else:
                    right_instances += []

                if out[3]:
                    classes += out[3]
                else:
                    classes += []

                if out[4]:
                    common_ancestors += out[4]
                else:
                    common_ancestors += []

                if out[5]:
                    left_ancestors += out[5]
                else:
                    left_ancestors += []

                if out[6]:
                    right_ancestors += out[6]
                else:
                    right_ancestors += []

                if out[7]:
                    left_wordnet += out[7]
                else:
                    left_wordnet += []

                if out[8]:
                    right_wordnet += out[8]
                else:
                    right_wordnet += []

                if out[9]:
                    all_pos_gv.update(out[9])
                else:
                    all_pos_gv.update()

                if out[10]:
                    all_neg_gv.update(out[10])
                else:
                    all_neg_gv.update()

            except Exception:
                break

        for proc in processes:
            p.join(2)
            # try:
            #     p.join(2)
            # except Exception:
            #     break

        # TO BE SAFE
        counter += 1
        if temporary_directory == 'temp/':
            divider = 500
        elif temporary_directory == 'temp_dev/':
            divider = 250
        else:
            divider = 1075

        if counter % divider == 0:
            print('\n\n')
            print(counter)
            print('\n\n')
            train_labels, x_train, y_train, x_train_ancestors, x_train_subpaths, x_train_wordnet = labels, (
            left_instances, right_instances), classes, common_ancestors, (left_ancestors, right_ancestors), (
                                                                                                   left_wordnet,
                                                                                                   right_wordnet)

            if temporary_directory == 'temp/':
                preprocess_what = 'train'
            elif temporary_directory == 'temp_test/':
                preprocess_what = 'test'
            elif temporary_directory == 'temp_sample/':
                preprocess_what = 'sample'
            else:
                preprocess_what = 'dev'    

            np.save(temporary_directory + 'drug_gene/' + preprocess_what + str(counter) + '_y_ck.npy', y_train)
            np.save(temporary_directory + 'drug_gene/' + preprocess_what + str(counter) + '_labels_ck.npy',
                    train_labels)
            np.save(temporary_directory + 'drug_gene/' + preprocess_what + str(counter) + '_x_words_ck.npy', x_train)
            np.save(temporary_directory + 'drug_gene/' + preprocess_what + str(counter) + '_x_wordnet_ck.npy',
                    x_train_wordnet)
            np.save(temporary_directory + 'drug_gene/' + preprocess_what + str(counter) + '_x_subpaths_ck.npy',
                    x_train_subpaths)
            np.save(temporary_directory + 'drug_gene/' + preprocess_what + str(counter) + "_x_ancestors_ck.npy",
                    x_train_ancestors)

    e = time.time()
    print(e - s)

    return labels, (left_instances, right_instances), classes, common_ancestors, (left_ancestors, right_ancestors), (
    left_wordnet, right_wordnet), all_neg_gv, all_pos_gv


def add_helper(q, base_dir, parser, entities, positive_entities, parsed_sentences, wordnet_sentences, k, sentence,
               entities_sentence):
    sentence_pairs = {}
    sentence_entities = entities_sentence
    # entity id; offset1; offset2; entity text; entity type.

    sentence_entities = sorted(sentence_entities, key=lambda x: x[1])

    support_dict = {}
    counter = 1
    for element in sentence_entities:
        support_dict[element[0]] = str(counter)
        counter += 1

    all_pairs = {}
    if k in relations_info:
        for pair in relations_info[k]:
            if pair[1] in all_pairs:
                all_pairs[pair[1]].append((pair[0], pair[2]))
            else:
                all_pairs[pair[1]] = []
                all_pairs[pair[1]].append((pair[0], pair[2]))

    comb = combinations([e[0] for e in sentence_entities], 2)

    for pair in comb:
        sentence_pairs[(sentence[0] + '.u' + pair[0] + '.e' + support_dict[pair[0]],
                        sentence[0] + '.u' + pair[1] + '.e' + support_dict[pair[1]])] = label_to_pair_type[
            'NO_RELATION']

    for entity in sentence_entities:
        if entity[0] in all_pairs:
            for matching_entities in sentence_entities:
                for pair in all_pairs[entity[0]]:
                    if pair[1] == matching_entities[0]:
                        sentence_pairs[(sentence[0] + '.u' + entity[0] + '.e' + support_dict[entity[0]],
                                        sentence[0] + '.u' + matching_entities[0] + '.e' + support_dict[
                                            matching_entities[0]])] = \
                            label_to_pair_type[pair[0]]

    if len(sentence_pairs) > 0:  # skip sentences without pairs
        sentence_entities = entities[sentence[0]]
        parsed_sentence = parsed_sentences[sentence[0]]
        wordnet_sentence = wordnet_sentences[sentence[0]]

        if parser == 'spacy':

            sentence_labels, sentence_we_instances, sentence_wn_instances, sentence_classes, pos_gv, neg_gv = \
                process_sentence_spacy(base_dir, parsed_sentence, sentence_entities, sentence_pairs,
                                       positive_entities, wordnet_sentence)

        else:
            sentence_labels = sentence_we_instances = sentence_classes = sentence_wn_instances = pos_gv = neg_gv = None

        sentence_ancestors, sentence_subpaths = get_ancestors(sentence_labels, sentence_entities)

        labels = sentence_labels
        left_instances = sentence_we_instances[0]
        right_instances = sentence_we_instances[1]
        classes = sentence_classes
        common_ancestors = sentence_ancestors
        left_ancestors = sentence_subpaths[0]
        right_ancestors = sentence_subpaths[1]

        left_wordnet = sentence_wn_instances[0]
        right_wordnet = sentence_wn_instances[1]

        q.put([labels, left_instances, right_instances, classes, common_ancestors, left_ancestors, right_ancestors,
               left_wordnet, right_wordnet, pos_gv, neg_gv])

    return

