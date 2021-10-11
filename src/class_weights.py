import numpy as np
import math

# labels_dict : {ind_label: count_label}
# mu : parameter to tune


def create_class_weight(labels_dict,mu=0.15):
    total = np.sum(list(labels_dict.values()))
    keys = labels_dict.keys()
    class_weight = dict()

    for key in keys:
        score = math.log(mu*total/float(labels_dict[key]))
        class_weight[key] = score if score > 1.0 else 1.0

    return class_weight


# random labels_dict
labels_dict = {0: 17288, 1: 5392, 2: 886, 3: 2003, 4: 1429, 5: 1330, 6: 972, 7: 1379, 8: 659, 9: 2250, 10: 921, 11: 29, 12: 13, 13: 25}

print(create_class_weight(labels_dict))


# label_to_pair_type = {'NO_RELATION': 0,
#                       'INHIBITOR': 1,
#                       'PART-OF': 2,
#                       'SUBSTRATE': 3,
#                       'ACTIVATOR': 4,
#                       'INDIRECT-DOWNREGULATOR': 5,
#                       'ANTAGONIST': 6,
#                       'INDIRECT-UPREGULATOR': 7,
#                       'AGONIST': 8,
#                       'DIRECT-REGULATOR': 9,
#                       'PRODUCT-OF': 10,
#                       'AGONIST-ACTIVATOR': 11,
#                       'AGONIST-INHIBITOR': 12,
#                       'SUBSTRATE_PRODUCT-OF': 13}
