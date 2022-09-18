import copy
import itertools
from io import StringIO
import math


class Queue(object):

    def __init__(self):
        self.queue = []

    def enqueue(self, b):
        self.queue.insert(0, b)

    def dequeue(self):
        return self.queue.pop()

    def isEmpty(self):
        return self.queue == []


def getheight(node):
    if not node:
        return 0
    else:
        return max(getheight(node.left), getheight(node.right)) + 1


def add_padding(str, pad_length_value):
    str = str.strip()
    return str.center(pad_length_value, ' ')


# sotre node , space and slashes in list first, then print out
def pretty_print(tree):
    output = StringIO()
    pretty_output = StringIO()

    current_level = Queue()
    next_level = Queue()
    current_level.enqueue(tree)
    depth = 0

    # get the depth of current tree
    # get the tree node data and store in list
    if tree:
        while not current_level.isEmpty():
            current_node = current_level.dequeue()
            output.write('%s ' % current_node.data if current_node else 'N ')
            next_level.enqueue(
                current_node.negative_child if current_node else current_node)
            next_level.enqueue(
                current_node.positive_child if current_node else current_node)

            if current_level.isEmpty():
                if sum([i is not None for i in next_level.queue]
                       ):  # if next level has node
                    current_level, next_level = next_level, current_level
                    depth = depth + 1
                output.write('\n')
    print('the tree printed level by level is :')
    print(output.getvalue())
    print("current tree's depth is %i" % (depth + 1))

    # add space to each node
    output.seek(0)
    pad_length = 3
    keys = []
    spaces = int(math.pow(2, depth))

    while spaces > 0:
        skip_start = spaces * pad_length
        skip_mid = (2 * spaces - 1) * pad_length

        key_start_spacing = ' ' * skip_start
        key_mid_spacing = ' ' * skip_mid

        keys = output.readline().split(' ')  # read one level to parse
        padded_keys = (add_padding(key, pad_length) for key in keys)
        padded_str = key_mid_spacing.join(padded_keys)
        complete_str = ''.join([key_start_spacing, padded_str])

        pretty_output.write(complete_str)

        # add space and slashes to middle layer
        slashes_depth = spaces
        print('current slashes depth is:')
        print(spaces)
        print("current levle's list is:")
        print(keys)
        spaces = spaces // 2
        if spaces > 0:
            pretty_output.write('\n')  # print '\n' each level

            cnt = 0
            while cnt < slashes_depth:
                inter_symbol_spacing = ' ' * (pad_length + 2 * cnt)
                symbol = ''.join(['/', inter_symbol_spacing, '\\'])
                symbol_start_spacing = ' ' * (skip_start - cnt - 1)
                symbol_mid_spacing = ' ' * (skip_mid - 2 * (cnt + 1))
                pretty_output.write(''.join([symbol_start_spacing, symbol]))
                for i in keys[1:-1]:
                    pretty_output.write(''.join([symbol_mid_spacing, symbol]))
                pretty_output.write('\n')
                cnt = cnt + 1

    print(pretty_output.getvalue())

class Node:
    def __init__(self, data, positive_child=None, negative_child=None):
        self.data = data
        self.positive_child = positive_child
        self.negative_child = negative_child


class Record:
    def __init__(self, illness, symptoms):
        self.illness = illness
        self.symptoms = symptoms


def parse_data(filepath):
    with open(filepath) as data_file:
        records = []
        for line in data_file:
            words = line.strip().split()
            records.append(Record(words[0], words[1:]))
        return records


class Diagnoser:
    def __init__(self, root):
        self.root = root

    def diagnos_helper(self, node, symptoms):
        if node.positive_child is None:
            return node.data

        if node.data in symptoms:
            return self.diagnos_helper(node.positive_child, symptoms)
        else:
            return self.diagnos_helper(node.negative_child, symptoms)

    def diagnose(self, symptoms):
        return self.diagnos_helper(self.root, symptoms)

    def calculate_success_rate(self, records):
        if len(records)== 0:
            return 0
        successful_diagnoses = 0
        for record in records:
            if self.diagnose(record.symptoms) == record.illness:
                successful_diagnoses += 1
        return float(successful_diagnoses / len(records))

    def all_illnesses(self):
        illnes_dict = {}
        self.all_illnesses_helper(self.root, illnes_dict)
        print(illnes_dict,"illnes_dict")
        return [k for k in sorted(illnes_dict, key=illnes_dict.get, reverse=True)]

    def all_illnesses_helper(self, node, illnes_dict):
        if node.positive_child is None:
            if node.data is not None:
                if node.data in illnes_dict:
                    illnes_dict[node.data] += 1
                else:
                    illnes_dict[node.data] = 1
        else:
            self.all_illnesses_helper(node.positive_child, illnes_dict)
            self.all_illnesses_helper(node.negative_child, illnes_dict)

    def paths_to_illness(self, illness):
        path_list = []
        self.paths_to_illness_helper(self.root, illness, path_list, [])
        return path_list

    def paths_to_illness_helper(self, node, illness, path_list, path):
        if node.positive_child is None:
            if node.data == illness:
                path_list.append(path)

        else:
            pathh = copy.copy(path)
            path.append(True)
            pathh.append(False)
            self.paths_to_illness_helper(node.positive_child, illness, path_list, path)
            self.paths_to_illness_helper(node.negative_child, illness, path_list, pathh)


def build_tree(records, symptoms):
    root = Node("_")

    # build tree with None leaves
    build_tree_helper_structure(root, symptoms)

    # fills illness in leaves
    build_tree_helper_disease_replace(root, [], records, symptoms)
    return root


def build_tree_helper_disease_replace(node, path, records, symptoms):
    if node.negative_child is None:
        illnes_ratte_dict = {}
        for record in records:
            for i in range(len(symptoms)):
                if (symptoms[i] in record.symptoms and symptoms[i] not in path) or (
                        symptoms[i] in path and symptoms[i] not in record.symptoms):
                    break
                if i == len(symptoms) - 1:
                    if record.illness in illnes_ratte_dict:
                        illnes_ratte_dict[record.illness] += 1
                    else:
                        illnes_ratte_dict[record.illness] = 1

        if len(illnes_ratte_dict) != 0:
            node.data = max(illnes_ratte_dict, key=lambda k: illnes_ratte_dict[k])
        else:
            node.data = None
    else:
        pathh = copy.copy(path)
        path.append(node.data)
        build_tree_helper_disease_replace(node.positive_child, path, records, symptoms)
        build_tree_helper_disease_replace(node.negative_child, pathh, records, symptoms)


def build_tree_helper_structure(node, symptoms):
    if len(symptoms) != 0:
        node.data = symptoms[0]
        node.positive_child = Node("_")
        node.negative_child = Node("_")
        build_tree_helper_structure(node.positive_child, symptoms[1:])
        build_tree_helper_structure(node.negative_child, symptoms[1:])


def optimal_tree(records, symptoms, depth):
    best_tree = None
    sccor = -1

    for symptomim in itertools.combinations(symptoms, depth):
        tree = build_tree(records, symptomim)
        rate = Diagnoser(tree).calculate_success_rate(records)
        if rate > sccor:
            sccor = rate
            best_tree = tree
    return best_tree

if __name__ == "__main__":
    # record1 = Record("influenza", ["cough", "fever"])
    # record2 = Record("influenza", ["cough"])
    # record3 = Record("cold", ["1","cough"])
    # #record4 = Record("cold", ["3", "cough"])
    # records = [record1,record3,record2]
    # pretty_print(optimal_tree(records, ["cough", "fever"], 2))

    ###################  general test ###########################################

    records = parse_data("Data/tiny_data.txt")
    sim = list({ill for j in records for ill in j.symptoms})
    tree=optimal_tree(records,sim,len(sim))
    pretty_print(tree)
    diagnoz = Diagnoser(tree)
    print(diagnoz.all_illnesses())
    print(diagnoz.paths_to_illness("strep"))
