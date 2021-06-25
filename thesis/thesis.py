import argparse
from collections import defaultdict
from requests import get
from json import dump
from Bio import pairwise2


class Thesis:
    miRna_dict = {}
    miRna_list = []
    pre_miRna_list = []
    pre_miRNa_dict = {}
    motifs = []
    miRna_motifs_dict = {}
    premiRna_motifs = {}
    min_percentage_identity = 80

    '''
    Class constructor 
    '''

    def __init__(self, miRna_file, motif_database, percentage_identity=80, organism=None, pre_miRna_file=None):
        """

        :param miRna_file: name of miRna file
        :param pre_miRna_file: name of pre-miRna file
        :param motif_database: name of motif database
        :param organism: the given organism from which the sequence is derived
        :param percentage_identity: minimum percent identity used for miRna and pre-miRna clustering
        """
        if miRna_file is None or motif_database is None:
            none = [i for i in [miRna_file, motif_database] if i is None]
            print(f'Needed parameter: {none}')
            return
        self.min_percentage_identity = percentage_identity
        self.run(miRna_file, pre_miRna_file, motif_database, organism)

    '''
    Main method in which runs the program
    '''

    def run(self, miRna_file, pre_miRna_file, motif_database, organism):
        """

        :param miRna_file: name of miRna file
        :param pre_miRna_file: name of pre-miRna file
        :param motif_database: name of motif database
        :param organism: the given organism from which the sequence is derived
        """
        print("Analysis for" + miRna_file)
        print("#Step 1 - Read from file")
        nr_file = 0
        self.motifs = self.read_from_file(motif_database, self.motifs, organism, nr_file)
        nr_file += 1
        self.miRna_list = self.read_from_file(miRna_file, self.miRna_list, organism, nr_file)
        print("#Step 2 - Group sequences")
        self.miRna_dict = self.check_identity(self.miRna_list)
        print("#Step 3 - Find motifs and remove duplicates")
        self.miRna_motifs()
        self.save_to_file(self.miRna_dict, miRna_file + "_before_reduction")
        self.miRna_dict = self.remove_motifs_duplicates(self.miRna_dict)
        self.save_to_file(self.miRna_dict, miRna_file + "_after_reduction")
        copy = self.miRna_dict
        print("#Step 4 - Download more information")
        self.miRna_dict = self.download_gene(self.miRna_dict)
        print("#Step 5 - Find homologies and group")
        self.miRna_dict = self.find_homologies(self.miRna_dict, organism)
        self.save_to_file(self.miRna_dict, miRna_file + "_final_miRna")
        print("#Step 6 - Count motifs")
        self.miRna_motifs_dict = self.count_motifs(self.miRna_dict)
        self.save_to_file(self.miRna_motifs_dict, miRna_file + "_count_motifs_miRNA")
        if pre_miRna_file is not None:
            print("Analysis for" + pre_miRna_file)
            print("#Step 1 - Read from file")
            self.pre_miRna_list = self.read_from_file(pre_miRna_file, self.pre_miRna_list, organism, nr_file)
            print("#Step 2 - Group sequences")
            self.pre_miRNa_dict = self.check_identity(self.pre_miRna_list)
            print(self.pre_miRNa_dict)
            print("#Step 3 - Find motifs and remove duplicates")
            self.pre_miRna_motifs(copy)
            self.save_to_file(self.pre_miRNa_dict, pre_miRna_file + "before_reduction")
            self.pre_miRNa_dict = self.remove_motifs_duplicates(self.pre_miRNa_dict)
            self.save_to_file(self.pre_miRNa_dict, pre_miRna_file + "after_reduction")
            print("#Step 5 - Download more information")
            self.pre_miRNa_dict = self.download_gene(self.pre_miRNa_dict)
            print("#Step 5 - Find homologies and group")
            self.pre_miRNa_dict = self.find_homologies(self.pre_miRNa_dict, organism)
            self.save_to_file(self.pre_miRNa_dict, pre_miRna_file + "_final_pre_miRna")
            print("#Step 6 - Count motifs")
            self.premiRna_motifs = self.count_motifs(self.pre_miRNa_dict)
            self.save_to_file(self.premiRna_motifs, pre_miRna_file + "_count_motifs_pre_miRNA")

    '''
    Read from file miRna, pre-miRna and motif database
    '''

    def read_from_file(self, file_name, list_to_save, organism, nr_file):
        """

        :param nr_file: 0 - miRna, 1 pre-miRna
        :param file_name: name of file in which are sequences of miRna, pre-miRna or motif database
        :param list_to_save: list in which sequences will be saved
        :param organism: given organism
        :return: list with saved sequences
        """
        f = open(file_name, "r")
        if nr_file == 0:
            lines = f.readlines()
            for line in lines:
                motif = line.split("\t")
                if motif[3].lower() == organism.lower():
                    list_to_save.append(motif)
        else:
            lines = f.readlines()
            for line in range(len(lines)):
                if lines[line].startswith(">"):
                    org = organism.lower().split('_')
                    result = all(elem in lines[line].lower() for elem in org)
                    if result:
                        list_to_save.append([lines[line].rstrip(), lines[line + 1].rstrip()])
            self.remove_duplicates(list_to_save)
        f.close()
        return list_to_save

    '''
    Check identity percentage between sequences and group them if they have more than min_percentage_identity value.
    '''

    def check_identity(self, list_seq):
        """

        :param list_seq: list of miRna or pre-miRna
        :return: dictionary with clustered sequences
        """
        dictionary = defaultdict(list)
        index = 0
        already_used = []
        for first in list_seq:
            similar = []
            if self.check_sequence(first[1], already_used):
                for second in list_seq:
                    if self.check_sequence(second[1], already_used):
                        if first != second:
                            alignments = pairwise2.align.globalxx(first[1], second[1], score_only=True)
                            divider = len(first[1]) if len(first[1]) < len(second[1]) else len(second[1])
                            if alignments / divider >= self.min_percentage_identity / 100:
                                if len(similar) == 0:
                                    similar.append(first)
                                    similar.append(second)
                                else:
                                    similar.append(second)
                if len(similar) == 0:
                    similar.append(first)
                dictionary[index] = similar
                index += 1
                already_used.extend(similar)
        return dictionary

    '''
    Search for motifs in miRna
    '''

    def miRna_motifs(self):
        for index in self.miRna_dict:
            result = []
            for seq in self.miRna_dict[index]:
                for motif in self.motifs:
                    if motif[4] in seq[1]:
                        result.append([seq[0], seq[1], motif[1], motif[4]])
            self.miRna_dict[index] = result

    '''
    Search for motifs in pre_miRna excluding motifs found in miRna
    '''

    def pre_miRna_motifs(self, miRna_dict):
        for index_p in self.pre_miRNa_dict:
            result = []
            for pre_miRna in self.pre_miRNa_dict[index_p]:
                for index_m in miRna_dict:
                    for miRna in miRna_dict[index_m]:
                        if miRna[1] in pre_miRna[1]:
                            pos1 = pre_miRna[1].find(miRna[1])
                            for motif in self.motifs:
                                if motif[4] in pre_miRna[1]:
                                    pos2 = pre_miRna[1].find(motif[4])
                                    if pos2 < pos1 or pos2 + len(motif[4]) > pos1 + len(miRna[1]):
                                        result.append([pre_miRna[0], pre_miRna[1], motif[1], motif[4]])
            self.pre_miRNa_dict[index_p] = result

    def count_motifs(self, dictionary):
        temp_dict = {}
        for motif in self.motifs:
            if motif[4] not in temp_dict:
                found = []
                count = 0
                for index in dictionary:
                    for name in dictionary[index]:
                        if name != "homologies":
                            for elem in dictionary[index][name]['motifs']:
                                if motif[4] == elem[3]:
                                    count += 1
                                    found.append([elem[0], elem[1]])
                if count != 0:
                    temp_dict[motif[4]] = {'count': count, 'found': found}
        return temp_dict

    @staticmethod
    def find_homologies(dictionary, organism):
        server = "https://rest.ensembl.org"
        already_used = []
        temp_dict = {}
        index = 0
        for elem in dictionary.keys():
            if elem not in already_used:
                homologies = []
                ext = f"/homology/id/{elem}"
                r = get(server + ext, headers={"Content-Type": "application/json"})
                if r.ok:
                    decoded = r.json()
                    for item in decoded['data']:
                        homologies = [i['target']['id'] for i in item['homologies']
                                      if i['target']['species'] == organism.lower()]
                already_used.append(elem)
                temp_dict[index] = {'homologies': [elem], elem: dictionary[elem]}
                for item in homologies:
                    if item in dictionary.keys():
                        homo = dictionary.get(item)
                        temp_dict[index][item] = homo
                        temp_dict[index]['homologies'].append(item)
                        already_used.append(item)
                index += 1
        return temp_dict

    @staticmethod
    def download_gene(dictionary):
        server = "https://rest.ensembl.org"
        temp_dict = {}
        for index in dictionary:
            for elem in dictionary[index]:
                if elem[2] in temp_dict:
                    temp_dict[elem[2]]["motifs"].append(elem)
                else:
                    functions = []
                    protein = []
                    ext = f"/xrefs/id/{elem[2]}?external_db=GO;all_levels=1"
                    r = get(server + ext, headers={"Content-Type": "application/json"})

                    if r.ok:
                        decoded = r.json()
                        for item in decoded:
                            if len(functions) == 0 or [item['display_id'], item['description']] not in functions:
                                functions.append([item['display_id'], item['description']])

                    ext = f"/lookup/id/{elem[2]}?expand=1"
                    r = get(server + ext, headers={"Content-Type": "application/json"})
                    if r.ok:
                        decoded = r.json()
                        for item in decoded["Transcript"]:
                            if "Translation" in item:
                                protein.append(item['Translation']['id'])
                    temp_dict[elem[2]] = {
                            "protein_id": protein, "motifs": [elem], "functions": functions
                    }
        return temp_dict

    '''
    Remove duplicates from list of miRna or pre-miRna
    '''

    @staticmethod
    def remove_duplicates(list):
        """

        :param list: list of miRna or pre-miRna
        :return: list without duplicates
        """
        for first in list:
            for second in list:
                if first[0] != second[0]:
                    if first[1] == second[1]:
                        list.remove(second)
        return list

    '''
    Save to file result in json format
    '''

    @staticmethod
    def save_to_file(dictionary, name):
        """

        :param dictionary: dictionary of miRNa or pre-miRNa sequences with motifs
        :param name: name of file in which dictionary will be saved
        """
        name = name.split('/')[2]
        with open(f"./result/{name}.json", "w", encoding='utf-8') as output_file:
            dump(dictionary, output_file, ensure_ascii=False, indent=4)

    '''
    Remove motif duplicates in pre-miRna dictionary and miRna dictionary
    '''

    @staticmethod
    def remove_motifs_duplicates(dictionary):
        """

        :param dictionary: dictionary with miRna or pre-miRna which have motifs
        :return: dictionary without duplicated motifs
        """
        for index in dictionary:
            for first in dictionary[index]:
                for second in dictionary[index]:
                    if first[0] != second[0] and first[3] == second[3] and first[2] == second[2]:
                        dictionary[index].remove(second)
                    else:
                        continue
        return dictionary

    '''
    Check if list already contains sequence
    '''

    @staticmethod
    def check_sequence(seq, already_used):
        """

        :param seq: sequence to check
        :param already_used: list with used sequences
        :return: True / False depending on whether the sequence is in list
        """
        for item in already_used:
            if seq == item[1]:
                return False
        return True


if __name__ == '__main__':
    # python thesis.py -m mature.fa -pre hairpin.fa -db ATtRACT_db.txt -o Homo_sapiens
    ap = argparse.ArgumentParser()
    ap.add_argument("-m", "--miRna", required=True,
                    help="fasta file with miRna sequences")
    ap.add_argument("-pre", "--pre_miRna", required=False,
                    help="fasta file with pre-miRNa sequences")
    ap.add_argument("-db", "--database", required=True,
                    help="database of motifs")
    ap.add_argument("-o", "--organism", required=False,
                    help="organism")
    ap.add_argument("-pi", "--percentage_identity", required=False, default=80,
                    help="The minimum identity percentage that is used to group sequences (default value = 80%%)")
    ap.add_argument("-f", "--functions", required=False,
                    help="Filter motifs with searched function")
    args = ap.parse_args()
    Thesis(miRna_file=args.miRna,
           pre_miRna_file=args.pre_miRna,
           motif_database=args.database,
           organism=args.organism,
           percentage_identity=args.percentage_identity
           )
