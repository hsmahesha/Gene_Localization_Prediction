import sys
import os
import csv
import copy
from collections import OrderedDict


#------------------------------------------------------------------------------#
# compute the accuray of the implemented nearest neighbor prediction algorithm 
#------------------------------------------------------------------------------#
def compute_accuracy(pred_dict, actual_pred_dict):
    succ_count = 0
    tot_count = 0 
    
    for key, value in pred_dict.items():
        act_loc = actual_pred_dict[key]
        if act_loc in value:
           succ_count += 1
        tot_count += 1

    accuracy = (float(succ_count) / float(tot_count)) * 100.0
    return accuracy 
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# construct dictionary data structure from actual prediction data list
#------------------------------------------------------------------------------#
def construct_actual_predict_dictionary(act_pred_data):
    act_pred_dict = {}
    for row in act_pred_data:
        act_pred_dict[row[0]] = row[1]

    return act_pred_dict
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# helper routines for predict_gene_localization()
#------------------------------------------------------------------------------#
def get_localization(NN):                                                       
    loc_dict = {}                                                               
    for n_key, n_value in NN.items():                                           
        loc = n_value["LOCALIZATION"]                                           
        if loc not in loc_dict:                                                 
           loc_dict[loc] = 0                                                    
        loc_dict[loc] += 1                                                      
    if len(loc_dict) == 0:                                                      
       return "none"                                                            
    loc_list = sorted(loc_dict.items(), key=lambda x: x[1], reverse=True)       
    max_val = loc_list[0][1]                                                    
    res_list = []                                                               
    i = 0                                                                       
    while i < len(loc_list) and loc_list[i][1] == max_val:                      
       res_list.append(loc_list[i][0])                                         
       i += 1                                                                  
    return res_list

def is_agree(n_key, t_list, n_list, feature):
    if feature != "INTERACTION":
       size = len(t_list)
       for i in range(0, size):
           if t_list[i] == n_list[i]:
              return True
    else:
       if n_key in t_list:
          return True

    return False

def get_nn_data_dict(NN, t_value, feature):
    new_NN = {}
    t_list = t_value[feature]
    for n_key, n_value in NN.items():
        n_list = n_value[feature]
        if is_agree(n_key, t_list, n_list, feature):
           new_NN[n_key] = copy.deepcopy(NN[n_key])
    return new_NN

def get_nearest_neighbours(t_key, t_value, PFL, NN):
    for feature in PFL:
        new_NN = get_nn_data_dict(NN, t_value, feature)
        if len(new_NN) != 0:
           NN = copy.deepcopy(new_NN)
    return NN
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# predict localization of genes from test data
#------------------------------------------------------------------------------#
def predict_gene_localization(tr_dict, ts_dict):
    os.system("clear")
    print("\n\nBe patient... I am learning and predicting...")

    PFL = ["COMPLEX", "CLASS", "INTERACTION", "MOTIF"]

    pred_dict = { }
    for key, value in ts_dict.items():
        NN = get_nearest_neighbours(key, value, PFL, tr_dict)
        pred_list = get_localization(NN)
        pred_dict[key] = pred_list

    return pred_dict
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# create appropriate dictionary data structures from data
#------------------------------------------------------------------------------#
def construct_data_dictionary(data, inter):
    data_dict = {}

    for row in data:
        gene = row[0]
        c_list = list(row[2:26])
        co_list = list(row[26:82])
        m_list = list(row[93:444])
        l = row[2959]

        gene_dict = {}

        gene_dict["GENE"] = gene
        gene_dict["CLASS"] = c_list
        gene_dict["COMPLEX"] = co_list
        gene_dict["MOTIF"] = m_list
        gene_dict["LOCALIZATION"] = l
        gene_dict["INTERACTION"] = list()

        data_dict[gene] = gene_dict

    for row in inter:
        gene1 = row[0]
        gene2 = row[1]
        if gene1 in data_dict:
           data_dict[gene1]["INTERACTION"].append(gene2)
        if gene2 in data_dict:
           data_dict[gene2]["INTERACTION"].append(gene1)

    return data_dict
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# close input data files
#------------------------------------------------------------------------------#
def close_files(train_data_file, train_inter_file, test_data_file, \
               test_inter_file):
    train_data_file.close()
    train_inter_file.close()
    test_data_file.close()
    test_inter_file.close()
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# read input data files
#------------------------------------------------------------------------------#
def read_file(ifile):
    try:
        data = list(csv.reader(ifile, delimiter=','))
    except IOError:
        print("\nError: Failed to read input file")
        sys.exit()
    return data

def read_files(train_data_file, train_inter_file, test_data_file, \
               test_inter_file):
    tr_data = read_file(train_data_file)
    tr_inter = read_file(train_inter_file)
    ts_data = read_file(test_data_file)
    ts_inter = read_file(test_inter_file)

    return tr_data, tr_inter, ts_data, ts_inter
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# open input data files
#------------------------------------------------------------------------------#
def open_file(fname, mode):
    try:
        ifile = open(fname, mode)
    except IOError:
        print("\nError: Failed to open input file " + fname)
        sys.exit()
    return ifile

def open_files():
    train_data_file = open_file("./data/train/Full_File.data", "r")
    train_inter_file = open_file("./data/train/Interactions_relation.data", "r")
    test_data_file = open_file("./data/test/Full_File.test", "r")
    test_inter_file = open_file("./data/test/Interactions_relation.test", "r")

    return train_data_file, train_inter_file, test_data_file, test_inter_file
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# write result to output result.txt file
#------------------------------------------------------------------------------#
def write_result(pred_dict, act_pred_dict, accuracy):
    rfile = open_file("./result.txt", "w")

    rfile.write("Gene No      Actual Loc\t\t\tPredicted Locs\n\n")

    for key, value in act_pred_dict.items():
        pred_locs = pred_dict[key]
        rfile.write(key + "      " + value + "\t\t\t" + str(pred_locs) + "\n")

    rfile.write("\n\n");
    rfile.write("Accuracy = " + str(accuracy) + "%") 
    rfile.write("\n\n")
    rfile.close()
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# main entry point api function
#------------------------------------------------------------------------------#
def main():
    # open input files
    train_data_file, train_inter_file, \
                     test_data_file, test_inter_file = open_files()
   
    # read input files
    tr_data, tr_inter, ts_data, ts_inter = \
                  read_files(train_data_file, \
                             train_inter_file, test_data_file, test_inter_file)             

    # close input files
    close_files(train_data_file, train_inter_file, test_data_file, \
                test_inter_file)

    # construct dict structures from input data
    tr_dict = construct_data_dictionary(tr_data, tr_inter)
    ts_dict = construct_data_dictionary(ts_data, ts_inter)

    # predict localization using nearest neightbor algorithm
    pred_dict = predict_gene_localization(tr_dict, ts_dict)

    # open actual prediction file and read it
    act_pred_file = open_file("./correct_result/LocalizationKey.txt", "r")
    act_pred_data = read_file(act_pred_file)
    act_pred_dict = construct_actual_predict_dictionary(act_pred_data)

    # compute accuracy of the model
    accuracy = compute_accuracy(pred_dict, act_pred_dict)

    # write result to output file
    write_result(pred_dict, act_pred_dict, accuracy)

    print("\n\n")
    print("...done with prediction, result is available in ./result.txt file");
    print("\n\n")
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# call to main entry point api function
#------------------------------------------------------------------------------#
main()
#------------------------------------------------------------------------------#
