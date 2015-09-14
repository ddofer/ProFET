#!/usr/bin/python3
"""ProFET is a program for feature extraction from fasta files.
4 main options:
1. Extract features from fasta files(s) (Requires data)
2. Train a classifier (Requires training data)
3. Extract classifier's performance (Requires trained model and testing set)
4. Predict new data sets (Needs a trained model and input data set)

Training data sets can be in one of the two formats:
a. Each class in a seperated file. For example:

b. All data in the same file, but class can be extract from sequence identifier.
For example:

"""

import argparse
import pickle
import json
import sys
import time
import pandas as pd
import numpy as np
import pprint as pp
from Bio.SeqIO.FastaIO import SimpleFastaParser
from multiprocessing import Pool, Manager
from os.path import basename, exists
from functools import lru_cache as memoized

#Internal imports
from FeatureGen import Get_Protein_Feat #To generate features
from Model_trainer import trainClassifier, load_data
from GetPredictorPerf import get_scores

# from IPython.core.debugger import Tracer #TO REMOVE!!!
# import warnings
# import traceback

# def warn_with_traceback(message, category, filename, lineno, file=None, line=None):
#     traceback.print_stack()
#     log = file if hasattr(file,'write') else sys.stderr
#     log.write(warnings.formatwarning(message, category, filename, lineno, line))
# warnings.showwarnings = warn_with_traceback
# warnings.simplefilter("always")
# np.seterr(all='raise')

def get_params():
    '''Parse arguments'''
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset', '-r', dest='dataset', action='append',
        type=str, help='The path to the data set (fasta file/s. May serve as training, testing, new set, or feature extracting')
    #parser.add_argument('--resultsDir', '-rs', dest='resultsDir', type=str,
    #help='The path to directory to write the results files')
    parser.add_argument('--extractfeat', '-f', dest='extractfeatures',
        action='store_true', default=False, help='Extract features')
    parser.add_argument('--performance', '-p', dest='performance',
        action='store_true', default=False, help='Print performance of the model on a given data (test) set')
    parser.add_argument('--classformat', '-c', dest='classformat',
        type=str, default='file', help='Defines the class of each '\
         'sequence, by on of the next options: \'dir\', \'file\', or \'id\'."')
    parser.add_argument('--outputmodel', '-o', dest='outputmodel',
        default=None, type=str, help='Filename for saving the model')
    parser.add_argument('--inputmodel', '-i', dest='inputmodel',
        default=None, type=str, help='Filename for saving the model')
    parser.add_argument('--outputfile', '-s', dest='outputfile',
        default=None, type=str, help='Filename for saving the '\
        'output (features, or predicted labels)')
    parser.add_argument('--classifier', '-t', dest='classifiertype',
        default='forest', help='The type of the classifier model. '\
        'Can be one of next options: forest')

    results = parser.parse_args()
    return results

@memoized(maxsize=32)
def get_label_from_filename(filename):
    '''Extract the label according to filename
    Remove paht, and extensions
    For example:
    '../Mamm_Organellas/Mammal_peroxisome9.features.csv' -> 'Mammal_peroxisome9'
    '''
    return basename(filename).split('.')[0]

def write_csv(features_dict, out_file=None):
    '''If out_file is not given, return features as a pandas.dataFrame object
    '''
    #Since different sequences can have different set of features(!) get the union of features:
    feature_names = set().union(*[tuple(val.keys()) for val in features_dict.values()])
    if out_file:
        f = open(out_file,'w')
    else:
        csv_str = ''
    header = sorted(feature_names)
    #Print feature names
    if out_file:
        f.write('accession\t' + '\t'.join(header) + '\n')
    else:
        csv_str += 'accession\t' + '\t'.join(header) + '\n'
    values_line_fmt = '%s\t' * len(header) + '%s\n'
    for acc, features in features_dict.items():
        #If feature doesn't exists, put a 0
        values_line = acc + '\t' + '\t'.join([str(features.get(f_name, 0)) for f_name in header]) + '\n'
        if out_file:
            f.write(values_line)
        else:
            csv_str += values_line
    if not out_file: return load_data(csv_str)

def write_features_to_file(features_dict, out_file):
    '''Write in json format'''
    with open(out_file, 'w') as f:
        f.write(json.dumps(features_dict, sort_keys=True, indent=4))

def load_features_file(features_file):
    '''Load json format file'''
    with open(features_file) as f:
        return json.load(f)

def update_dict_with_features(d, seq_id, seq, classname=None):
    '''Extract features and update the dict
    A worker function to be executed in paralel with multiprocessing.
    '''
    features = Get_Protein_Feat(seq)
    if classname:
        features.update({'classname': classname})
    d[seq_id] = features

@memoized(maxsize=32)
def extract_features(sequences_file, classformat=None, force=False):
    '''Given a fasta file with sequences, and an output file,
    Extract features, write it to file, and return the features as a dict
    of dicts
    Uses mutliprocessin to accelerate.
    '''
    features_dict_file = sequences_file + '_%s.pkl' % classformat
    if exists(features_dict_file) and not force:
        fh = open(features_dict_file, 'rb')
        features = pickle.load(fh)
        fh.close()
        return features
    pool = Pool()
    manager = Manager()
    features_dict = manager.dict()
    if classformat == 'file':
        classname = get_label_from_filename(sequences_file)
    else:
        classname = None
    with open(sequences_file) as f:
        for record in SimpleFastaParser(f):
            p = pool.apply_async(update_dict_with_features, args=(features_dict, record[0], record[1], classname))
            print('.', end='')
            sys.stdout.flush()
        pool.close()
        pool.join()
    with open(features_dict_file, 'wb') as out_fh:
        features_dict = dict(features_dict)
        pickle.dump(features_dict, out_fh, protocol=pickle.HIGHEST_PROTOCOL)
    return features_dict

def extract_datasets_features(datasets, classformat, output_file=None):
    '''For each file in dataset execut get features, and save to file, or save to one file if output file is given
    '''
    if output_file and exists(output_file):
        #Returns filename, maybe change to dataframe object
        return output_file
    all_features_dict = {}
    for filename in datasets:
        print('Generating features for file:', filename)
        features_dict = extract_features(filename, classformat)
        all_features_dict.update(features_dict)
        if not output_file:
            write_csv(features_dict, filename + '.features.csv')
        print('Done')
    if output_file:
        print('Writing all features to:', output_file)
        write_csv(all_features_dict, output_file)
    return all_features_dict

def train_model(dataset_files, classformat, outputmodel, classifiertype):
    '''Given set of files with fasta sequences, class format (e.g., file), 
    filename to save model and required model type (e.g., forest)
    Train the model and save it to the given file
    '''
    output_features_file='ProFET_features.csv'
    features_dict = extract_datasets_features(dataset_files, classformat, output_features_file)
    #features_df = load_data(output_features_file)
    print('Learning %s model' % classifiertype)
    model, label_encoder, scaler, feature_names = trainClassifier(output_features_file, classifiertype, kbest=0, alpha=False, optimalFlag=False, normFlag=True)
    #Save model and additional data to file
    pickle.dump((model, label_encoder, scaler, feature_names), open(outputmodel,'wb'), protocol=pickle.HIGHEST_PROTOCOL)
    print('Done')

def load_model_from_file(filename):
    '''Given a pickle filename return model, label encoder and a scaler'''
    #return pickle.load(open(filename, 'rb'))
    fh = open(filename, 'rb')
    model_obj = pickle.load(fh)
    fh.close()
    return model_obj

def remove_nan_vals(array):
    '''Remove rows with nan values and return the new array'''
    row_len = array.shape[1]
    while np.isnan(array.min()):
        ind = array.argmin()
        array = np.delete(array,(int(ind/row_len)), axis=0)
    return array

def class_data(inputmodel_file, dataset, outfile=None):
    '''Given a file of trained model, a dataset (can be multiple files) print for
    each sequence the accession and the predicted label.
    If outfile is given, print to this file instead of STDOUT
    '''
    #Load model
    model, label_encoder, scaler, model_feature_names = load_model_from_file(inputmodel_file)
    if outfile:
        f = open(outfile, 'w')
    else:
        f = sys.stdout
    for sequences_file in dataset:
        res = write_csv(extract_features(sequences_file, classformat=None))
        features = res[0]
        accessions = res[1]
        data_feature_names = res[-1]
        features = match_features(model_feature_names, features, data_feature_names)
        #Remove nan values
        features = remove_nan_vals(features)
        #Predict
        scaled_features = scaler.transform(features)
        labels_idx = model.predict(scaled_features)
        labels = label_encoder.inverse_transform(labels_idx)
        #Tracer()() #TO REMOVE!!!
        for acc, label in zip(accessions.values, labels):
            f.write('%s\t%s\n' % (acc, label))
    if outfile:
        f.close()

def match_features(main_features_list, features, minor_features_list):
    '''Given feature names list, and a features dataFrame,
    Add column with zeros for missing features from the dataFrame'''
    main_features_set = set(main_features_list)
    minor_features_set = set(minor_features_list)

    #Remove spare features in data
    for feature_name in minor_features_set - main_features_set:
        idx = np.where(minor_features_list==feature_name)[0][0]
        #Remove a column
        features = np.delete(features, idx, axis=1)
        #Remove name from list to keep indexing consictency
        #minor_features_list.remove(feature_name)
        minor_features_list = np.delete(minor_features_list, idx)
        #raise 'Data contains more features than the model!'

    additional_features = main_features_set - minor_features_set
    #Find missing features in data
    mask = np.in1d(main_features_list, list(additional_features), assume_unique=True)
    #Find the indexes
    idxs = sorted(np.where(mask))
    for idx in idxs:
        #idx = np.where(main_features_list==feature_name)[0][0]
        #Add zeros for missing feature (column idx)
        features = np.insert(features, idx, 0, axis=1)
    return features

def get_classifier_performance(inputmodel, dataset, classformat):
    '''
    '''
    #Load Model
    model, model_label_encoder, model_scaler, model_feature_names = load_model_from_file(inputmodel)
    true_labels = np.array([])
    predicted_labels = np.array([])

    for sequences_file in dataset:
        #Extract data set features and true lables
        features, part_true_labels, label_encoder, feature_names = write_csv(extract_features(sequences_file, classformat=classformat))
        true_labels = np.append(true_labels, model_label_encoder.transform(label_encoder.inverse_transform(part_true_labels)))
        features = match_features(model_feature_names, features, feature_names)
        
        #Predict labels
        scaled_features = model_scaler.transform(features)
        labels_encoded = model.predict(scaled_features)
        #labels = label_encoder.inverse_transform(labels_encoded)
        predicted_labels = np.append(predicted_labels, labels_encoded)
    #Generate scores
    results = get_scores(predicted_labels, true_labels, verbose=False)
    return results

def main():
    results = get_params()
    dataset = results.dataset
    extractfeatures = results.extractfeatures
    classformat = results.classformat
    outputmodel = results.outputmodel
    inputmodel = results.inputmodel
    outputfile = results.outputfile
    classifiertype = results.classifiertype
    performance = results.performance
    
    #Four options for program: (i) Extract features, (ii) Learn model, or given a model:
    # (iii) Classify new data, (iv) Extreact performance

    #(i) Extract features
    if extractfeatures:
        if not dataset or not outputfile:
            print('Extracting features requires a training dataset and '\
            'an ouptput file to be given')
            exit()
        else:
            print('Only extracting Features')
            features = extract_datasets_features(dataset, classformat=None, output_file=outputfile)
    #(ii) Learn model
    elif outputmodel:
        if not dataset:
            print('Training a model requires a training dataset will be'\
            'given')
            exit()
        else:
            train_model(dataset, classformat, outputmodel, classifiertype)
    #or given a model:
    elif inputmodel:
        # (iii) Classify new data
        if dataset and not performance:
            class_data(inputmodel, dataset, outputfile)
        #(iv) Extreact performance
        else:
            results = get_classifier_performance(inputmodel, dataset, classformat)
            pp.pprint(results)
    else:
        print('Got wrong combination of parameters, please refer '\
        'to help or README')
        exit()

if __name__=="__main__":
    main()
    # classformat = 'file'
    # inputmodel = 'chap_model'
    # model, label_encoder, scaler, model_feature_names = load_model_from_file(inputmodel)
    # sequences_file = '../Chap/test/negative.fasta'
    # res = write_csv(extract_features(sequences_file, classformat=None))
    # features = res[0]
    # accessions = res[1]
    # data_feature_names = res[-1]
    # features = match_features(model_feature_names, features, data_feature_names)
    
    # #Predict
    # scaled_features = scaler.transform(features)
    # labels_idx = model.predict(scaled_features)
    # labels = label_encoder.inverse_transform(labels_idx)

    # features1, labels1, label_encoder1, feature_names1 = write_csv(extract_features(sequences_file1, classformat=classformat))
    # # sequences_file2 = '../Mamm_Organellas/Mammal_melanosome_0.9.fasta'
    # # #features2, labels2, label_encoder2, feature_names2 = write_csv(extract_features(sequences_file, classformat=classformat))
    # sequences_file2 = '../Chap/test/positive.fasta'
    # dataset = [sequences_file1, sequences_file2]
    # results = get_classifier_performance(inputmodel, dataset, classformat)
#Example of execution:
#python3 ProFET.py --dataset ../Chap/train/positive.fasta --dataset ../Chap/train/negative.fasta --classformat file --outputmodel test_model
#python3 ProFET.py --dataset ../Chap/test/positive.fasta --dataset ../Chap/test/negative.fasta --inputmodel chap_model --performance
    
