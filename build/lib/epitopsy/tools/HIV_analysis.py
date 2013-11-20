# -*- coding: utf-8 -*-
"""
Created on Mon May 27 12:28:33 2013

@author: Christoph Wilms
"""

import os
import gzip
import shutil

import numpy as np

from Bio import AlignIO

from epitopsy.tools.style import style
from epitopsy.tools.dca_stuff import read_dca

hiv_dict = {}
hiv_dict["x_dca_new"] = 0.8
hiv_dict["gp120_struct_start"] = 33
hiv_dict["gp120_struct_end"] = 492
hiv_dict["gp41_struct_start"] = 512
hiv_dict["gp41_struct_end"] = 856
hiv_dict["wd_source"] = os.path.join(os.getenv("HOME"),
                              "projects",
                              "Project6-Direct_Information",
                              "interesting_proteins")

def get_HXB2_seq_map(ref_seq):
    '''
    Args:
        ref_seq -> HXB2 sequence with gaps

    Returns:
        A dictionary mapping the di coloumns (starting with 1, ..., N) to
        the HXB2 sequence (starting with 1, ..., N).
    '''
    counter = 1
    seq_map = {}
    for i, char in enumerate(ref_seq):
        if char != "-" and char != ".":
            seq_map[i+1] = counter
            counter += 1

    return seq_map


def map_dca(dca_list, seq_map, seq, seq_dist):
    '''
    Args:
        dca_list -> [i,j,di], sorted with respect to di
        seq_map -> maps the di coloums (1,...,N) to the sequence positions
                    (1, ..., N)
        seq -> sequence without gaps
        seq_dist -> minimum  distance seperating two residues to be relevant

    Returns:
        A list of list, containing the di coloumns (1,...,N) which refer to
        HXB2, the sequence positions (1,...,N) itself, the domain
        information and the di value.
    '''
    filtered_dca_list = []
    filtered_dca_list.append(["i", "i_seq", "i_domain", "i_resn",
                              "j", "j_seq", "j_domain", "j_resn",
                              "di", "seq_neighbor"])
    for i, line in enumerate(dca_list):
        i_pos = line[0]
        j_pos = line[1]
        di = line[2]
        if i_pos in seq_map and j_pos in seq_map:
            i_seq = seq_map[i_pos]
            i_domain = map_domain(i_seq)
            i_resn = seq[i_seq-1]
            j_seq = seq_map[j_pos]
            j_domain = map_domain(j_seq)
            j_resn = seq[j_seq-1]
            if np.abs(i_seq - j_seq) > seq_dist:
                seq_neighbor = 0
            else:
                seq_neighbor = 1
        
            new_line = [i_pos, i_seq, i_domain, i_resn,
                        j_pos, j_seq, j_domain, j_resn,
                        di, seq_neighbor]
            filtered_dca_list.append(new_line)

    return filtered_dca_list


def reduce_dca(dca_list, n_dis, include_neighbors):
    '''
    Args:
        dca_list -> [*,di,seq_neighbor], sorted with respect to di
        n_dis -> number of di values
        include_neighbors -> include neighbors or not (True/False)

    Returns:
        A list.
    '''
    header = dca_list[0]
    neighbor_pos = header.index('seq_neighbor')
    new_list = [header]
    counter = 0
    include = True
    for i_dca, dca_item in enumerate(dca_list[1:]):
        if not include_neighbors:
            seq_neighbor = dca_item[neighbor_pos]
            if seq_neighbor == 1:
                # it is a neighbor, do not include it
                include = False
            else:
                # not a neighbor include it
                include = True

        if include:
            new_list.append(dca_item)
            counter += 1

        if counter >= n_dis:
            break

    return new_list

def filter_dca(dca_list, seq_map, seq, seq_dist, n_dis):
    '''
    Args:
        dca_list -> [i,j,di], sorted with respect to di
        seq_map -> maps the di coloums (1,...,N) to the sequence positions
                    (1, ..., N)
        seq -> sequence without gaps
        seq_dist -> minimum  distance seperating two residues to be relevant
        n_dis -> number of di values to check

    Returns:
        A list of list, containing the di coloumns (1,...,N) which refer to
        HXB2, the sequence positions (1,...,N) itself, the domain
        information and the di value.
    '''
    filtered_dca_list = []
    filtered_dca_list.append(["i", "i_seq", "i_domain", "i_resn",
                              "j", "j_seq", "j_domain", "j_resn",
                              "di"])
    counter = 0.
    for i, line in enumerate(dca_list):
        i_pos = line[0]
        j_pos = line[1]
        di = line[2]
        if i_pos in seq_map and j_pos in seq_map:
            i_seq = seq_map[i_pos]
            i_domain = map_domain(i_seq)
            i_resn = seq[i_seq-1]
            j_seq = seq_map[j_pos]
            j_domain = map_domain(j_seq)
            j_resn = seq[j_seq-1]
            counter += 1
            if np.abs(i_seq - j_seq) > seq_dist:
                new_line = [i_pos, i_seq, i_domain, i_resn,
                            j_pos, j_seq, j_domain, j_resn,
                            di]
                filtered_dca_list.append(new_line)
        
        if counter >= n_dis:
            break

    return filtered_dca_list


def map_domain(x):
    '''
    Args:
        x -> position in the HXB2 sequence (1,...,N)
    '''
    if 1 <= x <= 32:
        domain = "SP"
    elif 132 <= x <= 156 :
        domain = "V1"
    elif 157 <= x <= 196:
        domain = "V2"
    elif 296 <= x <= 330:
        domain = "V3"
    elif 385 <= x <= 418:
        domain = "V4"
    elif 461 <= x <= 471:
        domain = "V5"
    elif 33 <= x <= 118 or 206 <= x <= 252 or 476 <= x <= 511:
        domain = "ID"
    elif 253 <= x <= 397 or 410 <= x <= 419 or 437 <= x <= 475:
        domain = "OD"
    elif 118 <= x <= 205 or 421 <= x <= 437:
        domain = "BS"
    elif 512 <= x <= 856:
        domain = "GP41"
    else:
        domain = "GP120"

    return domain

def get_color_for_domain(domain):
    if domain == "ID":
        color = "magenta"
    elif domain == "OD":
        color = "yellow"
    elif domain == "BS":
        color = "blue"
    elif domain == "V1":
        color = "red"
    elif domain == "V2":
        color = "green"
    elif domain == "V3":
        color = "orange"
    elif domain == "V4":
        color = "dirtyviolet"
    elif domain == "V5":
        color = "teal"
    else:
        color = "gray60"

    return color

