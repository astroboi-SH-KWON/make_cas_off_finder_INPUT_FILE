from os import listdir
from os.path import isfile, join
import pandas as pd
import openpyxl
from time import clock
import random
import math

class Utils:
    def __init__(self):
        self.ext_txt = ".txt"
        self.ext_dat = ".dat"
        self.ext_xlsx = ".xlsx"

    def read_excel_2_dataframe(self, path, f_name):
        return pd.read_excel(path + f_name + self.ext_xlsx)

    """
    make_txt_cas_off_finder_input : .txt input file for CAS_OFF_FINDER
    :param
        seq_set_pool
        init
            pam_seq
            spacr_len
            mis_mtch_num : maximum # of mismatch count
            sub_set_num : sub set # of seq_set_pool to make .txt input files separately
            path
            ref_srv_path
    """
    def make_txt_cas_off_finder_input(self, seq_set_pool, init):
        pam_seq = init[0]
        spacr_len = init[1]
        mis_mtch_num = init[2]
        sub_set_num = init[3]
        path = init[4]
        ref_srv_path = init[5]

        raw_each_sub_len = len(seq_set_pool) / sub_set_num
        each_sub_len = int(raw_each_sub_len)
        # * 1000) / 1000 : raw_each_sub_len 값이 반욜림되는 경우의 error 방지
        raw_last_sub_len_to_add = int(((raw_each_sub_len - each_sub_len) * sub_set_num) * 1000)/1000
        last_sub_len_to_add = math.ceil(raw_last_sub_len_to_add)
        print("total len : " + str(len(seq_set_pool)) + ", raw_each_sub_len : " + str(raw_each_sub_len))
        print("raw_last_sub_len_to_add : " + str(raw_last_sub_len_to_add) + " , last_sub_len_to_add : " + str(last_sub_len_to_add))

        for fname_idx in range(sub_set_num):
            if fname_idx == sub_set_num -1:
                each_sub_len += last_sub_len_to_add
            tmp_sub_set = set(random.sample(seq_set_pool, each_sub_len))
            seq_set_pool -= tmp_sub_set
            print("tmp_sub_set_" + str(fname_idx) + " len : " + str(len(tmp_sub_set)))
            print("rest total len : " + str(len(seq_set_pool)))
            print(" ")
            with open(path + str(fname_idx) + self.ext_txt, 'a') as f:
                f.write(ref_srv_path + "\n")
                f.write("N"*spacr_len + pam_seq + "\n")
                for data_str in tmp_sub_set:
                    f.write(data_str + " " + str(mis_mtch_num) +"\n")

