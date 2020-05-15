from time import clock

import Util
import Logic

############### start to set env ################
WORK_DIR = "D:/000_WORK/SongMyunJae_YuGooSang/20200417/WORK_DIR/"

PAM_SEQ = ['NGG']
ADD_SEQ1_LEN = 4
SPACER_LEN = 20
ADD_SEQ2_LEN = 3
CLVG_AFTER_PAM = 3

FILE_NAME_LIST = ["X", "Y"]

############### CAS OFF FINDER ################
CAS_OFF_EXCEL = "Chlorocebus_sabaeus_excel_20200424/Chlorocebus_sabaeus_filter_out_05_65_20200506/chlorochebus_sabaeus_w_cas9_scores_"
CAS_OFF_TXT = "Chlorocebus_sabaeus_excel_20200424/Chlorocebus_sabaeus_filter_out_05_65_20200506/chlorochebus_sabaeus_w_cas9_scores_re_off_target_seq_700.6646942.txt"
# CAS_OFF_TXT = "CAS_OFF_FINDER/INPUT/chlorochebus_sabaeus_cas_NGG_off_"
CAS_OFF_TXT_TOP_20 = "CAS_OFF_FINDER/INPUT/chlorochebus_sabaeus_top20_cas_NGG_off_"
MAX_MISMATCH = 3
SUB_SET_NUM = 100
REF_SRV_PATH = "Input/FASTA/chlorocebus_sabaeus_chr"
PAM_N = "NNN"
TOP_N = 20

# INITIAL_CAS_OFF = ['NGG', SPACER_LEN, MAX_MISMATCH, SUB_SET_NUM, WORK_DIR + CAS_OFF_TXT, REF_SRV_PATH]
INITIAL_CAS_OFF = ['NGG', SPACER_LEN, MAX_MISMATCH, 18, WORK_DIR + CAS_OFF_TXT_TOP_20, REF_SRV_PATH]

############### end setting env ################

def make_cas_off_finder_input_top_20_CAS_score():
    util = Util.Utils()
    logic = Logic.Logics()

    for i in range(1, 30):
        FILE_NAME_LIST.append(str(i))

    top_20_seq_dict = logic.get_top_N_seq_dict(FILE_NAME_LIST, WORK_DIR + CAS_OFF_EXCEL, TOP_N, PAM_N)

    # print(top_20_seq_dict)

    # for key ,val in top_20_seq_dict.items():
    #     print("key : " + key)
    #     print("val : " + str(val))

    cas_input_set = set()
    for val in top_20_seq_dict.values():
        for seq_str in val:
            cas_input_set.add(seq_str)
    print("cas_input_set len : " + str(len(cas_input_set)))
    util.make_txt_cas_off_finder_input(cas_input_set, INITIAL_CAS_OFF)

def make_cas_off_finder_input():
    util = Util.Utils()
    cas_input_set = set()
    for i in range(1, 30):
        FILE_NAME_LIST.append(str(i))

    total_len = 0
    for f_num in FILE_NAME_LIST:
        df_obj = util.read_excel_2_dataframe(WORK_DIR + CAS_OFF_EXCEL, f_num)

        data_obj = df_obj.to_dict()
        query_seq_dict = data_obj["order sgRNA Target sequence"]

        query_seq_dict_len = len(query_seq_dict)
        total_len += query_seq_dict_len
        print("file [" + str(f_num) + "] len : " + str(query_seq_dict_len))

        for query_idx in range(query_seq_dict_len):
            cas_input_set.add(query_seq_dict[query_idx] + PAM_N)

    util.make_txt_cas_off_finder_input(cas_input_set, INITIAL_CAS_OFF)

def make_random_cas_off_finder_input():
    cas_input_set = set()
    with open(WORK_DIR + CAS_OFF_TXT, "r") as f:
        while True:
            tmp_str = f.readline().replace("\n","")
            if tmp_str == "":
                break
            cas_input_set.add(tmp_str.split("\t")[1])

    with open(WORK_DIR + CAS_OFF_TXT_TOP_20 + "18.txt", "a") as new_f:
        new_f.write(REF_SRV_PATH + "\n")
        new_f.write("N" * SPACER_LEN + PAM_SEQ[0] + "\n")
        for new_str in cas_input_set:
            new_f.write(new_str + "NNN " + str(MAX_MISMATCH) +"\n")


start_time = clock()
print("start >>>>>>>>>>>>>>>>>>")
make_random_cas_off_finder_input()
# make_cas_off_finder_input_top_20_CAS_score()
# make_cas_off_finder_input()
print("::::::::::: %.2f seconds ::::::::::::::" % (clock() - start_time))
