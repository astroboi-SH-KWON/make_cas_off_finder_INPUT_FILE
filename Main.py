from time import clock

import Util

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
# CAS_OFF_TXT = "CAS_OFF_FINDER/INPUT/chlorochebus_sabaeus_cas_NGG_off_"
CAS_OFF_TXT_TOP_20 = "CAS_OFF_FINDER/INPUT/chlorochebus_sabaeus_top20_cas_NGG_off_"
MAX_MISMATCH = 3
SUB_SET_NUM = 100
REF_SRV_PATH = "Input/FASTA/chlorocebus_sabaeus_chr"
PAM_N = "NNN"

# INITIAL_CAS_OFF = ['NGG', SPACER_LEN, MAX_MISMATCH, SUB_SET_NUM, WORK_DIR + CAS_OFF_TXT, REF_SRV_PATH]
INITIAL_CAS_OFF = ['NGG', SPACER_LEN, MAX_MISMATCH, 18, WORK_DIR + CAS_OFF_TXT_TOP_20, REF_SRV_PATH]

########### MERGE OFF_TAGET RESULT #############
CNT_RESULT_PATH = "CAS_OFF_FINDER/Output/CountResult/chlorochebus_sabaeus_top20_cas_NGG_off_"
CNT_RESULT_EXT = "_result_count.txt"
FILE_CNT = 18
# FILE_CNT = SUB_SET_NUM

INITIAL_MRGE_OFF_TRGT = [WORK_DIR + CNT_RESULT_PATH, CNT_RESULT_EXT, FILE_CNT, PAM_N]
############### end setting env ################

def make_cas_off_finder_input_top_20_CAS_score():
    util = Util.Utils()
    for i in range(1,30):
        FILE_NAME_LIST.append(str(i))

    total_len = 0
    total_trnscrpt_cnt = 0
    top_20_seq_dict = {}
    total_seq_cnt = 0
    for f_num in FILE_NAME_LIST:
        df_obj = util.read_excel_2_dataframe(WORK_DIR + CAS_OFF_EXCEL, f_num)

        data_obj = df_obj.to_dict()
        transcript_id_dict = data_obj['Ensembl transcript ID']
        query_seq_dict = data_obj["order sgRNA Target sequence"]

        query_seq_dict_len = len(query_seq_dict)
        total_len += query_seq_dict_len
        print("file [" + str(f_num) + "] len : " + str(query_seq_dict_len))

        for query_idx in range(query_seq_dict_len):
            trnscrpt_id = transcript_id_dict[query_idx]
            if trnscrpt_id in top_20_seq_dict:
                if len(top_20_seq_dict[trnscrpt_id]) == 20:
                    continue
                else:
                    total_seq_cnt += 1
                    top_20_seq_dict[trnscrpt_id].append(query_seq_dict[query_idx] + PAM_N)
            else:
                total_trnscrpt_cnt += 1
                total_seq_cnt += 1
                top_20_seq_dict.update({trnscrpt_id: [query_seq_dict[query_idx] + PAM_N]})

    print("total_trnscrpt_cnt : " + str(total_trnscrpt_cnt))
    print("total_seq_cnt : " + str(total_seq_cnt))
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

start_time = clock()
print("start >>>>>>>>>>>>>>>>>>")
make_cas_off_finder_input_top_20_CAS_score()
# make_cas_off_finder_input()
print("::::::::::: %.2f seconds ::::::::::::::" % (clock() - start_time))
