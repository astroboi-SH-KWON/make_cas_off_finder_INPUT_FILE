

import Util

class Logics:
    def __init__(self):
        self.tmp = ""

    def get_top_N_seq_dict(self, f_name_list, path, top_n, pam_N):
        util = Util.Utils()

        total_len = 0
        total_trnscrpt_cnt = 0
        top_N_seq_dict = {}
        total_seq_cnt = 0
        for f_num in f_name_list:
            df_obj = util.read_excel_2_dataframe(path, f_num)

            data_obj = df_obj.to_dict()
            transcript_id_dict = data_obj['Ensembl transcript ID']
            query_seq_dict = data_obj["order sgRNA Target sequence"]

            query_seq_dict_len = len(query_seq_dict)
            total_len += query_seq_dict_len
            print("file [" + str(f_num) + "] len : " + str(query_seq_dict_len))

            for query_idx in range(query_seq_dict_len):
                trnscrpt_id = transcript_id_dict[query_idx]
                if trnscrpt_id in top_N_seq_dict:
                    if len(top_N_seq_dict[trnscrpt_id]) == top_n:
                        continue
                    else:
                        total_seq_cnt += 1
                        top_N_seq_dict[trnscrpt_id].append(query_seq_dict[query_idx] + pam_N)
                else:
                    total_trnscrpt_cnt += 1
                    total_seq_cnt += 1
                    top_N_seq_dict.update({trnscrpt_id: [query_seq_dict[query_idx] + pam_N]})

        print("total_trnscrpt_cnt : " + str(total_trnscrpt_cnt))
        print("total_seq_cnt : " + str(total_seq_cnt))
        return top_N_seq_dict