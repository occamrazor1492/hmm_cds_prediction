from collections import Counter

def counterATCG(filename):
    with open(filename, "r") as myfile:
        data = myfile.read().replace('\n','')

        test = Counter(data)
        print(filename , " : ", test)

    myfile.close()
    return test

cds_dic = counterATCG("my_own_cds.txt")
all_dic = counterATCG("temp.txt")
cds_count_sum = sum(cds_dic.values())
all_dic_sum = sum(all_dic.values())
non_cds_dic_sum = all_dic_sum-cds_count_sum

for k, v in cds_dic.items():
    print(k, ":")
    print("whole length: ",all_dic[k])
    print("CDS ",v)
    print("rate", v/cds_count_sum)
    temp = all_dic[k]-v
    print("NON CDS: ",temp)
    print("rate", temp/non_cds_dic_sum)
    print("----")
