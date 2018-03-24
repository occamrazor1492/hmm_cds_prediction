from collections import Counter

def counterATCG(filename):
    with open(filename, "r") as myfile:
        data = myfile.read().replace('\n','')

        test = Counter(data)
        print(filename , " : ", test)

    myfile.close()
    return test

cds_dic = counterATCG("xxxcds_cp001942.1.txt")
no_cds_dic = counterATCG("temp.txt")

for k, v in cds_dic.items():
    print(v)
    print(no_cds_dic[k])
    temp = no_cds_dic[k]-v
    print(temp)
