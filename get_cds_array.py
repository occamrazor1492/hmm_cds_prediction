import re
from bs4 import BeautifulSoup
from urllib.request import urlopen
from selenium import webdriver
import time
import functools
from pomegranate import *

dna = ['A','T','C','G']

# deal with the raw cds_xxxx.txt file which contain info about CDS
def get_cds_string(string):

    m = re.search('\d+\.\.\d+(?<!=(.gbkey=CDS.))', string)
    if m:
        found = m.group(0)
        return found
    # print("error in regx")


def get_cds_range(path):
    cds_range_two_d_array = []
    num_lines = sum(1 for line in open(path))
    # print("line number", num_lines)
    count = 0 # line number
    with open(path, "r") as ins:
        for line in ins:
            count += 1
            cds_range = []
            found = get_cds_string(line)
            # cds_range[0]: start  cds_range[1]: end
            try:
                cds_range = list(map(int, sorted(found.split('..'))))
                cds_range_two_d_array.append(cds_range)
            except Exception as e:
                pass
    print("total cds zone: ", len(cds_range_two_d_array))
    # print(cds_range_two_d_array)
    return cds_range_two_d_array

def get_sequence_length(genbank_id):
    url = 'https://www.ncbi.nlm.nih.gov/nuccore/'+str(genbank_id)
    # page = urlopen(url).read()
    # soup = BeautifulSoup(page, "html.parser")
    # # result = soup.findAll("pre", {"class": "genbank"})
    # # result = soup.get_text()
    # print(soup)
    # WebDriver driver = new RemoteWebDriver("http://localhost:9515", DesiredCapabilities.chrome());
    driver = webdriver.Chrome()
    driver.get(url)
    time.sleep(200)
    page = driver.page_source
    soup = BeautifulSoup(page, "html.parser")
    result = soup.findAll("pre", {"class": "genbank"})
    result = soup.get_text()
    m = re.search('(?!(LOCUS.*))\d+(?=\sbp)', result)
    if m:
        found = m.group(0)
        return int(found)
    print("error in regx")

def in_range(item, list):
    # list[0] is start, list[1] is End
    if item in range(list[0], list[1]+1):
        # print("yes")
        return True
    else:
        return False

# test
# data = get_cds_range('cds_cp001942.1.txt')
# print([in_range(898989, a_list) for a_list in data])

def write_label_file(fasta_file, cds_file, label_file, genbank_id):
    try:
        os.remove("temp.txt")
        os.remove("label_array.txt")
    except:
        pass
    # label array
    label_array = []
    emission_array = []
    flag = False
    # get sequence lenght from ncbi
    # sequence_lenght = get_sequence_length(genbank_id)
    # print("sequence length: ", sequence_lenght)
    # get cds zon [[],[],[]]
    cds_range = get_cds_range(cds_file)
    # delete first line
    with open(fasta_file, 'r') as fin:
        data = fin.read().splitlines(True)
    with open('temp.txt', 'w') as fout:
        fout.writelines(data[1:])

    with open('temp.txt') as f:
        word_count = 0
        while True:
            flag=False
            c = f.read(1)
            if (c in dna):
                word_count += 1
                emission_array.append(c)
                for a_range in cds_range:
                    if(in_range(word_count,a_range)):
                        label_array.append(1)
                        flag=True
                        break
                if(not flag):
                    label_array.append(0)
            if not c:
                print("End of file")
                break

    print ("total char count is : ", word_count)
    str1 = ''.join(str(e) for e in label_array)
    text_file = open("label_array.txt", "w")
    text_file.write(str1)
    text_file.close()
    print("label array length", len(label_array))
    return label_array, emission_array


def calculate_transition_probability(label_array):
    a1_0 = 0 # frequncy of 1 to 0 in state
    a0_1 = 0 # frequency of 0 to 1 in state
    a1_1 = 0 # 1 to 1
    a0_0 = 0 # 0 to 0
    for i in range(0, len(label_array)):
        if i+1 == len(label_array):
            break
        if label_array[i] == 1 and label_array[i+1] == 0:
            a1_0 += 1
        if label_array[i] == 0 and label_array[i+1] == 1:
            a0_1 += 1
        if label_array[i] == 1 and label_array[i+1] == 1:
            a1_1 += 1
        if label_array[i] == 0 and label_array[i+1] == 0:
            a0_0 += 1
    print("10:", a1_0, "01:", a0_1, "11:" ,a1_1, "00:", a0_0)
    return a1_0, a0_1, a1_1, a0_0

def cat_label_emission(label, emission):
    return str(label)+emission


def calculate_emission_probability(label_array, emission_array):
    b1_a = 0
    b1_t = 0
    b1_c = 0
    b1_g = 0

    b0_a = 0
    b0_t = 0
    b0_c = 0
    b0_g = 0

    if (len(label_array) != len(emission_array)):
        print("len error!")

    for i in range(0, len(label_array)):
        catString = cat_label_emission(label_array[i], emission_array[i])
        if catString == '1A':
            b1_a += 1
        elif catString == '1T':
            b1_t += 1
        elif catString == '1C':
            b1_c += 1
        elif catString == '1G':
            b1_g += 1
        elif catString == '0A':
            b0_a += 1
        elif catString == '0T':
            b0_t += 1
        elif catString == '0C':
            b0_c += 1
        else:
            b0_g += 1
    print('1a:',b1_a,'1t:',b1_t,'1c:', b1_c,'1g:', b1_g,'0a:', b0_a,'0t:', b0_t,'0c:', b0_c,'0g:', b0_g)

    return b1_a, b1_t, b1_c, b1_g, b0_a, b0_t, b0_c, b0_g


label_array, emission_array = write_label_file('CP001942.1.fasta', 'cds_number.txt', 'new.txt', 'cp001942')
a1_0, a0_1, a1_1, a0_0 = calculate_transition_probability(label_array)
b1_a, b1_t, b1_c, b1_g, b0_a, b0_t, b0_c, b0_g = calculate_emission_probability(label_array, emission_array)

# state transition
a01 = a0_1/(a0_0 + a0_1)
a10 = a1_0/(a1_0 + a1_1)

print('a01:',a01, 'a10', a10)


# emission
sum1 = b1_a + b1_t + b1_c + b1_g
sum0 = b0_a + b0_t + b0_c + b0_g

b1a = b1_a/sum1
b1t = b1_t/sum1
b1c = b1_c/sum1
b1g = b1_g/sum1

b0a = b0_a/sum0
b0t = b0_t/sum0
b0c = b0_c/sum0
b0g = b0_g/sum0

print("emission :   ",'1a:',b1a,'1t:',b1t,'1c:', b1c,'1g:', b1g,'0a:', b0a,'0t:', b0t,'0c:', b0c,'0g:', b0g)
model = HiddenMarkovModel()
