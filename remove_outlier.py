import numpy as np

def monotonic(x):
    dx = np.diff(x)
    return np.all(dx <= 0) or np.all(dx >= 0)

def f7(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


fname = 'JF4117_CDS.txt'
with open(fname) as f:
    content0 = f.readlines()
# you may also want to remove whitespace characters like `\n` at the end of each line
content1 = [x.strip() for x in content0]
content = list(map(int, content1))

# actual_cds_all_index
actual_cds_all_index = []
for x in range(0, len(content)-1, 2):
    for index in range(content[x], content[x+1]):
        actual_cds_all_index.append(index)


# print all index out see how much different
# predict
predict_cds_all_index = []
fname1 = 'test_cds_index.txt'
with open(fname1) as f:
    predict_cds0 = f.readlines()
# you may also want to remove whitespace characters like `\n` at the end of each line
predict_cds1 = [x.strip() for x in predict_cds0]
predict_cds = list(map(int, predict_cds1))
for x in range(0, len(predict_cds)-1, 2):
    for index in range(predict_cds[x], predict_cds[x+1]):
        # print(index)
        predict_cds_all_index.append(index)

# not able to predict
print("predict", len(set(predict_cds_all_index)))
print("acutl", len(set(actual_cds_all_index)))

# for i in sorted(list(set(predict_cds_all_index))):
#     print(i)

# print difference
# difference = set(predict_cds_all_index)-set(actual_cds_all_index)
# for i in list(difference):
#     print(i)
# count = 0
# for i in set(predict_cds_all_index):
#     if i in set(actual_cds_all_index):
#         count += 1
#         print(count)
# print("count", count)

print((set(content)).intersection(set(predict_cds)))

# print(len(set(predict_cds_all_index).intersection(set(actual_cds_all_index))))

# print GeneMark.hmm cds_index
# genemark_cds_all_index = []
# fname1_1 = 'genemark.txt'
# with open(fname1_1) as f:
#     predict_cds0_1 = f.readlines()
# # you may also want to remove whitespace characters like `\n` at the end of each line
# predict_cds1_1 = [x.strip() for x in predict_cds0_1]
# predict_cds_1 = list(map(int, predict_cds1_1))
# for x in range(0, len(predict_cds_1)-1, 2):
#     for index in range(predict_cds_1[x], predict_cds_1[x+1]):
#         # print(index)
#         genemark_cds_all_index.append(index)
#
# for i in set(genemark_cds_all_index):
#     print(i)



# finalCDS = []
# print(len(content))
# def remove(content):
#     newCDS = []
#     for k in range(3, len(content)-1):
#         if k%2 == 1:
#             if content[k-1] < content[k]:
#                 newCDS.append(content[k-2])
#                 newCDS.append(content[k-1])
#                 newCDS.append(content[k])
#                 newCDS.append(content[k+1])
#     return f7(newCDS)
#
# finalCDS= remove(content)
# # print(len(finalCDS))
# while(monotonic(finalCDS) == False):
#     finalCDS = remove(finalCDS)
#
# print(len(finalCDS))
# for i in finalCDS:
#     print(i)
