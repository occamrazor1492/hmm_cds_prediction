from pomegranate import *
import numpy as np

with open('JF411744.1.fasta', 'r') as myfile:
    data=myfile.read().replace('\n', '')

seq = list(data)

# seq = list('CGACTACTGACTACTCGCCGACGCGACTGCCGTCTATACTGCGCATACGGC')
total_len = 330611
d1 = DiscreteDistribution({'A': 0.30, 'C': 0.20, 'G': 0.20, 'T': 0.30})
d2 = DiscreteDistribution({'A': 0.36, 'C': 0.15, 'G': 0.15, 'T': 0.34})

s1 = State( d1, name='cds' )
s2 = State( d2, name='no_cds' )

# gmm = GeneralMixtureModel( [d1, d2] )
hmm = HiddenMarkovModel()

hmm.add_states(s1, s2)
# init , maybe is not 0.5 for each
hmm.add_transition( hmm.start, s1, 0.5 )
hmm.add_transition( hmm.start, s2, 0.5 )

hmm.add_transition( s1, s1, 305999/total_len )
hmm.add_transition( s1, s2, 273/total_len )
hmm.add_transition( s2, s1, 273/total_len )
hmm.add_transition( s2, s2, 24065/total_len )

hmm.bake()

# Inertia
# print(hmm.fit( seq, distribution_inertia=0.3, edge_inertia=0.25 ))
# Pseudocounts
for i in range(5):
    print(hmm.fit(seq, max_iterations=1, use_pseudocount=True ))
# print(hmm.to_json())

# gmm_predictions = gmm.predict( np.array(seq) )
hmm_predictions = hmm.predict( seq, algorithm='viterbi' )[1:-1]
# hmm_predictions = hmm.predict(seq)

# cds_index = []
# for k in range(1,len(hmm_predictions)):
#     if hmm_predictions[k-1] == 0 and hmm_predictions[k] == 1:
#         cds_index.append(k)
#
#     elif hmm_predictions[k-1] == 1 and hmm_predictions[k] == 0:
#         cds_index.append(k)
#
# for x in range(0, len(cds_index)-1 ,2):
#     print((cds_index[x+1] - cds_index[x])/2 + cds_index[x])


# fhmm = open('hmm_predictions.txt','w')
# fhmm.write(( ''.join( map( str, hmm_predictions ) )))

# print ("sequence: {}".format( ''.join( seq ) ))
# print ("gmm pred: {}".format( ''.join( map( str, gmm_predictions ) ) ))
# print ("hmm pred: {}".format( ''.join( map( str, hmm_predictions ) ) ))


# print ("hmm state 0: {}".format( hmm.states[0].name ))
# print ("hmm state 1: {}".format( hmm.states[1].name ))

#
# # print(hmm.predict_log_proba(seq))
#
# # print (hmm.predict_proba( seq ))
#
trans, ems = hmm.forward_backward( seq )
print (trans)


# v_logp, v_seq = hmm.viterbi( seq )
# m_logp, m_seq = hmm.maximum_a_posteriori( seq )
