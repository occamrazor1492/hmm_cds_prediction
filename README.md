# hmm_cds_prediction
# Predict CDS of Paramecium Bursaria Chlorella virus 1 base on Hidden Markov Model in Pomegranate

### Introduction
#### Relevant background to data
In the previous literature review, Encephalitozoon intestinalis is my target organism, which has the minimum amino acid quantity among all known Eukaryota. But things went south during experiment, the after counting the ATCG frequency in CP001942 (E.intestianlis), I found that the frequency of G in CP001942 CDS file is larger than that in the complete sequence file. According to Michael Hall, the reason that the G counts are higher in the CDS sequences than the superset of the total genome is because the coding sequences are often on the complementary DNA strand. So I decide to target another organism, Paramecium bursaria Chlorella virus 1 (NC_000852.5). 

Paramecium bursaria Chlorella virus 1(PBCV-1) is the type member of the genus Chlorovirus (family Phycodnaviridae) that infects certain chlorella-like green algae from freshwater sources; these viruses are found throughout the world. The chlorovirus host algae are normally symbionts of aquatic protists and in that state are resistant to virus infection.[1] Because of the nature of eukaryotic viruses and the nature of certain prokaryotic viruses, Chlorophytavirus is, in a sense, a type of virus that is intermediate between eukaryotic viruses and bacteriophages. The classification status is very special. The PBCV-1 virus infects only the so-called "external symbiotic" "Chlorella parasites" and "Ciliococcus chlorella" that are not free of protozoa, but does not infect the "internal symbiotic" chlorella. That is, the virus cannot enter the body of protozoa. This kind of chlorella can only live in animals because once they are separated, they are immediately attacked by the virus and die, which may have an evolutionary significance.[2]
#### Relevant background to method

Using HMM to predict CDS is the method that most adopted, since it is a sequence to sequence problem. 


![problem statement](https://i.imgur.com/jU72B1A.png)


The state diagram is like the image below, coding/non-coding state could go back to itself or can go to the other state. 


![state diagram](https://i.imgur.com/KiK54qN.png)


One of a good application of this method is a software named GeneMarkS[3], using statistics-based method (ab initio), doesn’t rely on any previous knowledge or experience parameters, it only used the features of sequence itself; so, the speed is better than the other two. The advantage of it is when we are lake of experimental data of an organism, it is a good choice, but the disadvantages is that this method tend to do over perdition, the false positive would be relative higher than others.    

But the method I would use is a combined of the statistics-based method (unsupervised learning) and supervised learning. I get the descriptions of HMM supervised learning training from this book[4]. I plan to use maximum likelihood estimation method. Generally, just calculate amino acid and CDS frequency of the knowing sequence with its CDS, and then get the three parameters of HMM model. 


What's more, I would use Pomegranate[5] for my supervised learning part , Pomegranate is an open source machine learning package for probabilistic modelling in Python; Three widely used probabilistic models implemented in pomegranate are general mixture models, hidden Markov models, and Bayesian networks; Pomegranate tends to abstract away the complexities of training models from their definition, which let user could focus on choose the correct training model instead of understanding the complex algorithm, in other words, it is like a black box. This approach trivially enables many useful learning strategies, such as out-of-core learning, mini batch learning, and semi-supervised learning, all could be called by few lines of code. Additionally, since it was written in Cpython, it could perform multithreaded parallelism.     

#### Hypothesis
* How accuracy in recognizing coding region. 

$$accuracy = \dfrac{C}{P} = 0.80 $$

C = correct prediction P = total prediction

* Fail rate of the predition

$$f = \dfrac{F}{P} = 0.10 $$
f = fail recognize rate F = fail recognize count
P = total actual count

### Method
#### Experimental approach

#####Get transition matrix, emission probability matrix, init matrix using Maximum likelihood method

Base on the method in this book[4], the formula is

•	The transition probability aij is:
If we assume that the frequency of sample transit from state i to j from time t to t+1 is Aij
  
$$ \widehat {a_{ij}} = \dfrac {A_{ij}}{\sum ^{N}_{j=1}A_{ij}}, i = 1,2,..N; j = 1,2,...N$$ 

$\hat{a_{01}} = \dfrac {A_{01}} {A_{00} + A_{01}}$

$\hat{a_{10}} = \dfrac {A_{10}} {A_{10} + A_{11}}$

$\hat{a_{11}} = \dfrac {A_{11}} {A_{10} + A_{11}}$

$\hat{a_{00}} = \dfrac {A_{00}} {A_{00} + A_{01}}$

• The emission probability bij is:
If we let Bjk stands for state j with observation value k, so the probability for state j with observation k, 

$$ \widehat {B_j}(k) = \dfrac {B_{ij}}{\sum ^{M}_{k = 1}B_{jk}}, j = 1,2,..N; k = 1,2,...N$$ 

$sum_1 = b_{1a} + b_{1t} + b_{1c} + b_{1g}$

$sum_0 = b_{0a} + b_{0t} + b_{0c} + b_{0g}$

$b_{1a} = \dfrac {b_{1a}}{sum_1}$ $b_{1t} = \dfrac{b_{1t}}{sum_1}$

$b_{1c} = \dfrac{b_{1c}}{sum_1}$ $b_{1g} = \dfrac{b_{1g}}{sum_1}$

$b_{0a} = \dfrac{b_{0g}}{sum_0}$ $b_{0t} = \dfrac{b_{0t}}{sum_0}$

$b_{0c} = \dfrac{b_{0c}}{sum_0}$ $b_{0g} = \dfrac{b_{0g}}{sum_0}$

• Initial state probability $\pi_i$ , the estimation $\widehat{\pi_i}$ would be 

$$\widehat{\pi_i}  =  \dfrac{Q}{S}$$
Q = frequency of $q_i$ as initial state in whole sample

S = sample number

The program for this part in `get_cds_array.py` 

```python
label_array, emission_array = write_label_file('complete_NC_000852.5.fasta', 'cds_number.txt', 'new.txt', 'NC_000852')

a1_0, a0_1, a1_1, a0_0 = calculate_transition_probability(label_array)

b1_a, b1_t, b1_c, b1_g, b0_a, b0_t, b0_c, b0_g = calculate_emission_probability(label_array, emission_array)
```
Input file is the complete sequence file, and the cds number I get from NCBI summary, also with the accession number of this target (crawler would get information from NCBI to compare with the actual number as error checking).

```python
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
```

Output is three parameter like this

```txt
total cds zone:  802
total char count is :  330611
label length 330611
```

Transiotion matrix

| state | count  |
|-------|--------|
| 10    | 273    |
| 01    | 273    |
| 11    | 305999 |
| 00    | 24065  |

Emission matrix

| State | Count |
|-------|-------|
| 1a    | 90520 |
| 1t    | 90935 |
| 1c    | 62838 |
| 1g    | 61979 |
| 0a    | 8768  |
| 0t    | 8254  |
| 0c    | 3728  |
| 0g    | 3589  |

ATCG counter in CDS and non-CDS state program in `countATCG_CDS.py`:

```txt
cds.txt  :  Counter({'T': 90935, 'A': 90520, 'C': 62838, 'G': 61979})

sequence.txt  :  Counter({'A': 99288, 'T': 99189, 'C': 66566, 'G': 65568})

C :
whole length:  66566
CDS  62838
rate 0.20517056733883607
NON CDS:  3728
rate 0.15316980977032746
----
T :
whole length:  99189
CDS  90935
rate 0.29690928325148885
NON CDS:  8254
rate 0.33912650478655654
----
A :
whole length:  99288
CDS  90520
rate 0.29555427854978583
NON CDS:  8768
rate 0.36024487448128517
----
G :
whole length:  65568
CDS  61979
rate 0.20236587085988925
NON CDS:  3589
rate 0.1474588109618308
----
```
#####Using Pomegranate to do prediction 
I do a blast of Paramecium bursaria Chlorella virus 1(PBCV-1) NC_000852,and decide to use JF411744.1 as test object. 
I follow this tutorial to do model training[6] 

```python
total_len = 330611 # Sequence length

# Amino acid count in each state, non-CDS and CDS
d1 = DiscreteDistribution({'A': 0.30, 'C': 0.20, 'G': 0.20, 'T': 0.30})
d2 = DiscreteDistribution({'A': 0.36, 'C': 0.15, 'G': 0.15, 'T': 0.34})
# Two state
s1 = State( d1, name='cds' )
s2 = State( d2, name='no_cds' )


hmm = HiddenMarkovModel()

hmm.add_states(s1, s2)
# init probability
hmm.add_transition( hmm.start, s1, init1)
hmm.add_transition( hmm.start, s2, init0)

# transition probability
hmm.add_transition( s1, s1, a11)
hmm.add_transition( s1, s2, a10)
hmm.add_transition( s2, s1, a01)
hmm.add_transition( s2, s2, a00)

hmm.bake()

# Pseudocounts
for i in range(10):
    print(hmm.fit(seq, max_iterations=1, use_pseudocount=True ))
    
hmm_predictions = hmm.predict( seq, algorithm='viterbi' )[1:-1]

trans, ems = hmm.forward_backward( seq )
print (trans)

```

The output of pseudocunts smoothing training improvement is 

```txt
Training improvement: 1783.2389817383373

Training improvement: -1.8283026292920113e-07

Training improvement: 1.8748687580227852e-07

Training improvement: 2.2462336346507072e-07

Training improvement: -8.731149137020111e-10
```

For the GeneMarkS I just choose the virus parameters and get the output. And then using regular expression like I did for the NCBI summary file to get the CDS zone number. 

```txt
   Gene    Strand    LeftEnd    RightEnd       Gene     Class
    #                                         Length
    1        -         512        1063          552        1
    2        +        1822        2217          396        2
    3        +        2288        3094          807        1
    4        -        3292        4701         1410        1
    5        +        4998        5735          738        1
    6        +        5768        6973         1206        1
    7        -        6970        8181         1212        1
```
```
512
1063
1822
2217
2288
3094
3292
4701
4998
5735
```







### Result

The results support my first hypothesis, the sensitivity rate is about 0.87. But I the CDS I fail to recognized is too huge, which doesn't support my hypothesis.  


*Chart*

In this chart we could find that my prediction miss predict some part with the actual one (in red circle), and GeneMarkS prediction just miss the very beginning part

![Compare 3 CDS](https://i.imgur.com/dJg3Zri.png)

This is another view of the same information, I use this formula `LOG(IFS(A1<340,340, A1<3400, 3400, A1< 34000, 34000, A1<340000, 340000))`; by this step, we could have a better visualization. But his image is only for my prediction and the actual CDS, not include GeneMarkS prediction. The grey point stands for miss predict that not exist in my prediction but exist in actual, golden point stands for those exist in my prediction but not in actual. All CDS I predict is 114269 and the actual CDS is 313522.    

![All index comparesion](https://i.imgur.com/9293WXJ.png)

This image below is to show from index 0 to 1000 what is the gap looks like

![Index < 1000](https://i.imgur.com/DXjRYIu.png)

index (1001-2000)
![Index (1001,2000)](https://i.imgur.com/XuEElik.png)

index (2000-3000)
![Index (2000,3000)](https://i.imgur.com/HfTcaXb.png)



*Analysis of Results*

* How accuracy in recognizing coding region. 

$$accuracy = \dfrac{C}{P} = \dfrac{99446}{114269} = 0.87 $$

C = correct prediction

P = total prediction

* Fail rate of the predition

$$f = \dfrac{F}{P} = \dfrac{313522-99446}{313522} = 0.68$$
f = fail to recognize rate

F = fail recognize count

P = total actual count



#### Further analysis 

For further analysis, first of all, choose another organism whose CDS has more gap, maybe this would be more suitable for using HMM. Secondly, using longer wordlist instead of only one amino acid may higher my accuracy. What's more, the visualization of result could be better with some tools like d3.js. Draw out the predicted CDS and actual CDS parallel, colour out those difference instead of putting them in a coordinate system. 









#### Reference
1. Cerny, R., Bauman, A. T., Roach, J. C., Lane, L. C., Agarkova, I. V., Wulser, K. W., ... & Jeannard, A. (2012). Paramecium bursaria Chlorella Virus 1 Proteome Reveals Novel Architectural and Regulatory Features of a Giant Virus.
2. Van Etten, J. L., Burbank, D. E., Kuczmarski, D., & Meints, R. H. (1983). Virus infection of culturable chlorella-like algae and dlevelopment of a plaque assay. Science, 219(4587), 994-996.
3. Lukashin, A. V., & Borodovsky, M. (1998). GeneMark. hmm: new solutions for gene finding. Nucleic acids research, 26(4), 1107-1115.
4. R. Durbin, S.Eddy, A.Krogh, G.Mitchison, "Parameter estimation for HMMs," in Biological Sequence Analysis, Cambridge, Cambridge University Press 1998, 2002, pp. 72-80.
5. Schreiber, J. (2017). Pomegranate: fast and flexible probabilistic modeling in python. arXiv preprint arXiv:1711.00137.
6. Jacob Schreiber (2017, Nov 30). Article title. Retrieved from https://github.com/jmschrei/pomegranate/blob/master/tutorials/Tutorial_3_Hidden_Markov_Models.ipynb
