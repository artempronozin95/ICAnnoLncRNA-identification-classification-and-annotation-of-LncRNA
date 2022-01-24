import pandas as pd
import sys
import glob

def Specificity(x,y):
    if x+y == 0:
       return float(0)
    else:
       return x/(x+y)
def Precision(x,y):
    if x+y ==0:
       return float(0)
    else:
       return x/(x+y)
def Sensitivity(x,y):
    if x+y == 0:
       return float(0)
    else:
       return x/(x+y)
def F1(x,y):
    if x+y == 0:
       return float(0)
    else:
       return 2*(x*y)/(x+y)
def TP_FP(q_c,q_n,t_c,t_n):
    TP = 0
    FP = 0
    FN = 0
    TN = 0
    for w in t_n.to_list():
        if w in q_n.to_list():
           TP = TP+1
        if w in q_c.to_list():
           FP = FP+1
    for w in t_c.to_list():
        if w in q_c.to_list():
           TN = TN+1
        if w in q_n.to_list():
           FN = FN+1
    print(TP , FN , FP, TN)
    Pre = Precision(TP,FP)
    Spe = Specificity(TN,FP)
    Sen = Sensitivity(TP,FN)
    F = F1(Pre,Sen)
    return Pre , Sen, Spe, F

dirr = glob.glob(sys.argv[1] + '/*/')
f1_comp = {}
for w in dirr:
    cds_lnc = pd.read_csv(str(w) + '/compare.csv', sep='\t', header=None)
    test = pd.read_csv(str(w) + '/results.csv', sep=',', header=None)
    test = test[[0,1]]
    results = open(str(w) + '/F1.csv', 'w', encoding='utf-8')
    cds_lnc_cod = cds_lnc[cds_lnc[1].str.contains("mrna")]
    cds_lnc_non = cds_lnc[cds_lnc[1].str.contains("lnc")]
    test_cod = test[test[1].str.contains("NonCoding")==False]
    test_non = test[test[1].str.contains("NonCoding")]
    print(len(cds_lnc_cod[0]), len(cds_lnc_non[0]), len(test_cod[0]), len(test_non[0]))
    pre, sen, spe, f = TP_FP(cds_lnc_cod[0], cds_lnc_non[0], test_cod[0], test_non[0])
    f1_comp[w] = f
    print(pre, sen, spe, f, sep='\t' , file=results)
best_model = open(sys.argv[1] + '/best_model.txt', 'w', encoding='utf-8')
print(max(f1_comp, key=f1_comp.get), file=best_model)