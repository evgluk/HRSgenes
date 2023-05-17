import pandas
a=pd.read_csv("crossref.csv")
a['hhidpn']=1000*a['HHID'] + a['PN']
a_out =a.drop_duplicates(subset='hhidpn')
b=pd.read_csv("PGS3.csv")
b['hhidpn']=1000*b['HHID'] + b['PN']
b_out =b.drop_duplicates(subset='hhidpn')
ab = pd.merge(a_out,b_out, on='hhidpn', how='inner')
ab['IID']=ab['LOCAL_ID']
a_out['IID']=a_out['LOCAL_ID']
fam=pd.read_csv("[...]/targetDR.fam",sep = ' ',names=['FID', 'IID','POS1','POS2','G','PH'])
fam.loc[fam.apply(lambda x: 'NA' in x.IID,axis =1),"IID"]=np.nan
fam.loc[:,"IID"] = pd.to_numeric(fam['IID'])
eur = pd.merge(ab,fam, on='IID', how='inner')
eur_out = eur[['FID','IID']]
eur_out.to_csv('[...]/data/LDPred_sets123/EUR.valid.sample', index=False, sep=" ", header=None)
