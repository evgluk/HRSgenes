import pandas
omni = pd.read_csv('[...]/data/HumanOmni2.5-8v1_C.csv', sep =',', skiprows = 7)
omni = omni.drop(list(range(2379855,2379879)))
omnineg = omni[omni.apply(lambda x: '-' in x.RefStrand,axis =1)].copy()
flip_out = omnineg[['Name']]
flip_out.to_csv("[...]/data/fliplist.txt", index=False, header=False)

