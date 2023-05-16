import pandas as pd
# HRSnew-PLINK-EA2
ea2_snps = pd.read_csv('[...]/data/EA2_excl_23andMe_HRS.meta',sep='\t', skiprows=1, names=['SNP', 'CHR','BP','A1old','A2old', 'Freq','EAF_SE','EAF_min', 'EAF_max', 'N' ,'Z' , 'P','Direction' ,'HetISq' ,'HetChiSq' ,'HetDf' ,'HetPVal' ,'BETAold', 'SE'])
# case in which we have negative beta - flip
ea2_snps['A1'] = ea2_snps.A1old
ea2_snps.loc[ea2_snps.BETAold<0, "A1"] = ea2_snps.A2old[ea2_snps.BETAold<0]
ea2_snps['A2'] = ea2_snps.A2old
ea2_snps.loc[ea2_snps.BETAold<0, "A2"] = ea2_snps.A1old[ea2_snps.BETAold<0]
ea2_snps['BETA'] = ea2_snps.BETAold
ea2_snps.loc[ea2_snps.BETAold<0, "BETA"] = ea2_snps.BETAold[ea2_snps.BETAold<0]*(-1) 
cleaned_out2 = ea2_snps[['SNP','CHR','BP','A1','A2','BETA','P','SE']]
cleaned_out2.to_csv('[...]/data/cleaning/forprsice/EA2_HRS2_bf.txt', index=False, sep=" ")
ea2_snps_signeg = ea2_snps[(ea2_snps['P']<=5e-8) & (ea2_snps['BETAold']<0)]
ea2_snps_sigpos = ea2_snps[(ea2_snps['P']<=5e-8) & (ea2_snps['BETAold']>=0)]
ea2_snps_nsneg = ea2_snps[(ea2_snps['P']>5e-8) & (ea2_snps['BETAold']<0)]
ea2_snps_nspos = ea2_snps[(ea2_snps['P']>5e-8) & (ea2_snps['BETAold']>=0)]
cleaned_out_signeg = ea2_snps_signeg[['SNP','CHR','BP','A1','A2','BETA','P','SE']]
cleaned_out_signeg.to_csv('[...]/data/cleaning/forprsice/sigsplit/EA2_HRS2_bf_signeg.txt', index=False, sep=" ")
cleaned_out_sigpos = ea2_snps_sigpos[['SNP','CHR','BP','A1','A2','BETA','P','SE']]
cleaned_out_sigpos.to_csv('[...]/data/cleaning/forprsice/sigsplit/EA2_HRS2_bf_sigpos.txt', index=False, sep=" ")
cleaned_out_nsneg = ea2_snps_nsneg[['SNP','CHR','BP','A1','A2','BETA','P','SE']]
cleaned_out_nsneg.to_csv('[...]/data/cleaning/forprsice/sigsplit/EA2_HRS2_bf_nsneg.txt', index=False, sep=" ")
cleaned_out_nspos = ea2_snps_nspos[['SNP','CHR','BP','A1','A2','BETA','P','SE']]
cleaned_out_nspos.to_csv('[...]/data/cleaning/forprsice/sigsplit/EA2_HRS2_bf_nspos.txt', index=False, sep=" ")
