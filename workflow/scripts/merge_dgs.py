#!/usr/bin/python3

import pandas as pd
import pyranges as pr
import os


outdir = snakemake.params.outdir
samples = snakemake.params.samples
gtf = snakemake.params.gtf

df_samples = pd.read_csv(samples,sep='\t')
samples_exp = df_samples['experiment'].values.tolist()
samples_acc = df_samples['accession'].values.tolist()
samples_id = df_samples['id'].values.tolist()

col_names = ["seqnames1","start1","end1","width1","strand1","region_id1","seqnames2","start2",
            "end2","width2","strand2","region_id2","n_reads","dg_id","score","gene_id.A","gene_id.B",
            "gene_name.A","gene_name.B","gene_type.A","gene_type.B","ambig.A","ambig.B","ambig_list.A",
            "ambig_list.B","cis","p_val","gene_count.A","gene_count.B","samples_id"]

df_result = pd.DataFrame(columns=col_names)

for i in range(len(samples_acc)):

    dg_file = os.path.join(outdir,samples_exp[i],samples_acc[i]+'_dg.tsv')
    df_dg = pd.read_csv(dg_file,sep='\t')
    df_dg['samples_id'] = samples_id[i]
    df_result = pd.concat([df_result,df_dg],ignore_index=True)

# If ambig list is not empty, separate them into different lines
df_result['ambig_list.A'] = df_result['ambig_list.A'].str.split(',')
df_result['ambig_list.B'] = df_result['ambig_list.B'].str.split(',')
df_result = df_result.explode('ambig_list.A',ignore_index=True)
df_result = df_result.explode('ambig_list.B',ignore_index=True)

# Reformat df by removing ambig columns
df_result["gene_name.A"] = df_result["ambig_list.A"]
df_result["gene_name.B"] = df_result["ambig_list.B"]
df_result = df_result[["seqnames1","start1","end1","width1","strand1","region_id1","seqnames2","start2",
                        "end2","width2","strand2","region_id2","n_reads","dg_id","score","gene_id.A","gene_id.B",
                        "gene_name.A","gene_name.B","gene_type.A","gene_type.B","cis","p_val","gene_count.A","gene_count.B","samples_id"]]

# Get gene info from gtf
df_gtf = pr.read_gtf(gtf).df
id_name_type = df_gtf[['gene_id','gene_name','gene_biotype']].drop_duplicates(ignore_index=True)

# Reannotate df_result (gene_id & gene_biotype)
for j in range(len(df_result)):

    name1 = df_result.iloc[j]['gene_name.A']
    name2 = df_result.iloc[j]['gene_name.B']
    
    if pd.notna(name1):
        id1 = id_name_type[id_name_type['gene_name']==name1].iloc[0]['gene_id']
        type1 = id_name_type[id_name_type['gene_name']==name1].iloc[0]['gene_biotype']
        df_result.loc[j,"gene_id.A"] = id1
        df_result.loc[j,"gene_type.A"] = type1
    
    if pd.notna(name2):
        id2 = id_name_type[id_name_type['gene_name']==name2].iloc[0]['gene_id']
        type2 = id_name_type[id_name_type['gene_name']==name2].iloc[0]['gene_biotype']
        df_result.loc[j,"gene_id.B"] = id2
        df_result.loc[j,"gene_type.B"] = type2

df_result.sort_values(by=['seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2'],inplace=True,ascending=[True,True,True,True,True,True])
df_result.to_csv(snakemake.output.dg,sep='\t',index=None)