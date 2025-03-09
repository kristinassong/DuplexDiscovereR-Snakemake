#!/usr/bin/python3

import pandas as pd
import pyranges as pr
from pybedtools import BedTool
import pybedtools
from multiprocessing import Pool
import ast
import sys
import os


outdir = sys.argv[2]
samples = sys.argv[1]
gtf = sys.argv[3]

# Set temporary directory for large files
pybedtools.set_tempdir(sys.argv[4])

df_samples = pd.read_csv(samples,sep='\t')
samples_exp = df_samples['experiment'].values.tolist()
samples_acc = df_samples['accession'].values.tolist()
samples_id = df_samples['id'].values.tolist()

col_names = ["seqnames1","start1","end1","width1","strand1","region_id1","seqnames2","start2",
            "end2","width2","strand2","region_id2","n_reads","dg_id","score","gene_id.A","gene_id.B",
            "gene_name.A","gene_name.B","gene_type.A","gene_type.B","ambig.A","ambig.B","ambig_list.A",
            "ambig_list.B","cis","p_val","gene_count.A","gene_count.B","samples_id"]

df_result = pd.DataFrame(columns=col_names)

# merge all samples and add sample id
for i in range(len(samples_acc)):

    dg_file = os.path.join(outdir,samples_exp[i],samples_acc[i]+'_dg.tsv')
    df_dg = pd.read_csv(dg_file,sep='\t')
    df_dg['samples_id'] = samples_id[i]
    df_result = pd.concat([df_result,df_dg],ignore_index=True)


###### Already well-annotated ######
df_result_anno1 = df_result[(df_result["ambig.A"]==0) & (df_result["ambig.B"]==0)]
df_result_anno2 = df_result[(df_result["ambig.A"].isna()) & (df_result["ambig.B"].isna())]
df_result_anno3 = df_result[(df_result["ambig.A"].isna()) & (df_result["ambig.B"]==0)]
df_result_anno4 = df_result[(df_result["ambig.A"]==0) & (df_result["ambig.B"].isna())]
df_result_uniq = pd.concat([df_result_anno1,df_result_anno2,df_result_anno3,df_result_anno4],ignore_index=True)


###### Handle rows containing ambig. genes ######
df_ambig = df_result[(df_result["ambig.A"]==1) | (df_result["ambig.B"]==1)].reset_index(drop=True)
df_ambig.reset_index(inplace=True)
df_ambig.rename(columns={"index":"org_index"},inplace=True)

# GTF
df_gtf = pr.read_gtf(gtf).df
df_gtf_exon = df_gtf[df_gtf['Feature']=='exon'][['Chromosome','Start','End','Strand','gene_id','gene_name','gene_biotype']] # reflect coco-corrected exons
df_gtf_transcript = df_gtf[df_gtf['Feature']=='transcript'][['Chromosome','Start','End','Strand','gene_id','gene_name','gene_biotype']]


def annotate(rna,gtf,type):

    pr_rna = pr.PyRanges(rna)
    pr_gtf = pr.PyRanges(gtf)
    
    annotated = pr_rna.join(pr_gtf,how='left',report_overlap=True).df
    annotated = annotated[['Chromosome','Start','End','Strand','Strand_b','org_index','ambig_list.'+type,'ambig.'+type,'gene_id','gene_name','gene_biotype','Overlap']].drop_duplicates(ignore_index=True)
    annotated = annotated[annotated['Strand'] == annotated['Strand_b']]
    annotated = annotated.groupby(by=['org_index','Overlap']).agg({
        'Chromosome':'first','Start':'first','End':'first','Strand':'first','ambig_list.'+type:'first','ambig.'+type:'first',
        'gene_id':';'.join,'gene_name':';'.join,'gene_biotype':';'.join
    })
    annotated.reset_index(inplace=True)
    annotated.org_index = annotated.org_index.astype(int)
    annotated.Overlap = annotated.Overlap.astype(int)
    annotated.sort_values(by=['org_index','Overlap'],inplace=True,ascending=[True,False])
    annotated.drop_duplicates(subset=['org_index'],keep="first",inplace=True)
    
    return annotated


def merge(df_org,df_annotated,type):
    # merge new annotations to original output
    df_org = df_org.merge(df_annotated,on='org_index',how='left')
    df_org.loc[df_org['ambig.'+type+'_x'] == 1, 'gene_id.'+type] = df_org['gene_id']
    df_org.loc[df_org['ambig.'+type+'_x'] == 1, 'gene_name.'+type] = df_org['gene_name']
    df_org.loc[df_org['ambig.'+type+'_x'] == 1, 'gene_type.'+type] = df_org['gene_biotype']
    return df_org


# ambig.genes in RNA1
df_ambigA = df_ambig[df_ambig['ambig.A']==1][["seqnames1","start1","end1","strand1","ambig_list.A",'ambig.A',"org_index"]]
df_ambigA.rename(columns={"seqnames1":"Chromosome","start1":"Start","end1":"End","strand1":"Strand"},inplace=True)
ambigA_annotated_run1 = annotate(df_ambigA,df_gtf_exon,'A') # annotate by coco-corrected exons

# second run for rows not annotated by coco-corrected exons
ambigA_annotated_run2 = ambigA_annotated_run1[ambigA_annotated_run1['Overlap']<0][['Chromosome','Start','End','Strand','ambig_list.A','ambig.A','org_index']]
ambigA_annotated_run2 = annotate(ambigA_annotated_run2,df_gtf_transcript,'A')
ambigA_annotated_run2 = ambigA_annotated_run2[ambigA_annotated_run2['Overlap']>=0] # rows annotated by transcripts

ambigA_annotated_run1 = ambigA_annotated_run1[ambigA_annotated_run1['Overlap']>=0] # rows annotated by coco-corrected exons
ambigA_annotated = pd.concat([ambigA_annotated_run1,ambigA_annotated_run2],ignore_index=True)

# merge new annotations to original output
df_ambig = merge(df_ambig,ambigA_annotated,'A')
df_ambig = df_ambig[['org_index','seqnames1','start1','end1','width1','strand1','region_id1','seqnames2','start2','end2','width2','strand2','region_id2','n_reads',
                     'dg_id','score','gene_id.A','gene_id.B','gene_name.A','gene_name.B','gene_type.A','gene_type.B','ambig.B','ambig_list.B','samples_id']]


# ambig.genes in RNA2
df_ambigB = df_ambig[df_ambig['ambig.B']==1][["seqnames2","start2","end2","strand2","ambig_list.B",'ambig.B',"org_index"]]
df_ambigB.rename(columns={"seqnames2":"Chromosome","start2":"Start","end2":"End","strand2":"Strand"},inplace=True)
ambigB_annotated_run1 = annotate(df_ambigB,df_gtf_exon,'B')

ambigB_annotated_run2 = ambigB_annotated_run1[ambigB_annotated_run1['Overlap']<0][['Chromosome','Start','End','Strand','ambig_list.B','ambig.B','org_index']]
ambigB_annotated_run2 = annotate(ambigB_annotated_run2,df_gtf_transcript,'B')
ambigB_annotated_run2 = ambigB_annotated_run2[ambigB_annotated_run2['Overlap']>=0]

ambigB_annotated_run1 = ambigB_annotated_run1[ambigB_annotated_run1['Overlap']>=0]
ambigB_annotated = pd.concat([ambigB_annotated_run1,ambigB_annotated_run2],ignore_index=True)

# merge new annotations to original output
df_ambig = merge(df_ambig,ambigB_annotated,'B')
df_ambig = df_ambig[['seqnames1','start1','end1','width1','strand1','region_id1','seqnames2','start2','end2','width2','strand2','region_id2','n_reads',
                     'dg_id','score','gene_id.A','gene_id.B','gene_name.A','gene_name.B','gene_type.A','gene_type.B','samples_id']]


###### Combine well-annotated with newly annotated interactions ######
df_result_uniq = df_result_uniq[['seqnames1','start1','end1','width1','strand1','region_id1','seqnames2','start2','end2','width2','strand2','region_id2','n_reads',
                                 'dg_id','score','gene_id.A','gene_id.B','gene_name.A','gene_name.B','gene_type.A','gene_type.B','samples_id']]
df_result = pd.concat([df_result_uniq,df_ambig],ignore_index=True)
df_result['dg_id'] = df_result['dg_id'].astype(str)
df_result['dg_id'] = df_result['samples_id']+'-'+df_result['dg_id']
df_result.drop(columns=['region_id1','region_id2','samples_id','score'],inplace=True)
df_result.to_csv(os.path.join(outdir,'all_dgs_before_merge.tsv'),sep='\t',index=None)


###### Merge RNA-RNA interactions ######
# Multiprocessing by chromosome1
chromo1_list = set(df_result['seqnames1'].values.tolist())
os.makedirs(os.path.join(outdir,'tmp'),exist_ok=True)

def merge_interactions(chromo1,df_result=df_result):
    print(chromo1)
    df_result = df_result[df_result['seqnames1']==chromo1] # select chromosome
    rna1_bed = BedTool.from_dataframe(df_result[['seqnames1','start1','end1','dg_id','n_reads','strand1']]).sort()
    rna1_merged = rna1_bed.merge(s=True,c=[4,5,6],o="collapse,sum,distinct",delim="|").to_dataframe()
    rna2_bed = BedTool.from_dataframe(df_result[['seqnames2','start2','end2','dg_id','n_reads','strand2']]).sort()
    rna2_merged = rna2_bed.merge(s=True,c=[4,5,6],o="collapse,sum,distinct",delim="|").to_dataframe()

    # create dict to store merged groups
    rna1_merged_groups = rna1_merged[rna1_merged['name'].str.contains('|', regex=False)]
    rna2_merged_groups = rna2_merged[rna2_merged['name'].str.contains('|', regex=False)]

    if len(rna1_merged_groups)>0 and len(rna2_merged_groups)>0:
        rna1_merged_groups_dict = {}
        for idx, row in rna1_merged_groups.iterrows():
            ids = row['name'].split('|')
            rna1_merged_groups_dict[f"group{idx}"] = ids

        rna2_merged_groups_dict = {}
        for idx, row in rna2_merged_groups.iterrows():
            ids = row['name'].split('|')
            for id in ids:
                rna2_merged_groups_dict[id] = f"group{idx}"

        df_to_merge = pd.DataFrame(columns=['id'])

        for key, value in rna1_merged_groups_dict.items():
            group_list = []
            for v in value:
                group_list.append(rna2_merged_groups_dict.get(v,'absent'))
            df = pd.DataFrame({'id':value,'rna2_merged_group':group_list})
            df = df[df['rna2_merged_group'] != 'absent']
            if len(df)>0:
                df = df.groupby("rna2_merged_group")['id'].apply('|'.join).reset_index()
                df = df[df['id'].str.contains('|', regex=False)].reset_index(drop=True)
                df = df[['id']]
                df_to_merge = pd.concat([df_to_merge,df],ignore_index=True)

        if len(df_to_merge)>0:
            df_to_merge['id'] = df_to_merge['id'].str.split('|')

            for i in df_to_merge.id.values.tolist():
                # i : list of dg ids to merge
                rows_extracted = pd.DataFrame(columns=df_result.columns)
                i = str(i)
                i = ast.literal_eval(i)
                for id in i:
                    rows_extracted = pd.concat([rows_extracted,df_result[df_result['dg_id']==id]],ignore_index=True)
                rows_extracted_rna1_bed = BedTool.from_dataframe(rows_extracted[['seqnames1','start1','end1','dg_id','n_reads','strand1']]).sort()
                rows_extracted_rna1_merged = rows_extracted_rna1_bed.merge(s=True,c=[4,5,6],o="collapse,sum,distinct",delim="|").to_dataframe()
                rows_extracted_rna2_bed = BedTool.from_dataframe(rows_extracted[['seqnames2','start2','end2','dg_id','n_reads','strand2']]).sort()
                rows_extracted_rna2_merged = rows_extracted_rna2_bed.merge(s=True,c=[4,5,6],o="collapse,sum,distinct",delim="|").to_dataframe()
                
                if len(rows_extracted_rna1_merged)>0 and len(rows_extracted_rna2_merged)>0:
                    validate1 = rows_extracted_rna1_merged[rows_extracted_rna1_merged['name'].str.contains('|', regex=False)]['name'].values.tolist()
                    validate1 = [tuple(s1.split("|")) for s1 in validate1]
                    validate2 = rows_extracted_rna2_merged[rows_extracted_rna2_merged['name'].str.contains('|', regex=False)]['name'].values.tolist()
                    validate2 = [tuple(s2.split("|")) for s2 in validate2]
                    if len(validate1)>0 and len(validate2)>0:
                        final = []
                        for v1 in validate1:
                            for v2 in validate2:
                                if len(set(v1) & set(v2))>1:
                                    final.append(frozenset(set(v1) & set(v2)))
                        
                        if len(final)>=1:
                            for ids in final:
                                final_rows_extracted = pd.DataFrame(columns=df_result.columns)
                                for el in ids:
                                    final_rows_extracted = pd.concat([final_rows_extracted,df_result[df_result['dg_id']==el]],ignore_index=True)
                                    df_result = df_result[df_result['dg_id']!=el] # remove rows to be merged from original data
                                final_rows_extracted_rna1_bed = BedTool.from_dataframe(final_rows_extracted[['seqnames1','start1','end1','dg_id','n_reads','strand1']]).sort()
                                final_rows_extracted_rna1_merged = final_rows_extracted_rna1_bed.merge(s=True,c=[4,5,6],o="collapse,sum,distinct",delim="|").to_dataframe()
                                final_rows_extracted_rna2_bed = BedTool.from_dataframe(final_rows_extracted[['seqnames2','start2','end2','dg_id','n_reads','strand2']]).sort()
                                final_rows_extracted_rna2_merged = final_rows_extracted_rna2_bed.merge(s=True,c=[4,5,6],o="collapse,sum,distinct",delim="|").to_dataframe()

                                new_merged_row = pd.DataFrame({'seqnames1':[final_rows_extracted_rna1_merged.loc[0,'chrom']],
                                                            'start1':[final_rows_extracted_rna1_merged.loc[0,'start']],
                                                            'end1':[final_rows_extracted_rna1_merged.loc[0,'end']],
                                                            'width1':[int(final_rows_extracted_rna1_merged.loc[0,'end'])-int(final_rows_extracted_rna1_merged.loc[0,'start'])+1],
                                                            'strand1':[final_rows_extracted_rna1_merged.loc[0,'strand']],
                                                            'seqnames2':[final_rows_extracted_rna2_merged.loc[0,'chrom']],
                                                            'start2':[final_rows_extracted_rna2_merged.loc[0,'start']],
                                                            'end2':[final_rows_extracted_rna2_merged.loc[0,'end']],
                                                            'width2':[int(final_rows_extracted_rna2_merged.loc[0,'end'])-int(final_rows_extracted_rna2_merged.loc[0,'start'])+1],
                                                            'strand2':[final_rows_extracted_rna2_merged.loc[0,'strand']],
                                                            'n_reads':[final_rows_extracted_rna2_merged.loc[0,'score']],
                                                            'dg_id':[final_rows_extracted_rna2_merged.loc[0,'name']],
                                                            'gene_id.A': [final_rows_extracted.loc[0,'gene_id.A']],
                                                            'gene_id.B': [final_rows_extracted.loc[0,'gene_id.B']],
                                                            'gene_name.A': [final_rows_extracted.loc[0,'gene_name.A']],
                                                            'gene_name.B': [final_rows_extracted.loc[0,'gene_name.B']],
                                                            'gene_type.A': [final_rows_extracted.loc[0,'gene_type.A']],
                                                            'gene_type.B': [final_rows_extracted.loc[0,'gene_type.B']]})

                                df_result = pd.concat([df_result,new_merged_row],ignore_index=True)
    df_result.to_csv(os.path.join(outdir,'tmp','chr'+str(chromo1)+'_tmp_all_dgs.tsv'),sep='\t',index=None)
    return df_result


with Pool(processes=len(chromo1_list)) as pool:
    results = pool.map(merge_interactions, chromo1_list)

df_final = pd.concat(results)

# final output
df_final.sort_values(by=['seqnames1', 'start1', 'end1', 'seqnames2', 'start2', 'end2'],inplace=True,ascending=[True,True,True,True,True,True])
df_final.to_csv(os.path.join(outdir,'all_dgs.tsv'),sep='\t',index=None)
