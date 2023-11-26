
import numpy as np
import pandas as pd
import os,sys


# The main config file
configfile: "config/config.yaml"

drugs_we_care_about = config['drugs_we_care_about']
sir_drugs_we_care_about = ['SIR'+drug for drug in drugs_we_care_about]

model_name = ['raisin']

kmer_len = [25] #config['kmer_len_preprocessing']
nfeats = [1000]
model_types = ['xgb'] #['xgb','svm','ann']
n_top_feats = [10]

db_names = ['homolog', 'overexpression'] #['variant','overexpression','homolog']

rule all:
    input:
    # Combine the CARD files into one fasta to make one db
        #"external_tools/card-fastas/combined.fasta"
    # --------------------------------------------------------------------------
    # Make BLAST database
        #expand("external_tools/card-dbs/{db}_db/{db}_db.ndb",
        #    db=db_names),
    # --------------------------------------------------------------------------
    # Make BLAST query files
        #expand("results/wgs_standard/{k}mer/models/{g}/{p}/nested_cv/{f}feats/annotation/top_{ntf}_feats_blast.query",
        #    k=kmer_len,
        #    g=model_name,
        #    p=sir_drugs_we_care_about,
        #    f=nfeats,
        #    ntf=n_top_feats),
    # --------------------------------------------------------------------------
    # BLAST the queries against the dbs
        #expand("results/wgs_standard/{k}mer/top_feats/{g}/annotation_card/{db}_db/{p}/{p}_blast_output_raw.txt",
        #    k=kmer_len,
        #    g=model_name,
        #    p=sir_drugs_we_care_about,
        #    db=db_names),
    # --------------------------------------------------------------------------
    # Format the output and add card details
        #expand("results/wgs_standard/{k}mer/top_feats/{g}/annotation_card/{db}_db/{p}/{p}_annotated.tsv",
        #    k=kmer_len,
        #    g=model_name,
        #    p=sir_drugs_we_care_about,
        #    db=db_names),
    # --------------------------------------------------------------------------
    # Format the output and add card details
        expand("results/wgs_standard/{k}mer/top_feats/{g}/annotation_card/combined_db/{p}/{p}_annotated.tsv",
            k=kmer_len,
            g=model_name,
            p=sir_drugs_we_care_about),

rule reorder_card_fastas:
    output:
        o = "external_tools/card-fastas/combined.fasta"
    params:
        db1 = "external_tools/card-fastas/nucleotide_fasta_protein_homolog_model_variants.fasta",
        db2 = "external_tools/card-fastas/nucleotide_fasta_protein_overexpression_model_variants.fasta",
        db3 = "external_tools/card-fastas/nucleotide_fasta_protein_variant_model_variants.fasta"
    run:
        from Bio import Seq, SeqIO

        with open(output.o, 'a') as fout:
            for file in [params.db1, params.db2, params.db3]:
                with open(file) as fin:
                    for record in SeqIO.parse(fin, "fasta"):
                        contig_seq = record.seq
                        contig_id = record.id
                        contig_desc = record.description

                        #print(contig_seq)
                        #print(contig_id)
                        #print(contig_desc)


                        pieces = contig_desc.split("|")
                        reorder = [pieces[2],pieces[4],pieces[1],pieces[0],pieces[3]]
                        #print(pieces)
                        bloop = '|'.join(reorder)
                        #print(bloop)
                        record.description = bloop
                        record.id = bloop

                        SeqIO.write(record, fout, "fasta")





rule makeblastdb_card:
    input:
        "external_tools/card-fastas/combined.fasta"
    output:
        #"external_tools/card-dbs/{db}_db/{db}_db.ndb"
        "external_tools/card-dbs/combined_db/combined_db.ndb"
    params:
        #db_path = "external_tools/card-dbs/{db}_db/{db}_db"
        db_path = "external_tools/card-dbs/combined_db/combined_db"
    shell:
        "makeblastdb -in external_tools/card-fastas/combined.fasta -dbtype nucl -out external_tools/card-dbs/combined_db/combined_db"

'''        
"makeblastdb -in external_tools/card-fastas/nucleotide_fasta_protein_{wildcards.db}_model_variants.fasta -dbtype nucl -out external_tools/card-dbs/{wildcards.db}_db/{wildcards.db}_db"
'''

rule make_blast_query:
    input:
        imptfeats = "results/wgs_standard/{k}mer/models/{g}/{p}/nested_cv/{f}feats/xgb/top_{ntf}_feats.npy"
    output:
        q = "results/wgs_standard/{k}mer/models/{g}/{p}/nested_cv/{f}feats/annotation/top_{ntf}_feats_blast.query"
    run:
        output.q = str(output).split('.')[0]+'.query'
        top_feats = np.load(input.imptfeats, allow_pickle=True)
        if(isinstance(top_feats[0],bytes)):
            top_feats = [i.decode('utf-8') for i in top_feats]
        assert(len(top_feats[0])==int(wildcards.k))
        with open(output.q,'a') as fh:
            for feat in top_feats:
                fh.write(">{}\n".format(feat))
                fh.write(feat+"\n")


rule blast_top_feats:
    input:
        #db_flagfile = "external_tools/card-dbs/{db}_db/{db}_db.ndb",
        db_flagfile = "external_tools/card-dbs/combined_db/combined_db.ndb",
        query = expand("results/wgs_standard/{{k}}mer/models/{{g}}/{{p}}/nested_cv/{f}feats/annotation/top_{ntf}_feats_blast.query",f=[1000],ntf=[10])
    output:
        #"results/wgs_standard/{k}mer/top_feats/{g}/annotation_card/{db}_db/{p}/{p}_blast_output_raw.txt"
        "results/wgs_standard/{k}mer/top_feats/{g}/annotation_card/combined_db/{p}/{p}_blast_output_raw.txt"
    params:
        #db_path = "external_tools/card-dbs/{db}_db/{db}_db"
        db_path = "external_tools/card-dbs/combined_db/combined_db"
    shell:
        "blastn \
        -task blastn-short \
        -query {input.query} \
        -dust no \
        -word_size 25 \
        -max_target_seqs 50000 \
        -outfmt 6 \
        -out {output} \
        -db {params.db_path}"

#-word_size 25
#-ungapped
#-perc_identity 100

#-evalue 100000

rule format_blast_output:
    input:
        #f = "results/wgs_standard/{k}mer/top_feats/{g}/annotation_card/{db}_db/{p}/{p}_blast_output_raw.txt",
        f = "results/wgs_standard/{k}mer/top_feats/{g}/annotation_card/combined_db/{p}/{p}_blast_output_raw.txt",
        aro = "external_tools/card-ontology/aro.tsv"
    output:
        #o = "results/wgs_standard/{k}mer/top_feats/{g}/annotation_card/{db}_db/{p}/{p}_annotated.tsv",
        o = "results/wgs_standard/{k}mer/top_feats/{g}/annotation_card/combined_db/{p}/{p}_annotated.tsv",
    run:
        # Load input txt file (tab separated) and add the blast column names
        # blast names: https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
        blast_column_names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
        blast_df = pd.read_csv(input.f,header=None,sep='\t', names=blast_column_names)

        # if nothing was found, just save empty df
        if(blast_df.empty):
            print("no hits for drug, saving blank")
            blast_df.to_csv(output.o,sep="\t",index=False)
        else:
            # Insert column 0 with the drug name
            drug = wildcards.p.split("SIR")[-1]
            new_col = [drug]
            blast_df.insert(loc=0,column="antimicrobial",value=drug)

            # Extract "ARO:###" from column "sseqid" and put it in a new column called Accession
            # Few dont have an ARO number, so instead take the ARO_Name, and manually look it up later
            def extract_ARO(row):
                return row['sseqid'].split("|")[0]
                '''
                try:
                    return row['sseqid'].split("|")[2]
                except:
                    print(row['sseqid'])
                    return row['sseqid'].split("|")[-1]
                '''
            blast_df['Accession'] = blast_df.apply(lambda row: extract_ARO(row), axis=1)

            # Load the card table with descriptions
            desc_df = pd.read_csv(input.aro,sep='\t')

            # Combine the two dataframes on the "Accession" columns ("ARO:###")
            df = pd.concat([blast_df, desc_df], axis=1)
            df = pd.merge(blast_df,desc_df,how="left",on="Accession")

            # Reorder columns for easier reading
            column_order = ['antimicrobial', 'qseqid', 'Accession', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'ID', 'Name', 'CARD Short Name', 'Description', 'sseqid', ]

            df = df[column_order]
            #print(df)

            # Save
            df.to_csv(output.o,sep="\t",index=False)


