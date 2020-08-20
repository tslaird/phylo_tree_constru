import re
genome_ids, = glob_wildcards("input_files/{genome}.fa")
CRP_ids = ["L2","L3","L4","L5","L6","L14","L16","L18","L22","L24","S3","S8","S10","S17","S19"]

rule all:
    input:
        expand("CRP_subset/{CRP_ids}.fa", CRP_ids = CRP_ids),
        expand("CRP_subset/{CRP_ids}.msa",CRP_ids = CRP_ids),
        expand("CRP_subset/{CRP_ids}.msa.trim",CRP_ids = CRP_ids),
        expand("CRP_subset/{genome_ids}.msa.concat",genome_ids = genome_ids ),
        expand("CRP_subset/{genome_ids}.msa.concat_edit",genome_ids = genome_ids ),
        "master_alignment_file",
        "master_alignment_file.treefile",
        "master_alignment_file.nexus"


rule get_proteins:
    output: "CRP_subset/{CRP_ids}.fa"
    threads: 1
    run:
        CRP_ids = ["L2","L3","L4","L5","L6","L14","L16","L18","L22","L24","S3","S8","S10","S17","S19"]
        all_files = os.listdir("input_files")
        for i in all_files:
            with open("input_files/"+i) as file:
                genomic_proteins = file.read()
                match = re.findall(">.+ribosomal_protein_"+wildcards.CRP_ids+"\D[\s\S]+?\n{2}", genomic_proteins)[0]
                with open("CRP_subset/"+wildcards.CRP_ids+".fa", 'a') as outfile:
                    outfile.write(match)

rule align_proteins:
    input: "CRP_subset/{CRP_ids}.fa"
    output:
        "CRP_subset/{CRP_ids}.msa",
    shell:
        ''' muscle -in {input} -out {output} '''

rule trim_proteins:
    input: "CRP_subset/{CRP_ids}.msa"
    output:
        "CRP_subset/{CRP_ids}.msa.trim",
    shell:
        ''' trimal -in {input} -automated1 > {output} '''

rule trim_done:
    input: expand("CRP_subset/{CRP_ids}.msa", CRP_ids= CRP_ids)
    output: 'done.txt'
    shell:
        '''touch done.txt'''

rule concatenate_proteins:
    input:
        rules.trim_done.output
    output:
        "CRP_subset/{genome_ids}.msa.concat",
    run:
        CRP_ids = ["L2","L3","L4","L5","L6","L14","L16","L18","L22","L24","S3","S8","S10","S17","S19"]
        for i in CRP_ids:
            with open("CRP_subset/"+i+".msa.trim") as file:
                CRP_proteins = file.read()
                query= wildcards.genome_ids.split(".")[0]
                match = re.findall(r">.+"+query+r"[\s\S]+?\n?(?=\n>|$)", CRP_proteins)[0].split('\n',1)[1]
                if not os.path.exists("CRP_subset/"+wildcards.genome_ids+".msa.concat"):
                    with open("CRP_subset/"+wildcards.genome_ids+".msa.concat", 'w') as outfile:
                        outfile.write(">"+wildcards.genome_ids+"_concatenated_ribosomal_proteins\n")
                with open("CRP_subset/"+wildcards.genome_ids+".msa.concat", 'a') as outfile:
                    outfile.write(match)

rule edit_concatenation:
    input: "CRP_subset/{genome_ids}.msa.concat"
    output: "CRP_subset/{genome_ids}.msa.concat_edit"
    run:
        with open(str(input)) as file:
            file_text = file.read()
        header= file_text.split('\n')[0]
        sequence = ''.join(file_text.split('\n')[1::])
        out_text= re.sub('\n','', sequence)
        out_text='\n'.join([out_text[i:i+80] for i in range(0, len(out_text), 80)])
        out_text= '\n'.join([header, out_text])+"\n"
        with open(str(output), 'w') as outfile:
            outfile.write(out_text)

rule make_master_alignment:
    input: expand("CRP_subset/{genome_ids}.msa.concat_edit",genome_ids = genome_ids )
    output: "master_alignment_file"
    shell:
        ''' cat {input} >> master_alignment_file '''

rule make_tree:
    input: "master_alignment_file"
    output: "master_alignment_file.treefile"
    shell:
        ''' iqtree -s {input} '''

rule newick2nexus:
    input: "master_alignment_file.treefile"
    output: "master_alignment_file.nexus"
    run:
        with open(str(input)) as file:
            file_text=file.read()
        all_names= re.findall("(?<=\(|,)\w.+?(?=:)", file_text)
        out_str='#NEXUS\nbegin trees;\n\ttranslate\n'
        for i in range(0,len(all_names)):
            file_text = re.sub(all_names[i],str(i+1),file_text)
            if i == max(range(0,len(all_names))):
                out_str+="\t\t"+str(i+1)+"\t"+"_".join(all_names[i].split('_')[0:2])+";\n"
            else:
                out_str+="\t\t"+str(i+1)+"\t"+"_".join(all_names[i].split('_')[0:2])+",\n"
        combined_text= out_str+"\t\ttree UNTITLED = "+file_text+"\nend;"
        print(combined_text)
        with open(str(output),'w') as outfile:
            outfile.write(combined_text)
