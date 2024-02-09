#! /usr/bin
# -*- coding: utf-8 -*-

#activate IGv2 conda environment

import sys 
import pysam
import os
import vcf
import glob
from itertools import islice
import random
import re
import json


def get_file_paths():
    # Get the directory where AIM_pipeline.py is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Construct paths for required files
    aim_positions = os.path.join(script_dir, "ALLAIMPostionshg38.txt")
    aimfile = os.path.join(script_dir, "Capture96AIMPos.txt")
    aim_string = os.path.join(script_dir, "Capture96AIMString.txt")
    aim_1kgp_sup_gt = os.path.join(script_dir, "1KGP_sup_AIM_GT.txt")
    aim_1kgp_sup_pop_code = os.path.join(script_dir, "1KGP_sup_Pop_code.txt")
    extraparams = os.path.join(script_dir, "extraparams")
    mainparams = os.path.join(script_dir, "mainparams")
    
    return aim_positions, aimfile, aim_string, aim_1kgp_sup_gt, aim_1kgp_sup_pop_code, extraparams, mainparams


def process_vcf_files(vcf_file,aim_positions_file,sample_name,sample_path):
    new_directory = os.path.join(sample_path,'{}_Ancestry'.format(sample_name))
    if not os.path.exists(new_directory):
        os.makedirs(new_directory)
    else:
        # Remove existing directory and create a new one
        os.rmdir(new_directory)
        os.makedirs(new_directory)
    os.chdir(new_directory)
    for i in range(1, 23):
        with open('AIMPosition{}'.format(i), 'w') as outfile:
            cmd1 = "awk -v i={} '$1 == i' {} | awk '{{ print \"chr\"$1\" \"$2 }}' > AIMPosition{}".format(i, aim_positions_file, i)
            os.system(cmd1)

            cmd2 = "vcftools --vcf {} --positions AIMPosition{} --recode --out SampleAIMInfo{}".format(vcf_file, i, i)
            os.system(cmd2)

def concatenate_vcf_files(directory_path, sample_name):
    input_files = [file for file in os.listdir(directory_path) if file.endswith('.recode.vcf')]
    output_file = "Sample_ALLAIM.vcf"

    writer = vcf.Writer(open(output_file, 'w'), vcf.Reader(open(os.path.join(directory_path, input_files[0]), 'r')))

    # Iterate through each input VCF file
    for file in input_files:
        reader = vcf.Reader(open(os.path.join(directory_path, file), 'r'))
        for record in reader:
            # Check if the record belongs to the desired sample
            if record.samples[0].sample == sample_name:
                # Write the record to the output file
                writer.write_record(record)

    # Close the output VCF file
    writer.close()
    [os.remove(file) for file in os.listdir(directory_path) if file != 'Sample_ALLAIM.vcf']


def read_aim_pos(aim_file):
    chr_no=[]
    aim_pos=[]
    all_aim_pos=[]
    with open(aim_file,'r') as f1:
        for line in f1:
            chr_no.append(line.rstrip().split('_')[0])
            aim_pos.append(line.rstrip().split('_')[1])
    for i,j in zip(chr_no,aim_pos):    
        all_aim_pos.append([i,j])
    return(all_aim_pos,aim_pos)


def get_aim_coverage(aim_str,aim_pos,bam_file):
    chrno=[]
    chrstartpos=[]
    chrendpos=[]
    with open (aim_str) as pos:
        for line in pos:
            lines=line.replace('\n','')
            chrno.append(lines.split('_')[0])
            chrstartpos.append(lines.split('_')[1])
            chrendpos.append(lines.split('_')[2])
    samfile = pysam.AlignmentFile(bam_file, "rb" )
    all_aim_with_any_cov=[]
    high_cov_column=[]
    my_coverage=[]
    for (a,b,c,d) in zip(chrno,chrstartpos,chrendpos,aim_pos):
        for pileupcolumn in samfile.pileup(a,int(b),int(c)):
            if  pileupcolumn.pos==int(d):
                all_aim_with_any_cov.append(pileupcolumn.pos)
                if pileupcolumn.n>=10:
                    high_cov_column.append(pileupcolumn.pos)
                    my_coverage.append(pileupcolumn.n)
    return(high_cov_column)


def read_Vcf(vcf_file):
    samplename=[]
    ChrPosGT={} # dict of pos and GT
    all_chr=[]
    all_pos=[]
    all_GT=[]
    with open(vcf_file,'r') as vcf:
        for line in vcf:
            if line.startswith('##'):
                continue
            elif line.startswith('#'):
                sampleinfo=line.rstrip().split('\t')
                samplename=sampleinfo[9:]                                 
            elif line.startswith('chr'):
                AllInfo=line.rstrip().split('\t')
                all_pos.append(AllInfo[1])
                ChrPos='_'.join(AllInfo[0:2])
                all_chr.append(ChrPos)
                all_GT.append(AllInfo[9].split(':')[0])
        ChrPosGT={str(y): str(z) for y,z in zip(all_pos,all_GT)}
    return(samplename,ChrPosGT)


def add_sample_code(ref_sample_no,samplename):
    sample_code=[int(ref_sample_no) + 1]
    sample_code_dict = dict(zip(samplename,sample_code))
    return(sample_code_dict)


def convert_to_structure(aim_pos,ChrPosGT,high_cov_column,sample_code_dict):
    AllGT=[]
    AllHap1=[]
    AllHap2=[]
    SamAllHap1=[]
    SamAllHap2=[]
    
    for n in aim_pos:
        if str(n) in map(str, high_cov_column):
            if str(n) in ChrPosGT:
                AllGT.append(ChrPosGT[str(n)])
            else:
                AllGT.append('0/0')
        else:
            AllGT.append('-9/-9')

    for n in AllGT:
        splittedN=n.split('/')
        AllHap1.append(splittedN[0])
        AllHap2.append(splittedN[1])
    for p in AllHap1:
        p=p.replace("1","2")
        p=p.replace("0","1")
        SamAllHap1.append(p)
    for p in AllHap2:
        p=p.replace("1","2")
        p=p.replace("0","1")
        SamAllHap2.append(p)
    for k,v in sample_code_dict.items():
        SamAllHap1.insert(0,str(v))
        SamAllHap2.insert(0,str(v))
    
    SamAllHap1 = SamAllHap1[:1] + ['1', '0'] + SamAllHap1[1:]
    SamAllHap2 = SamAllHap2[:1] + ['1', '0'] + SamAllHap2[1:]
    return(SamAllHap1,SamAllHap2)



def out_ref_GT(ref_GT,samplehap1,samplehap2,directory_path):
    newsam1=' '.join([str(elem) for elem in samplehap1])
    newsam2=' '.join([str(elem) for elem in samplehap2])
    with open(os.path.join(directory_path, '1KGP_SuperPop_with_SampleGT.txt'), 'w') as newfile:
        with open(ref_GT, 'r') as file:
            for line in file:
                newfile.write(line)
            newfile.write(newsam1)
            newfile.write('\n')
            newfile.write(newsam2)
            newfile.write('\n')


def run_structure(input_file, main_params, extra_params):
    command = 'structure -i {} -m {} -e {}'.format(input_file, main_params, extra_params)
    return os.system(command)

def structure_result_interpretation(myfile, PopNo, ClnoOld, PopCfile, SampleNo ,high_cov_column):
    InputPCode = []
    InputPop = []
    InputPopnew = []
    Clno = int(ClnoOld) + 1
    number = [str(i) for i in range(1, Clno)]
    mymatrix = []
    NoLines = int(PopNo) + 2
    PopCode = []
    Popvalue = []
    AllPopValue = []
    FinalPopValue = []
    PopDict = {}
    Finaldict = None
    ResultLine = None
    Mysampledict = {}

    with open(PopCfile, 'r') as codefile:
        for line in codefile:
            line = line.strip('\n')
            InputPCode.append(line.split('\t')[1])
            InputPop.append(line.split('\t')[0])
    InputPopnew = [InputPop for i in range(int(ClnoOld))]

    with open(myfile, 'r') as oldfile:
        for line in oldfile:
            line = line.rstrip()
            if "Given    Inferred Clusters                     Number of" in line:
                mymatrix.append(" ".join(islice(oldfile, NoLines)))
            if "        Label (%Miss) Pop" in line:
                continue
            if line.startswith('%s' % (int(SampleNo) + 1)):
                ResultLine = line
    ResultLine = ResultLine.strip()
    ResultLine = re.sub(r"\s+", ' ', ResultLine)
    ResultLine = ResultLine.split(':')[1]
    ResultLine = ResultLine.split(' ')[1:int(Clno)]
    ResultLine = [float(n) for n in ResultLine]

    for n in mymatrix: 
        char = n.split("\n")
    del (char[0:2])
    del (char[-1])
    for line in char:
        line = line.strip(' ')
        PopCode.append(line.split(':')[0])
        Popvalue.append(line.split(':')[1])
    for n in Popvalue:
        n = re.sub(r"\s+", ' ', n)
        n = n.split(' ')
        AllPopValue.append(n[1:int(Clno)])

    FinalPopValue = list(map(list, zip(*AllPopValue)))
    PopDict = {str(a): [b, c] for a, b, c in zip(number, InputPopnew, FinalPopValue)}

    for key, value in PopDict.items():
        PopDict[key] = dict(zip(* value))

    Finaldict = {key: dict(sorted(val.items(), key=lambda ele: ele[1], reverse=True)) for key, val in PopDict.items()}

    for count, value in enumerate(ResultLine, start=1):
        Mysampledict[count] = value
    Mysampledict = {k: v for k, v in sorted(Mysampledict.items(), key=lambda item: item[1], reverse=True)}
    
    max_ancestry = max(Mysampledict, key=Mysampledict.get)
    max_percentage = Mysampledict[max_ancestry]
    
    ancestry_labels = {
        "1": "European",
        "2": "American",
        "3": "African",
        "4": "East Asian",
        "5": "South Asian"
    }
    
    inferred_ancestry_label = ancestry_labels.get(max_ancestry, max_ancestry)

    with open("Sample_ancestry_result.txt", "w") as file:
        sys.stdout = file
        print("\nInferred ancestry: {}\n".format(ancestry_labels.get(str(inferred_ancestry_label), str(inferred_ancestry_label))))
        print("Likelihood of this sample belonging to different ancestries:")
        
        sorted_ancestries = sorted(Mysampledict.items(), key=lambda item: item[1], reverse=True)
        for k, v in sorted_ancestries:
            print("{}: {:.3f}".format(ancestry_labels.get(str(k), str(k)), v))
        
        num_aims_used = len(high_cov_column)
        num_aims_missing = 96 - num_aims_used
        
        print("\nNumber of AIMs used: {}/96 ".format(num_aims_used))
        print("Number of AIMs missing: {} \n".format(96 - num_aims_used))
    
        print("Reduced usage of AIMs reduces the accuracy of inferred ancestry.\n")
        sys.stdout = sys.__stdout__  # Reset stdout to its default value



def main():
    sample_name= sys.argv[1] #give an input necessary
    sample_path= sys.argv[2] #directory need to give where all sample are there
    aim_positions, aimfile, aim_str, ref_GT_file, PopCodeFILE, extra_params, main_params = get_file_paths()
    vcf_input = os.path.join(sample_path,'variants',"snvs_from_ccs.vcf")
    process_vcf_files(vcf_input,aim_positions,sample_name,sample_path)
    directory_path = os.path.join(sample_path,"{}_Ancestry".format(sample_name))
    concatenate_vcf_files(directory_path, sample_name)
    aim_pos_reading =read_aim_pos(aimfile)
    bam_file = "{}/alignments/ccs_to_ref_phased.sorted.bam".format(sample_path)
    vcf_file = "{}/{}_Ancestry/Sample_ALLAIM.vcf".format(sample_path,sample_name)
    ref_length=2504
    high_cov_pos=get_aim_coverage(aim_str,aim_pos_reading[1],bam_file)
    vcf_reading=read_Vcf(vcf_file)
    adding_sample_code=add_sample_code(ref_length,vcf_reading[0])
    converting_to_structure=convert_to_structure(aim_pos_reading[1],vcf_reading[1],high_cov_pos,adding_sample_code)
    final_output=out_ref_GT(ref_GT_file,converting_to_structure[0],converting_to_structure[1],directory_path)
    input_file = "{}/{}_Ancestry/1KGP_SuperPop_with_SampleGT.txt".format(sample_path,sample_name)
    structure_output=run_structure(input_file,main_params,extra_params)
    inputfile="{}/{}_Ancestry/OutfileK5withSuperPop_new_f".format(sample_path,sample_name)
    PopNumber=5
    ClusterNo=5
    SampleNo=2504
    structure_result_interpretation(inputfile,PopNumber,ClusterNo,PopCodeFILE,SampleNo,high_cov_pos)

if __name__ == "__main__":
    main()
