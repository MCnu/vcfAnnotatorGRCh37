'''
VCF annotator written by M. Chas Summner (MCnu)

INPUT:
VCF of human variants aligned to GRCh37 reference genome

OUTPUT:
TSV of annotated variants containing OUTCOLS as solumns
see README or OUTCOLS constant below for descriptions
'''
'''
IMPORTS
'''
import vcf
import csv
import requests
import argparse
from tqdm import tqdm
from collections import OrderedDict

'''
INPUTS
'''

'''
CONSTANTS
'''
## address for GRCh37 VEP API
REQSERVER = "https://grch37.rest.ensembl.org"

## list of strings to be used for output dataframe columns and loop dict keys
OUTCOLS = [
    "CHROM", ## Chromosome
    "GRCH37_POS", ## Position on chromosome from GRCh37 gneome assembly
    "REF", ## Reference Allele
    "ALT", ## List of Alternative Alleles
    "TOTCOV", ## Total read coverage of position
    "VARREADS", ## List of read count supporting each alt 
    "PERCVARREADS", ## Percentage of total reads supporing each alt
    "VCF_AAF", ## ALT allele frequency in vcf samples (from pyvcf)
    "VCF_FILTER", ## List of potential fitlers in vcf
    "VCF_VARTYPE", ## Variant type (from pyvcf)
    "VCF_VARSUBTYPE", ## Variant subtype (from pyvcf)
    "HGVS_GENEID", ## ID of gene containing variant, from hgvs
    "HGVS_GENESYMBOL", ## Name or symbol of gene contianing variant, from hgvs
    "HGVS_VARTYPE", ## Type of variant query used for hgvs notation
    "HGVS_EFFECT", ## Most Severe Consequence reported by hgvs
    "HGVS_AAF", ## Alternative allele frequency, from hgvs 
    "HGVS_COLOCVAR", ## rsID of other known variants at same position
]

## dicitonary to translate chromsome to accession number
HGVS_CHROMDICT_GRCH37 = {
    "1": "NC_000001.10",
    "2": "NC_000002.11",
    "3": "NC_000003.11",
    "4": "NC_000004.11",
    "5": "NC_000005.9",
    "6": "NC_000006.11",
    "7": "NC_000007.13",
    "8": "NC_000008.10",
    "9": "NC_000009.11",
    "10": "NC_000010.10",
    "11": "NC_000011.9",
    "12": "NC_000012.11",
    "13": "NC_000013.10",
    "14": "NC_000014.8",
    "15": "NC_000015.9",
    "16": "NC_000016.9",
    "17": "NC_000017.10",
    "18": "NC_000018.9",
    "19": "NC_000019.9",
    "20": "NC_000020.10",
    "21": "NC_000021.8",
    "22": "NC_000022.10",
    "23": "NC_000023.10",
    "24": "NC_000024.9",
    "x": "NC_000023.10",
    "y": "NC_000024.9",
    "X": "NC_000023.10",
    "Y": "NC_000024.9",
    "M": "NC_012920.1",
    "m": "NC_012920.1",
    "MT": "NC_012920.1",
    "mt": "NC_012920.1",
}
'''
FUNCTIONS
'''
## input record from vcf.Reader
## returns dicitonary with formatted annotations
def recordAnnotator(record):
    ## initialize ordered dict for record and query info
    ## values will be extracted into a single row for tsv
    loop_dict = OrderedDict.fromkeys(OUTCOLS)
    ## fill dict with parsed vcf info
    loop_dict["CHROM"] = record.CHROM
    loop_dict["GRCH37_POS"] = record.POS
    loop_dict["REF"] = record.REF
    loop_dict["ALT"] = record.ALT
    loop_dict["TOTCOV"] = record.INFO["TC"]
    loop_dict["VARREADS"] = record.INFO["TR"]
    ## calculate the percent of reads per variant
    loop_dict["PERCVARREADS"] = [
        round((tr / record.INFO["TC"]) * 100, 3) for tr in record.INFO["TR"]
    ]
    loop_dict["VCF_AAF"] = record.aaf
    loop_dict["VCF_FILTER"] = record.FILTER
    loop_dict["VCF_VARTYPE"] = record.var_type
    loop_dict["VCF_VARSUBTYPE"] = record.var_subtype
    ## initialize hgvs list with space for each alt allele lookup
    loop_dict["HGVS_GENEID"] = ["NA"] * len(record.ALT)
    loop_dict["HGVS_GENESYMBOL"] = ["NA"] * len(record.ALT)
    loop_dict["HGVS_VARTYPE"] = ["NA"] * len(record.ALT)
    loop_dict["HGVS_EFFECT"] = ["NA"] * len(record.ALT)
    loop_dict["HGVS_AAF"] = ["NA"] * len(record.ALT)
    loop_dict["HGVS_COLOCVAR"] = ["NA"] * len(record.ALT)
    ## use vcf info as extension for VEP hgvs Query
    for i, currALT in enumerate(record.ALT):
        ## convert current ALT to string for hgvs notation
        currALT = str(currALT)
        ## translate CHROM if single number or shortened string is used
        if record.CHROM in HGVS_CHROMDICT_GRCH37.keys():
            currCHROM = HGVS_CHROMDICT_GRCH37[record.CHROM]
        ## or just pass the vcf CHROM if not in dict
        else:
            currCHROM = record.CHROM
        ## determine type of mutation and convert to hgvs notation
        if loop_dict["VCF_VARTYPE"] == "snp":
            REQEXT = (
                f"/vep/human/hgvs/{currCHROM}:g.{record.POS}{record.REF}>{currALT}?"
            )
            loop_dict["HGVS_VARTYPE"][i] = "substitution"
        elif loop_dict["VCF_VARTYPE"] == "indel":
            if len(record.REF) == len(currALT):
                REQEXT = (
                    f"/vep/human/hgvs/{currCHROM}:g.{record.POS}_{record.POS + len(record.REF) - 1}delins{currALT}?"
                )
                loop_dict["HGVS_VARTYPE"][i] = "deletion-insertion"
            elif len(record.REF) > len(currALT):
                REQEXT = (
                    f"/vep/human/hgvs/{currCHROM}:g.{record.POS+len(currALT)-1}_{record.POS + len(record.REF) - 1}del?"
                )
                loop_dict["HGVS_VARTYPE"][i] = "deletion"
            elif len(record.REF) < len(currALT):
                REQEXT = (
                    f"/vep/human/hgvs/{currCHROM}:g.{record.POS+len(record.REF)-1}_{record.POS+len(record.REF)}ins{currALT[len(record.REF):]}?"
                )
                loop_dict["HGVS_VARTYPE"][i] = "insertion"
            else:
                raise Exception("Unsupported indel subtype:", loop_dict)
        else:
            raise Exception("Unsupported Variant Type:", loop_dict)
        req = requests.get(
            REQSERVER + REQEXT, headers={"Content-Type": "application/json"}
        )
        ## for ok queries, assign values to dict
        if req.ok:
            req = req.json()
            ## check existence before attemping assignment
            if "most_severe_consequence" in req[0]:
                loop_dict["HGVS_EFFECT"] = req[0]["most_severe_consequence"]
            if "colocated_variants" in req[0]:
                if "id" in req[0]["colocated_variants"][0]:
                    loop_dict["HGVS_COLOCVAR"][i] = req[0]["colocated_variants"][0]["id"]
                if "frequencies" in req[0]["colocated_variants"][0]:
                    loop_dict["HGVS_AAF"][i] = req[0]["colocated_variants"][0]["frequencies"]
            if "transcript_consequences" in req[0]:
                if "gene_id" in req[0]["transcript_consequences"][0]:
                    loop_dict["HGVS_GENEID"][i] = req[0]["transcript_consequences"][0]["gene_id"]
                if "gene_symbol" in req[0]["transcript_consequences"][0]:
                    loop_dict["HGVS_GENESYMBOL"][i] = req[0]["transcript_consequences"][0]["gene_symbol"]
            if (
                "intergenic_consequences" not in req[0]
                and "transcript_consequences" not in req[0]
            ):
                print("Query did not contain consequences:", loop_dict)
        else:
            print("Following VCF record could not be queried:", loop_dict)
    return loop_dict

## input paths to target vcf and output tsv
## reads in vcf, runs recordAnnotator() on each record
## writes row for each record
def main(input_vcf, output_tsv):
    ## initialize vcf reader
    try:
        vcf_reader = vcf.Reader(open(input_vcf, "r"))
    except:
        raise ValueError("Provide a valid vcf as input vcf")
    ## initialize output tsv
    with open(output_tsv, "x") as out_file:
        tsv_writer = csv.writer(out_file, delimiter="\t")
        ## write column names
        tsv_writer.writerow(OUTCOLS)
        ## loop through each record to annotate and write to tsv
        for record in tqdm(vcf_reader):
            loop_dict = recordAnnotator(record)
            tsv_writer.writerow(loop_dict.values())

if __name__ == "__main__":
    ## take in two arguments, input vcf and output tsv
    parser = argparse.ArgumentParser(prog="hgvsVEPQuery_vcfAnnotator_grch37")
    parser.add_argument("-i", "--input_vcf")
    parser.add_argument("-o", "--output_tsv")
    args = parser.parse_args()

    main(args.input_vcf, args.output_tsv)
