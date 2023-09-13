import os, sys, getopt, json
import subprocess

def main(argv):
    abdfile=''
    snvfile=''
    HTMLreport=''
    
    try:
        opts, args = getopt.getopt(argv,'-h:-a:-s:-o:', ["help","abundance=", "snv=", "HTMLreport="])
    except getopt.GetoptError:
        print("report.py -a <abundance> -s <snvfile> -o <output HTML file>")
        sys.exit(2)

    for opt, arg in opts:
        print(arg)
        if opt == "-h":
            print("report.py -a <abundance> -s <snvfile> -o <output HTML file>")
            sys.exit()
        elif opt in ("-a", "--abundance"):
            abdfile=arg
        elif opt in ("-s", "--snv"):
            snvfile=arg
        elif opt in ("-o", "--HTMLreport"):
            HTMLreport=arg

    outdir = os.path.dirname(os.path.abspath(HTMLreport)) + "/"
    scriptdir = os.path.dirname(os.path.abspath( __file__ )) + "/"

    ###preparations
    cmd="cp -r " + scriptdir + "HTML/css " + scriptdir + "HTML/js "  + outdir
    subprocess.call(cmd, shell=True)
    cmd="cp -r " + scriptdir + "HTML/img/workflow.png " + outdir + "/img"
    subprocess.call(cmd, shell=True)

    abddict={}
    AB=open(abdfile,"r")
    line = AB.readline()
    while(line):
        if line.startswith("ID"):
            pass
        else:
            arr=line.strip().split("\t")
            barcode=arr[0]
            abddict[barcode]=line
        line = AB.readline()
    AB.close()

    snvdict={}
    SS=open(snvfile,"r")
    line = SS.readline()
    while(line):
        if line.startswith("ID"):
            pass
        else:
            arr=line.strip().split("\t")
            barcode=arr[0]
            snvdict[barcode]=line
        line = SS.readline()
    SS.close()

    ###write report
    HTML=open(HTMLreport,"w")
    
    dec_text="<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\"\n\"http://www.w3.org/TR/html4/loose.dtd\">\n"
    HTML.write(dec_text)
    dec_text="<html>\n<head>\n<meta http-equiv=\"Content-Type\" content=\"text/html; charset=gb2312\">\n<link href=\"css/report.css\" rel=\"stylesheet\" type=\"text/css\">\n<script language=\"javascript\" src=\"js/jquery.js\"></script>\n<script language=\"javascript\" src=\"js/report.js\"></script>\n<title>LRTK report for multiple related samples </title>\n<style type=\"text/css\">\nbody{\n\tposition:relative;\n\toverflow:auto;\n}\n.note {color: #f00;}\n\n</style>\n</head>\n"
    HTML.write(dec_text)
    
    ###overview
    dec_text="<body>\n<div class=seq_cn1>\n<div class=\"top_btn1\"></div>\n"
    HTML.write(dec_text)
    dec_text="<!-- Part I Overview-->\n<div class=\"h0\"><a name=\"overview\"></a>Part I Overview</div>\n<hr class=\"clear-contentunit\" />\n<h1><a name=\"overview.bio1\"></a> Standard Bioinformatics Analysis</h1>\n<p>Linked-Read ToolKit (LRTK), is a unified and versatile toolkit to analyze both metagenome and human genome linked-read sequencing data derived from any of the three major linked-read sequencing platforms. LRTK delivers a suite of utilities to perform data simulation, format conversion, data preprocessing, read cloud assembly, barcode-aware alignment, reconstruction of long fragments, and variant detection and phasing. LRTK automatically produces HTML reports to summarize the quality statistics as part of its pipeline, and generates publication-ready visualizations.</p><br/>\n<br/>\n<div class=\"tc\"><a href=\"./img/workflow.png\" target=\"_blank\"><img style=\"width:720px\" src=\"./img/workflow.png\" /></a></div>\n<br/>\n<hr class=\"clear-contentunit-thin\"/>\n<br/>\n"
    HTML.write(dec_text)

    ###bioinformatics analysis
    dec_text="<!--Part 2 bioinformatics analysis-->\n\n<!--sequencing QC-->\n<div class=\"h0\"><a name=\"result\"></a>Part II Bioinformatics Analysis Results </div>\n<hr class=\"clear-contentunit\" /> \n<h1><a name=\"result.abundance\"></a>2.1 Taxonomic abundance</h1>\n<p>We developed a reference genome-based framework to quantify the abundance of species. In brief, we align decontaminated metagenome reads to the reference genome using a barcode aware alignment pattern. We adopt a tiered alignment strategy to detect species. In the initial alignment, we map reads to the complete reference genome. We then divided the reference genome into 1 k bp windows and calculated the number of mapped reads for each window. The mapped reads are classified into two types: unique mapping reads (U) and multiple mapping reads (M). For each window W, we calculate the number of unique RC(U) and multiple RC(M) reads, respectively. For each multiple mapping read, we calculate the weight ω in terms of the number of aligned locations. To ensure proper estimation of genome abundance, we calculate the abundance as the mean of the n% most closely covered windows. After initial classification of microbial genomes, we extract the corresponding genome sequence and align the reads to the new subset reference genome sequence to refine the alignment. The final relative abundance is calculated again by using the same  processing approach.</p>\n<hr class=\"clear-contentunit-thin\"/>\n"
    HTML.write(dec_text)

    dec_text="<div class=\"break\"></div>\n<!-- Content unit - One column -->\nData Statistics<br/><table>\n<tr><th >ID</th><th>Minimum value</th><th>Lower Quartile</th><th>Median</th><th>Upper Quartile</th><th>Maximum value</th></tr>\n"
    HTML.write(dec_text)
    for key, value in abddict.items():
        arr=value.strip().split("\t")
        dec_text="<tr><td>" + arr[0] + "</td><td>" + arr[1] + "</td><td>" + arr[2] + "</td><td>" + arr[3] + "</td><td>" + arr[4] + "</td><td>" + arr[5] + "</td></tr>"
        HTML.write(dec_text)
        
    dec_text="</table>"
    HTML.write(dec_text)

    ###ABD
    dec_text="<hr class=\"clear-contentunit-thin\" /><br/>Data distribution<br/>\n<div class=\"tc\"><a href=\"./img/abundance.jpeg\" target=\"_blank\"><img style=\"width:720px\" src=\"./img/abundance.jpeg\" /></a></div>\n"
    HTML.write(dec_text)

    ###PCA
    dec_text="<hr class=\"clear-contentunit-thin\" />\n<h1><a name=\"result.abundance\"></a> 2.1.1 Dimensionality Reduction </h1>\n<br/>Data distribution<br/>\n<div class=\"tc\"><a href=\"./img/pca.jpeg\" target=\"_blank\"><img style=\"width:720px\" src=\"./img/pca.jpeg\" /></a></div>\n"
    HTML.write(dec_text)

    ###Network
    #dec_text="<hr class=\"clear-contentunit-thin\" />\n<h1><a name=\"result.abundance\"></a> 2.1.2 Network </h1>\n<br/>Data distribution<br/>\n<div class=\"tc\"><a href=\"./img/network.jpeg\" target=\"_blank\"><img style=\"width:720px\" src=\"./img/network.jpeg\" /></a></div>\n"
    #HTML.write(dec_text)

    ###SNV
    dec_text="<hr class=\"clear-contentunit-thin\"/>\n<div class=\"break\"></div>\n<h1><a name=\"result.SNVs\"></a> 2.2 Single Nucleotide Variations</h1>\nLRTK supports several popular variant detection tools and promotes best practice for linked-read sequencing. For SNV and INDEL calling, LRTK provides FreeBayes, SAMtools and GATK.\n<hr class=\"clear-contentunit-thin\" />\nData Statistics<br/>\n<table>\n<tr><th>ID</th><th>SNV number</th><th>Homozygous rate</th><th>Heterozygous rate</th></tr>\n" 
    HTML.write(dec_text)
    for key, value in abddict.items():
        arr=value.strip().split("\t")
        dec_text="<tr><td>" + arr[0] + "</td><td>" + arr[1] + "</td><td>" + arr[2] + "</td><td>" + arr[3] + "</td></tr>"
        HTML.write(dec_text)

    dec_text="</table>"
    HTML.write(dec_text)

    ###SNV number
    dec_text="<hr class=\"clear-contentunit-thin\"/>\n<br/>Data distribution<br/>\n<div class=\"tc\"><a href=\"./img/snv.jpeg\" target=\"_blank\"><img style=\"width:720px\" src=\"./img/snv.jpeg\" /></a></div>\n"
    HTML.write(dec_text)

    ###clustering
    dec_text="<hr class=\"clear-contentunit-thin\"/>\n<h1><a name=\"result.SNVs\"></a> 2.2.1 Unsupervised Clustering </h1>\n<br/>Data distribution<br/>\n<div class=\"tc\"><a href=\"./img/clustering.jpeg\" target=\"_blank\"><img style=\"width:720px\" src=\"./img/clustering.jpeg\"/></a></div>\n"
    HTML.write(dec_text)

    ##allelic imbalance
    dec_text="<h1><a name=\"result.SNVs\"></a> 2.2.2 Mirrored allelic imbalance</h1>\n<br/>Data distribution<br/>\n<div class=\"tc\"><a href=\"./img/MAI.jpeg\" target=\"_blank\"><img style=\"width:720px\" src=\"./img/MAI.jpeg\" /></a></div>"
    HTML.write(dec_text)

    ###end
    dec_text="</div>\n</body>\n</html>\n"
    HTML.write(dec_text)

if __name__ == "__main__":
        main(sys.argv[1:])

