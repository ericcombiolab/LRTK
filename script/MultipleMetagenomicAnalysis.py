import os, sys, getopt, json
import subprocess

def main(argv):
    infile=''
    analysis_dir=''
    report_dir=''
    database_dir=''

    try:
        opts, args = getopt.getopt(argv,'-h:-i:-a:-d:-o:', ["help", "infile=", "analysis=", "database=","report="])
    except getopt.GetoptError:
        print("MultipleMetagenomicAnalysis.py -i <infile> -a <analysis dir> -o <report dir>")
        sys.exit(2)

    for opt, arg in opts:
        if opt == "-h":
            print("MultipleMetagenomicAnalysis.py -i <infile> -a <analysis dir> -o <report dir>")
            sys.exit()
        elif opt in ("-a", "--analysis"):
            analysis_dir=arg
        elif opt in ("-d", "--database"):
            database_dir=arg
        elif opt in ("-i", "--infile"):
            infile=arg
        elif opt in ("-o", "--report"):
            report_dir=arg

    scriptdir = os.path.dirname(os.path.abspath( __file__ )) + "/"

    ###preparations
    single_dict={}
    pairwise_dict={}
    group_dict={}
    abd_dict={}
    snv_dict={}
    bam_dict={}
    ref_dict={}

    merge_abd_file = report_dir + "/ABD/Abundance.merge.xls" 
    merge_abd_top_file = report_dir + "/ABD/Abundance.merge.top10.xls"
    merge_abd_top_figure = report_dir + "/ABD/Abundance.merge.top10.jpeg"    
    merge_abd_stat_file = report_dir + "/ABD/Abundance.stastistics.xls"
    merge_pca_file = report_dir + "/ABD/Meta.PCA.xls"
    merge_pca_figure = report_dir + "/ABD/Meta.PCA.jpeg"

    merge_snv_stat_file = report_dir + "/SNV/SNV.statistics.xls"
    merge_snv_top_file  = report_dir + "/SNV/SNV.top10.xls"
    merge_snv_top_figure = report_dir + "/SNV/SNV.top10.jpeg"

    os.makedirs(report_dir + "/ABD/")
    os.makedirs(report_dir + "/SNV/")
    os.makedirs(report_dir + "/img/")

    SI=open(infile,"r")
    line = SI.readline()
    while(line):
        if line.startswith("#"):
            pass
        else:
            arr=line.strip().split("=")
            if(arr[0].startswith("S")):
                single_dict[arr[0]]=arr[1]
                abd_dict[arr[0]]=analysis_dir + "/" + arr[1] + "/BAM/" + "/" + "UHGG.coverage_average.iter.xls"
                snv_dict[arr[0]]=analysis_dir + "/" + arr[1] + "/SNV/" + "/" + arr[1] + ".smallvariants.vcf"
                bam_dict[arr[0]]=analysis_dir + "/" + arr[1] + "/BAM/" + "/" + "/align.round2.bam"
                ref_dict[arr[0]]=analysis_dir + "/" + arr[1] + "/BAM/" + "/" + "/UHGG.reference.subset.fa"
            elif(arr[0].startswith("P")):
                pairwise_dict[arr[0]]=arr[1]
            elif(arr[0].startswith("G")):
                group_dict[arr[0]]=arr[1]

        line = SI.readline()

    SI.close()

    ###abundance
    genomes={}
    info={}
    statistics={}

    for barcode in single_dict.keys():
        abdfile=abd_dict[barcode]
        tmp=[]
        with open(abdfile,"r") as f:
            for line in f:
                a=line.strip().split()
                genomes[a[0]] = 1
                info.setdefault(barcode, {})[a[0]] = eval(a[1])
                tmp.append(eval(a[1]))
            
        values = sorted(tmp, reverse=False)
        lq = int(0.25 * (len(values) + 1)) - 1
        hq = int(0.75 * (len(values) + 1)) - 1
        mq = int(0.5 * (len(values) + 1)) - 1

        statistics[barcode] = {
            "min": values[0],
            "LowQuantile": values[lq],
            "median": values[mq],
            "HighQuantile": values[hq],
            "max": values[-1]
        }

    with open(merge_abd_file, "w") as f:
        f.write("ID")
        for barcode in sorted(single_dict.keys()):
            f.write("\t" + barcode)
        
        f.write("\n")

        for g in sorted(genomes.keys()):
            f.write(g)
            for barcode in sorted(single_dict.keys()):
                value = info.get(barcode, {}).get(g, 0)
                f.write("\t" + str(value))

            f.write("\n")

    with open(merge_abd_stat_file, "w") as f:
        f.write("ID\tMinimum value\tLower Quartile\tMedian\tUpper Quartile\tMaximum value\n")
        for barcode in statistics.keys():
            stats = statistics[barcode]
            f.write("\t".join([barcode, str(stats["min"]), str(stats["LowQuantile"]), str(stats["median"]), str(stats["HighQuantile"]), str(stats["max"])]))
            f.write("\n")

    with open(merge_pca_file, "w") as f:
        f.write("ID\tGroup\n")
        for barcode,values in group_dict.items():
            tmp=values.strip().split(",")
            for i in range(len(tmp)):
                f.write("\t".join([str(tmp[i]), barcode]))
                f.write("\n")

    ###abd figure
    abd_hash = {}
    abd_info = {}
    abd_record = {}
    abd_values = []
    abd_header = []
    
    with open(merge_abd_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("ID"):
                abd_header = line.split("\t")
                continue
            a = line.split("\t")
            abd_record[a[0]]=1
            sum_val = 0
            count = 0
            mean = 0
            
            for i in range(1, len(a)):
                sum_val += float(a[i])
                count += 1
                abd_info.setdefault(abd_header[i], {})[a[0]] = float(a[i])
            
            mean_value = sum_val / count
            abd_values.append(mean_value)
            
            if mean_value not in abd_hash:
                abd_hash[mean_value] = a[0]
            else:
                abd_hash[mean_value] = abd_hash[mean_value] + "," + a[0]
        
    sorted_values = sorted(abd_values, reverse=True)
    
    with open(merge_abd_top_file, 'w') as fw:
        fw.write("ID\tabundance\tspecies\n")
        for i in range(10):
            ele = abd_hash[sorted_values[i]]
            if ele in abd_record:
                for j in range(1, len(abd_header)):
                    fw.write(abd_header[j] + "\t" + str(abd_info[abd_header[j]][ele]) + "\t" + ele + "\n")
            else:
                tmp = ele.split(",")
                for m in range(len(tmp)):
                    for j in range(1, len(abd_header)):
                        fw.write(abd_header[j] + "\t" + str(abd_info[abd_header[j]][tmp[m]]) + "\t" + tmp[m] + "\n")

    ###
    cmd="Rscript " + scriptdir + "/PlotAbundance.R "+ merge_abd_top_file + " " + merge_abd_top_figure
    subprocess.call(cmd, shell=True)
    print(cmd)
    cmd="cp "+ " " + merge_abd_top_figure + " " + report_dir + "/img/abundance.jpeg"
    subprocess.call(cmd, shell=True)

    ###
    cmd="Rscript " + scriptdir + "/PlotPCA.R " + merge_abd_file + " " + merge_pca_file + " " + merge_pca_figure
    subprocess.call(cmd, shell=True)
    print(cmd)
    cmd="cp "+ " " + merge_pca_figure + " " + report_dir + "/img/pca.jpeg"
    subprocess.call(cmd, shell=True)

    ###SNV profile
    ###database
    ref_file = database_dir + "/UHGG/UHGG.reference.txt"
    genome_length = {}
    with open(ref_file, 'r') as file:
        for line in file:
            line = line.strip()
            fields = line.split('\t')
            b = fields[0].split('_')
            genome = b[0] + '_' + b[1]
            if genome in genome_length:
                genome_length[genome] += int(fields[1])
            else:
                genome_length[genome] = int(fields[1])

    ###SNV table
    snv_statistics = {}
    ###SNV figure
    snvcount = {}
    for barcode in single_dict.keys():
        snvfile=snv_dict[barcode]
        total_snv = 0
        homo_snv  = 0
        hete_snv  = 0

        with open(snvfile, 'r') as file:
            for line in file:
                line = line.strip()
                if line.startswith('#'):
                    continue
                
                if "AF1=1" in line:
                    homo_snv += 1
                else:
                    hete_snv += 1

                total_snv += 1

                fields = line.split('\t')
                b = fields[0].split('_')
                genome = b[0] + '_' + b[1]
                if genome in snvcount:
                    if barcode in snvcount[genome]:
                        snvcount[genome][barcode] = snvcount[genome][barcode] + 1
                    else:
                        snvcount.setdefault(genome, {})[barcode] = 1
                else:
                    snvcount.setdefault(genome, {})[barcode] = 1

            snv_statistics[barcode] = {
                "homo": homo_snv / total_snv,
                "hete": hete_snv / total_snv,
                "snv": total_snv
            }

    ###SNV table
    with open(merge_snv_stat_file, 'w') as file:
        file.write("ID\tSNV number\tHomozygous rate\tHeterozygous rate\n")
        for barcode in snv_statistics.keys():
            snv  = snv_statistics[barcode]["snv"]
            homo = snv_statistics[barcode]["homo"]
            hete = snv_statistics[barcode]["hete"]
            file.write(f"{barcode}\t{snv}\t{homo}\t{hete}\n")

    ###SNV figure
    with open(merge_snv_top_file, 'w') as fw:
        fw.write("ID\tSNVnumber\tspecies\n")
        for i in range(10):
            genome=abd_hash[sorted_values[i]]
            if genome in abd_record:
                for j in range(1, len(abd_header)):
                    if(genome in snvcount) and (abd_header[j] in snvcount[genome]):
                        snvrate=snvcount[genome][abd_header[j]]/genome_length[genome]
                        fw.write(abd_header[j] + "\t" + str(snvrate) + "\t" + genome + "\n")
                    else:
                        snvrate=0
                        fw.write(abd_header[j] + "\t" + str(snvrate) + "\t" + genome + "\n")
            else:
                tmp = genome.split(",")
                for m in range(len(tmp)):
                    for j in range(1, len(abd_header)):
                        if(genome in snvcount) and (abd_header[j] in snvcount[genome]):
                            snvrate=snvcount[genome][abd_header[j]]/genome_length[genome]
                            fw.write(abd_header[j] + "\t" + str(snvrate) + "\t" + tmp[m] + "\n")
                        else:
                            snvrate=0
                            fw.write(abd_header[j] + "\t" + str(snvrate) + "\t" + genome + "\n")

    ###
    cmd="Rscript " + scriptdir + "/PlotSNVnumber.R "+ merge_snv_top_file + " " + merge_snv_top_figure
    subprocess.call(cmd, shell=True)
    print(cmd)
    cmd="cp "+ " " + merge_snv_top_figure + " " + report_dir + "/img/snv.jpeg"
    subprocess.call(cmd, shell=True)


    ###SNV allelic imbalance
    ###database
    
    tmp_vcf_list = report_dir + "/SNV/vcf.xls"
    tmp_profile_list = report_dir + "/SNV/profile.xls"
    
    with open(tmp_vcf_list, 'w') as fw:
        for barcode in single_dict.keys():
            snvfile=snv_dict[barcode]
            fw.write(snvfile + "\n")
    
    for i in range(1):
        genome=abd_hash[sorted_values[i]]
        genomesnvdir=report_dir + "/SNV/" + genome
        os.makedirs(genomesnvdir)
        
        if genome in abd_record:
            cmd="python " + scriptdir + "/VCF2SITES.py -i " + tmp_vcf_list + " -g " + genome + " -o " + genomesnvdir + "/SNV." + genome + ".bed"
            subprocess.call(cmd, shell=True)
            print(cmd)
            
            with open(tmp_profile_list, 'w') as fw:
                for barcode in single_dict.keys():
                    cmd="bam-readcount -l " + genomesnvdir + "/SNV." + genome + ".bed" + " -f " + ref_dict[barcode] + " " + bam_dict[barcode] + " -q 20 -b 13 1> " + genomesnvdir + "/" + genome + "." + barcode + ".basecount.xls 2>/dev/null"
                    subprocess.call(cmd, shell=True)
                    print(cmd)
                    
                    cmd="python " + scriptdir + "/BASECOUNT2PROFILE.py -i " + genomesnvdir + "/" + genome + "." + barcode + ".basecount.xls" + " -o "+ genomesnvdir + "/" + genome + "." + barcode + ".profile.xls"
                    subprocess.call(cmd, shell=True)
                    print(cmd)
                    
                    fw.write(genomesnvdir + "/" + genome + "." + barcode + ".profile.xls" + "\n")

            cmd="python " + scriptdir +  "/PROFILE2BAF.py -i " + tmp_profile_list + " -o " + genomesnvdir + "/" + genome + ".baf.xls" + " -g " + genome
            subprocess.call(cmd, shell=True)
            print(cmd)

            outfile=genomesnvdir + "/ABD." + genome + ".xls"
            with open(outfile, 'w') as fw:
                fw.write("ID")
                for barcode in single_dict.keys():
                    fw.write("\t" +barcode)

                fw.write("\n")
                
                fw.write(genome)
                for barcode in single_dict.keys():
                    fw.write("\t" + str(abd_info[barcode][genome]))

                fw.write("\n")

            cmd="Rscript " + scriptdir + "/PlotclusterSNV.R " + genomesnvdir + "/" + genome + ".baf.xls " + genomesnvdir + "/ABD." + genome + ".xls " + genomesnvdir + "/" + "NUMBER." + genome + ".jpeg " + genomesnvdir + "/CLUSTER." + genome + ".jpeg"
            subprocess.call(cmd, shell=True)
            print(cmd)
            
            cmd="cp "+ " " + genomesnvdir + "/CLUSTER." + genome + ".jpeg" + " " + report_dir + "/img/clustering.jpeg"
            subprocess.call(cmd, shell=True)

    ###pairwise
    for i in range(1):
        genome=abd_hash[sorted_values[i]]
        for barcode in pairwise_dict.values():
            genomemaidir=report_dir + "/SNV/" + genome + "/" + barcode
            os.makedirs(genomemaidir)
            tmp_vcf_list = genomemaidir + "/vcf.sub.xls"
            tmp_profile_list = genomemaidir + "/profile.sub.xls"
            tmp_genome_length = genomemaidir + "/genome_length.xls"
            
            with open(tmp_vcf_list, 'w') as fw:
                names=barcode.strip().split("_vs_")
                snvfile=snv_dict[names[0]]
                fw.write(snvfile + "\n")
                snvfile=snv_dict[names[1]]
                fw.write(snvfile + "\n")

            print(genome)
            cmd="python "+ scriptdir + "/VCF2SITES.py -i " + tmp_vcf_list + " -g " + genome + " -o " + genomemaidir + "/SNV." + genome + ".bed"
            subprocess.call(cmd, shell=True)
            print(cmd)

            barcode=names[0]
            cmd="bam-readcount -l " + genomemaidir + "/SNV." + genome + ".bed" + " -f " + ref_dict[barcode] + " " + bam_dict[barcode] + " -q 20 -b 13 1> "+ genomemaidir + "/" + genome + "." + barcode + ".basecount.xls 2>/dev/null"
            subprocess.call(cmd, shell=True)
            print(cmd)

            cmd="python "+ scriptdir + "/BASECOUNT2PROFILE.py -i "+ genomemaidir + "/" + genome + "." + barcode + ".basecount.xls" + " -o " + genomemaidir + "/" + genome + "." + barcode + ".profile.xls"
            subprocess.call(cmd, shell=True)
            print(cmd)

            barcode=names[1]
            cmd="bam-readcount -l " + genomemaidir + "/" + "SNV." + genome + ".bed" + " -f " + ref_dict[barcode] + " " + bam_dict[barcode] + " -q 20 -b 13 1> "+ genomemaidir + "/" + genome + "." + barcode + ".basecount.xls 2>/dev/null"
            subprocess.call(cmd, shell=True)
            print(cmd)

            cmd="python "+ scriptdir + "/BASECOUNT2PROFILE.py -i "+ genomemaidir + "/" + genome + "." + barcode + ".basecount.xls" + " -o "+ genomemaidir + "/" + genome + "." + barcode + ".profile.xls"
            subprocess.call(cmd, shell=True)
            print(cmd)

            with open(tmp_profile_list, 'w') as fw:
                fw.write(genomemaidir + "/" + genome + "." + names[0] + ".profile.xls" + "\n")
                fw.write(genomemaidir + "/" + genome + "." + names[1] + ".profile.xls" + "\n")

            cmd="python "+ scriptdir + "/PROFILE2BAF.py -i " + tmp_profile_list + " -o " + genomemaidir + "/" + genome + ".baf.xls" + " -g " + genome
            subprocess.call(cmd, shell=True)
            print(cmd)

            ###draw
            ref_file = database_dir + "/UHGG/UHGG.reference.txt"
            with open(ref_file, 'r') as file:
                with open(tmp_genome_length, 'w') as fw:
                    for line in file:
                        line = line.strip()
                        fields = line.split('\t')
                        b = fields[0].split('_')
                        refgenome = b[0] + '_' + b[1]
                        if refgenome == genome:
                            fw.write(str(fields[1]) + "\n")
            
            cmd="Rscript "+ scriptdir + "/PlotMAI.R " + genomemaidir + "/" + genome + ".baf.xls " + tmp_genome_length  + " " + genomemaidir + "/pairwise." + genome + ".jpeg " + names[0] + " " + names[1] 
            subprocess.call(cmd, shell=True)
            print(cmd)
            cmd="cp "+ " " + genomemaidir + "/pairwise." + genome + ".jpeg " + " " + report_dir + "/img/MAI.jpeg"
            subprocess.call(cmd, shell=True)
        
if __name__ == "__main__":
        main(sys.argv[1:])

