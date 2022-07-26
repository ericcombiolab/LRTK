import sys, getopt, json

def main(argv):
	FQjson=''
	BAMjson=''
	VARjson=''
	HTMLreport=''

	try:
		opts, args = getopt.getopt(argv,'-h:-f:-b:-m:-o:', ["help","FQjson=", "BAMjson=", "VARjson=", "HTMLreport="])
	except getopt.GetoptError:
		print("report.py -f <input FQ json file> -b <input BAM json file> -m <input variants json file> -o <output HTML file>")
		sys.exit(2)
	
	for opt, arg in opts:
		print(arg)
		if opt == "-h":
			print("report.py -f <input FQ json file> -b <input BAM json file> -m <input variants json file> -o <output HTML file>")
			sys.exit()
		elif opt in ("-f", "--FQjson"):
			FQjson=arg
		elif opt in ("-b", "--BAMjson"):
			BAMjson=arg
		elif opt in ("-m", "--VARjson"):
			VARjson=arg
		elif opt in ("-o", "--HTMLreport"):
			HTMLreport=arg

#	print(FQjson)
	FQ=open(FQjson,"r")
	FQdata=json.load(FQ)
	FQ.close()
	
#	print(FQdata['read1_before_filtering']["quality_curves"])
	
	BAM=open(BAMjson,"r")
	BAMdata=json.load(BAM)
	BAM.close()

	VCF=open(VARjson,"r")
	VCFdata=json.load(VCF)
	VCF.close()
#	print(VCFdata)
	HTML=open(HTMLreport,"w")
	##Header
	dec_text = "<html><head><meta http-equiv=\"content-type\" content=\"text/html;charset=utf-8\" />\n<title>Linked-read report at 2021-04-09      09:32:08 </title><script src='http://opengene.org/plotly-1.2.0.min.js'></script>\n\n<script type=\"text/javascript\">\n\tfunction showOrHide(divname) {\n\t\tdiv = document.getElementById(divname);\n\t\tif(div.style.display == 'none')\n\t\t\tdiv.style.display = 'block';\n\t\telse\n\t\t\tdiv.style.display = 'none';\n\t}\n</script>\n\n";
	HTML.write(dec_text)
	dec_text = "<style type=\"text/css\">\ntd {border:1px solid #dddddd;padding:5px;font-size:12px;}\ntable {border:1px solid #999999;padding:2x;border-collapse:collapse; width:800px}\n.col1 {width:320px; font-weight:bold;}\n.adapter_col {width:500px; font-size:10px;}\nimg {padding:30px;}\n#menu {font-family:Consolas, 'Liberation Mono', Menlo, Courier, monospace;}\n#menu a {color:#0366d6; font-size:18px;font-weight:600;line-height:28px;text-decoration:none;font-family:-apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol'}\na:visited {color: #999999}\n.alignleft {text-align:left;}\n.alignright {text-align:right;}\n.figure {width:800px;height:600px;}\n.header {color:#ffffff;padding:1px;height:20px;background:#000000;}\n.section_title {color:#ffffff;font-size:14px;padding:7px;text-align:left;background:#70b7c8; margin-top:10px;}\n.subsection_title {font-size:12px;padding:5px;margin-top:10px;text-align:left;color:#70b7c8}\n#container {text-align:center;padding:3px 3px 3px 10px;font-family:Arail,'Liberation Mono', Menlo, Courier, monospace;}\n.menu_item {text-align:left;padding-top:5px;font-size:18px;}\n.highlight {text-align:left;padding-top:30px;padding-bottom:30px;font-size:20px;line-height:35px;}\n#helper {text-align:left;border:1px dotted #fafafa;color:#777777;font-size:12px;}\n#footer {text-align:left;padding:15px;color:#ffffff;font-size:10px;background:#396e47;font-family:Arail,'Liberation Mono', Menlo, Courier, monospace;}\n.kmer_table {text-align:center;font-size:8px;padding:2px;}\n.kmer_table td{text-align:center;font-size:8px;padding:0px;color:#ffffff}\n.sub_section_tips {color:#999999;font-size:10px;padding-left:5px;padding-bottom:3px;text-align:left;}\n</style>\n\n"
	HTML.write(dec_text)	
	dec_text = "</head><body><div id='container'><h1 style='text-align:left;'><a href='https://github.com/CicyYeung/LRTK' target='_blank' style='color:#56461f;text-decoration:none;'>Linked-read report</a</h1>\n<div style='font-size:10px;font-weight:normal;text-align:left;color:#666666;padding:5px;'>Created by <a href=' https://github.com/CicyYeung/LRTK' style='color:#1F77B4'>Linked Read ToolKits</a> v1.0, a unified and versatile toolkit for analyzing Linked-Read sequencing data</div>\n";
	HTML.write(dec_text)
	
	###section1
	readlength_before = len(FQdata['read1_before_filtering']["quality_curves"]["A"])
	readlength_after  = len(FQdata['read1_after_filtering']["quality_curves"]["A"])
	dec_text = "<div class='section_div_p1'>\n<div class='section_title' onclick=showOrHide('result')><a name='result'>Quality control on sequencing reads <font color='#56461f' > (click to show/hide) </font></a></div>\n<div id='result'>\n<div id='detection_result'>\n<table class='summary_table' style='width:800px'>\n"
	HTML.write(dec_text)
	dec_text =  ' '.join(["<tr><td class='col1'>Number of total reads before filtering:</td><td class='col2'>", str(FQdata['summary']['before_filtering']['total_reads']),"K </td></tr>\n"])
	HTML.write(dec_text)
	dec_text =  ' '.join(["<tr><td class='col1'>Read length before filtering:</td><td class='col2'>", str(readlength_before),"bp </td></tr>\n"])
	HTML.write(dec_text)
	dec_text =  ' '.join(["<tr><td class='col1'>Q20 before filtering:</td><td class='col2'>", str(round(FQdata['summary']['before_filtering']['q20_rate'] *100,2)),"% </td></tr>\n"])
	HTML.write(dec_text)
	dec_text =  ' '.join(["<tr><td class='col1'>Q30 before filtering:</td><td class='col2'>", str(round(FQdata['summary']['before_filtering']['q30_rate'] * 100, 2)),"% </td></tr>\n"])
	HTML.write(dec_text)
	dec_text =  ' '.join(["<tr><td class='col1'>GC content before filtering:</td><td class='col2'>", str(round(FQdata['summary']['before_filtering']['gc_content'],2)),"% </td></tr>\n"])
	HTML.write(dec_text)
	dec_text =  ' '.join(["<tr><td class='col1'>Number of total reads after filtering:</td><td class='col2'>", str(FQdata['summary']['after_filtering']['total_reads']),"K </td></tr>\n"])
	HTML.write(dec_text)
	dec_text =  ' '.join(["<tr><td class='col1'>Read length after filtering:</td><td class='col2'>", str(readlength_after),"bp </td></tr>\n"])
	HTML.write(dec_text)
	dec_text =  ' '.join(["<tr><td class='col1'>Q20 after filtering:</td><td class='col2'>", str(round(FQdata['summary']['after_filtering']['q20_rate'] *100,2)),"% </td></tr>\n"])
	HTML.write(dec_text)
	dec_text =  ' '.join(["<tr><td class='col1'>Q30 after filtering:</td><td class='col2'>", str(round(FQdata['summary']['after_filtering']['q30_rate'] * 100, 2)),"% </td></tr>\n"])
	HTML.write(dec_text)
	dec_text =  ' '.join(["<tr><td class='col1'>GC content after filtering:</td><td class='col2'>", str(round(FQdata['summary']['after_filtering']['gc_content'],2)),"% </td></tr>\n"])
	HTML.write(dec_text)
	dec_text="</table>\n</div>\n</div>\n</div>\n\n"
	HTML.write(dec_text)

	###section2
	readmat = [i+1 for i in range(readlength_after)]
	dec_text = "<div class='section_div_p2'>\n<div class='section_title' onclick=showOrHide('pbsq')><a name='summary'> Per base sequencing quality scores<font color='#56461f' > (click to show/hide) </font></a></div>\n<div id='pbsq' style='display:none'>\n<div id='pbsq_figure'>\n<div class='figure' id='plot_pbsq1' style='height:400px;width:800px'></div>\n<script type=\"text/javascript\">\n"
	HTML.write(dec_text)
	dec_text = "".join(["\t","var data=[{x:",str(readmat),",y:",str(FQdata['read1_after_filtering']["quality_curves"]["A"]),",name: 'A',mode:'lines',line:{color:'rgba(128,128,0,1.0)', width:1}},\n"])
	HTML.write(dec_text)
	dec_text = "".join(["\t","{x:",str(readmat),",y:",str(FQdata['read1_after_filtering']["quality_curves"]["T"]),",name: 'T',mode:'lines',line:{color:'rgba(128,0,128,1.0)', width:1}},\n"])
	HTML.write(dec_text)
	dec_text = "".join(["\t","{x:",str(readmat),",y:",str(FQdata['read1_after_filtering']["quality_curves"]["C"]),",name: 'C',mode:'lines',line:{color:'rgba(0,255,0,1.0)', width:1}},\n"])
	HTML.write(dec_text)
	dec_text = "".join(["\t","{x:",str(readmat),",y:",str(FQdata['read1_after_filtering']["quality_curves"]["G"]),",name: 'G',mode:'lines',line:{color:'rgba(0,0,255,1.0)', width:1}},\n"])
	HTML.write(dec_text)
	dec_text = "".join(["\t","{x:",str(readmat),",y:",str(FQdata['read1_after_filtering']["quality_curves"]["mean"]),",name: 'mean',mode:'lines',line:{color:'rgba(20,20,20,1.0)', width:1}},\n"])
	HTML.write(dec_text)
	dec_text = "\t];\n\tvar layout={title:'Base quality for read1', xaxis:{title:'position'}, yaxis:{title:'quality'}};\n\tPlotly.newPlot('plot_pbsq1', data, layout);\n</script>\n\n"
	HTML.write(dec_text)

	read2length_after = len(FQdata['read2_after_filtering']["quality_curves"]["A"])
	readmat = [i+1 for i in range(read2length_after)]
	dec_text = "<div class='figure' id='plot_pbsq2' style='height:600px;width:800px'></div>\n<script type=\"text/javascript\">\n"
	HTML.write(dec_text)
	dec_text = "".join(["\t","var data=[{x:",str(readmat),",y:",str(FQdata['read2_after_filtering']["quality_curves"]["A"]),",name: 'A',mode:'lines',line:{color:'rgba(128,128,0,1.0)', width:1}},\n"])
	HTML.write(dec_text)
	dec_text = "".join(["\t","{x:",str(readmat),",y:",str(FQdata['read2_after_filtering']["quality_curves"]["T"]),",name: 'T',mode:'lines',line:{color:'rgba(128,0,128,1.0)', width:1}},\n"])
	HTML.write(dec_text)
	dec_text = "".join(["\t","{x:",str(readmat),",y:",str(FQdata['read2_after_filtering']["quality_curves"]["C"]),",name: 'C',mode:'lines',line:{color:'rgba(0,255,0,1.0)', width:1}},\n"])
	HTML.write(dec_text)
	dec_text = "".join(["\t","{x:",str(readmat),",y:",str(FQdata['read2_after_filtering']["quality_curves"]["G"]),",name: 'G',mode:'lines',line:{color:'rgba(0,0,255,1.0)', width:1}},\n"])
	HTML.write(dec_text)
	dec_text = "".join(["\t","{x:",str(readmat),",y:",str(FQdata['read2_after_filtering']["quality_curves"]["mean"]),",name: 'mean',mode:'lines',line:{color:'rgba(20,20,20,1.0)', width:1}},\n"])
	HTML.write(dec_text)
	dec_text = "\t];\n\tvar layout={title:'Base quality for read2', xaxis:{title:'position'}, yaxis:{title:'quality'}};\n\tPlotly.newPlot('plot_pbsq2', data, layout);\n"
	HTML.write(dec_text)
	dec_text = "</script>\n</div>\n</div>\n</div>\n\n"
	HTML.write(dec_text)
	
	###section3
	dec_text = "<div class='section_div_p3'>\n<div class='section_title' onclick=showOrHide('pbgc')><a name='summary'> Per base GC-content <font color='#56461f' > (click to show/hide) </font></a></div>\n<div id='pbgc' style='display:none'>\n<div id='pbgc_figure'>\n<div class='figure' id='plot_pbgc1' style='height:400px;width:800px'></div>\n<script type=\"text/javascript\">\n"
	HTML.write(dec_text)
	dec_text = "".join(["\t","var data=[{x:",str(readmat),",y:",str(FQdata['read1_after_filtering']["content_curves"]["A"]),",name: 'A',mode:'lines',line:{color:'rgba(128,128,0,1.0)', width:1}},\n"])
	HTML.write(dec_text)
	dec_text = "".join(["\t","{x:",str(readmat),",y:",str(FQdata['read1_after_filtering']["content_curves"]["T"]),",name: 'T',mode:'lines',line:{color:'rgba(128,0,128,1.0)', width:1}},\n"])
	HTML.write(dec_text)
	dec_text = "".join(["\t","{x:",str(readmat),",y:",str(FQdata['read1_after_filtering']["content_curves"]["C"]),",name: 'C',mode:'lines',line:{color:'rgba(0,255,0,1.0)', width:1}},\n"])
	HTML.write(dec_text)
	dec_text = "".join(["\t","{x:",str(readmat),",y:",str(FQdata['read1_after_filtering']["content_curves"]["G"]),",name: 'G',mode:'lines',line:{color:'rgba(0,0,255,1.0)', width:1}},\n"])
	HTML.write(dec_text)
	dec_text = "".join(["\t","{x:",str(readmat),",y:",str(FQdata['read1_after_filtering']["content_curves"]["N"]),",name: 'N',mode:'lines',line:{color:'rgba(255, 0, 0, 1.0)', width:1}},\n"])
	HTML.write(dec_text)
	dec_text = "".join(["\t","{x:",str(readmat),",y:",str(FQdata['read1_after_filtering']["content_curves"]["GC"]),",name: 'GC',mode:'lines',line:{color:'rgba(20,20,20,1.0)', width:2}},\n"])
	HTML.write(dec_text)
	dec_text = "\t];\n\tvar layout={title:'GC-content for read1', xaxis:{title:'position'}, yaxis:{title:'base content ratios'}};\n\tPlotly.newPlot('plot_pbgc1', data, layout);\n</script>\n\n"
	HTML.write(dec_text)

	read2length_after = len(FQdata['read2_after_filtering']["content_curves"]["A"])
	readmat = [i+1 for i in range(read2length_after)]
	dec_text = "<div class='figure' id='plot_pbgc2' style='height:400px;width:800px'></div>\n<script type=\"text/javascript\">\n"
	HTML.write(dec_text)
	dec_text = "".join(["\t","var data=[{x:",str(readmat),",y:",str(FQdata['read2_after_filtering']["content_curves"]["A"]),",name: 'A',mode:'lines',line:{color:'rgba(128,128,0,1.0)', width:1}},\n"])
	HTML.write(dec_text)
	dec_text = "".join(["\t","{x:",str(readmat),",y:",str(FQdata['read2_after_filtering']["content_curves"]["T"]),",name: 'T',mode:'lines',line:{color:'rgba(128,0,128,1.0)', width:1}},\n"])
	HTML.write(dec_text)
	dec_text = "".join(["\t","{x:",str(readmat),",y:",str(FQdata['read2_after_filtering']["content_curves"]["C"]),",name: 'C',mode:'lines',line:{color:'rgba(0,255,0,1.0)', width:1}},\n"])
	HTML.write(dec_text)
	dec_text = "".join(["\t","{x:",str(readmat),",y:",str(FQdata['read2_after_filtering']["content_curves"]["G"]),",name: 'G',mode:'lines',line:{color:'rgba(0,0,255,1.0)', width:1}},\n"])
	HTML.write(dec_text)
	dec_text = "".join(["\t","{x:",str(readmat),",y:",str(FQdata['read2_after_filtering']["content_curves"]["N"]),",name: 'N',mode:'lines',line:{color:'rgba(255, 0, 0, 1.0)', width:1}},\n"])
	HTML.write(dec_text)
	dec_text = "".join(["\t","{x:",str(readmat),",y:",str(FQdata['read2_after_filtering']["content_curves"]["GC"]),",name: 'GC',mode:'lines',line:{color:'rgba(20,20,20,1.0)', width:2}},\n"])
	HTML.write(dec_text)
	dec_text = "\t];\n\tvar layout={title:'GC-content for read2', xaxis:{title:'position'}, yaxis:{title:'base content ratios'}};\n\tPlotly.newPlot('plot_pbgc2', data, layout);\n"
	HTML.write(dec_text)
	dec_text = "</script>\n</div>\n</div>\n</div>\n\n"
	HTML.write(dec_text)

	###section4
	dec_text = "<div class='section_div_p4'>\n<div class='section_title' onclick=showOrHide('summary')><a name='summary'>Barcode-aware read alignment <font color='#56461f' > (click to show/hide) </font></a></div>\n<div id='summary'>\n<div class='subsection_title' onclick=showOrHide('general')>Summary Table</div>\n<div id='general'>\n<table class='summary_table'>\n"	
	HTML.write(dec_text)
	dec_text="<tr><td class='col1'> Number of barcodes before error correction:</td><td class='col2'>" + "NA" + "</td></tr>\n"
	HTML.write(dec_text)
	dec_text="<tr><td class='col1'> Number of barcodes after error correction:</td><td class='col2'>" + str(BAMdata['TotalBarcodes']) + "</td></tr>\n"
	HTML.write(dec_text)
	dec_text="<tr><td class='col1'> Percentage of validated barcodes:</td><td class='col2'>" + "NA" + "</td></tr>\n"
	HTML.write(dec_text)
	dec_text="".join(["<tr><td class='col1'> Duplication rate:</td><td class='col2'>", str(round(float(BAMdata['DuplicationRate']) * 100,2))," %</td></tr>\n"])
	HTML.write(dec_text)
	dec_text="<tr><td class='col1'> Average phasical coverage of long fragments (C<sub>F</sub>):</td><td class='col2'>"+str(round(float(BAMdata['C_F']),2)) +" X </td></tr>\n"
	HTML.write(dec_text)
	dec_text="<tr><td class='col1'> Average short-read coverage per fragment (C<sub>R</sub>):</td><td class='col2'>" + str(round(float(BAMdata['C_R']),2)) + " X </td></tr>\n"
	HTML.write(dec_text)
	dec_text="<tr><td class='col1'> Mean depth:</td><td class='col2'>" + str(round(float(BAMdata['C']),2)) +" X </td></tr>\n"
	HTML.write(dec_text)
	dec_text="<tr><td class='col1'> Average fragment length (u<sub>FL</sub>):</td><td class='col2'>" + str(round(float(BAMdata['U_FL']),2)) + " bp </td></tr>\n"
	HTML.write(dec_text)
	dec_text="<tr><td class='col1'> Weighted average fragment length (Wu<sub>FL</sub>):</td><td class='col2'>" + str(round(int(BAMdata['WU_FL']),2)) +" bp</td></tr>\n"
	HTML.write(dec_text)
	dec_text="<tr><td class='col1'> Number of fragments per barcode (N<sub>F/P</sub>):</td><td class='col2'>" + str(round(float(BAMdata['NF_P']),2)) + "</td></tr>\n"
	HTML.write(dec_text)
	dec_text="<tr><td class='col1'> Median insert size:</td><td class='col2'>" + str(round(float(BAMdata['InsertSize']),2)) + " bp</td></tr>\n"
	HTML.write(dec_text)
	dec_text="</table>\n</div>\n</div>\n</div>\n\n"
	HTML.write(dec_text)

	###section5
	nfp_array=[i for i in range(50)]
	dec_text="<div class='section_div_p5'>\n<div class='section_title' onclick=showOrHide('nfp')><a name='summary'> Number of fragment per barcode <font color='#56461f' > (click to show/hide) </font></a></div>\n<div id='nfp' style='display:none'>\n<div id='nfp_figure'>\n<div class='figure' id='nfp' style='height:400px;width:800px;'></div>\n</div>\n<script type=\"text/javascript\">\n"
	HTML.write(dec_text)
	dec_text="".join(["\tvar data=[{\n\tx:",str(nfp_array),",\n\ty:",str(BAMdata['density_nfp']) ,",\n\tname: 'Percentages (%)  ',type:'bar',line:{color:'rgba(128,0,128,1.0)', width:1}\n\t}];","\n\tvar layout={title:'number of fragments per barcode distribution', xaxis:{title:'number of fragment'}, yaxis:{title:'percentages (%)'}};\n\tPlotly.newPlot('nfp', data, layout);\n"])
	HTML.write(dec_text)
	dec_text="</script>\n</div>\n</div>\n\n"
	HTML.write(dec_text)
	
	###section6
	fl_array=[i for i in range(100)]
	dec_text="<div class='section_div_p6'>\n<div class='section_title' onclick=showOrHide('fragment_length')><a name='summary'> Fragment length <font color='#56461f' > (click to show/hide) </font></a></div>\n<div id='fragment_length' style='display:none'>\n<div id='fragment_length'>\n<div class='figure' id='fragment_length_figure' style='height:400px;width:800px;'></div>\n</div>\n<script type=\"text/javascript\">\n"
	HTML.write(dec_text)
	dec_text="".join(["\tvar data=[{\n\tx:",str(fl_array),",\n\ty:",str(BAMdata['density_fl']),",\n\tname: 'Percent (%)  ',type:'bar',line:{color:'rgba(128,0,128,1.0)', width:1}\n\t}];\n\tvar layout={title:'Fragment length distribution', xaxis:{title:'fragment length (kb)'}, yaxis:{title:' percentages (%)'}};\n\tPlotly.newPlot('fragment_length_figure', data, layout);\n"])
	HTML.write(dec_text)
	dec_text="</script>\n</div>\n</div>\n\n"
	HTML.write(dec_text)

	###section7
	cr_array=[i for i in range(100)]
	dec_text="<div class='section_div_p7'>\n<div class='section_title' onclick=showOrHide('coverage_rate')><a name='summary'> Coverage of reads per fragment <font color='#56461f' > (click to show/hide) </font></a></div>\n<div id='coverage_rate' style='display:none'>\n<div id='coverage_rate_figure'>\n<div class='figure' id='plot_coverage_rate' style='height:400px;width:800px;'></div>\n</div>\n<script type=\"text/javascript\">\n"
	HTML.write(dec_text)
	dec_text="".join(["\tvar data=[{\n\tx:",str(cr_array),",\n\ty:",str(BAMdata['density_cr']),",\n\tname: 'percentages (%)  ',type:'bar',line:{color:'rgba(128,0,128,1.0)', width:1}\n\t}];\n\tvar layout={title:'Coverage of read per fragment', xaxis:{title:'coverage rate (X)'}, yaxis:{title:' percentages (%)'}};\n\tPlotly.newPlot('plot_coverage_rate', data, layout);\n"])
	HTML.write(dec_text)
	dec_text="</script>\n</div>\n</div>\n\n"
	HTML.write(dec_text)

	###section10
	var_array=[i+50 for i in range(450)]
	dec_text="<div class='section_div_p10'>\n<div class='section_title' onclick=showOrHide('result_2')><a name='result_2'>Genetic variant (SNVs, INDELs, SVs)<font color='#56461f' > (click to show/hide) </font></a></div>\n<div id='result_2'>\n<div id='detection_result_2'>\n<table class='summary_table_2' style='width:800px'> \n"	
	HTML.write(dec_text)
	dec_text="".join(["<tr><td class='col1'>Number of SNVs / small INDELs / large SVs:</td><td class='col2'><B>", str(VCFdata['SNVcount']), " / ", str(VCFdata['INDELcount']), " / ", str(VCFdata['SVcount'])," <B></td></tr>\n"])
	HTML.write(dec_text)
	dec_text="</table>\n</div>\n\n"
	HTML.write(dec_text)
	
	dec_text="<div id='kmer_hits_figure_2'>\n<div class='figure' id='plot_kmer_hits_2' style='height:500px;width:800px'></div>\n</div>\n<script type=\"text/javascript\">\n"
	HTML.write(dec_text)
	dec_text="".join(["\tvar data=[{\n","\tx:",str(var_array),",\n","\ty:",str(VCFdata['density_del']),",\n","\tname: 'Percentages (%)  ',type:'bar',line:{color:'rgba(128,0,128,1.0)', width:1}\n","\t}];\n","\tvar layout={title:'Deletion size distribution', xaxis:{title:'Deletion size (bp)'}, yaxis:{title:'Percentages (%)'}};\n","\tPlotly.newPlot('plot_kmer_hits_2', data, layout);\n"])
	HTML.write(dec_text)
	dec_text="</script>\n\n<div id='kmer_hits_figure_3'>\n<div class='figure' id='plot_kmer_hits_3' style='height:500px;width:800px'></div>\n</div>\n<div id='kmer_hits_figure_3'>\n<script type=\"text/javascript\">\n"
	HTML.write(dec_text)
	dec_text="".join(["\tvar data=[{\n","\tx:",str(var_array),",\n","\ty:",str(VCFdata['density_ins']),",\n","\tname: 'Percentages (%)  ',type:'bar',line:{color:'rgba(128,0,128,1.0)', width:1}\n","\t}];\n","\tvar layout={title:'Insertion size distribution', xaxis:{title:'Insertion size (bp)'}, yaxis:{title:'Percentages (%)'}};\n","\tPlotly.newPlot('plot_kmer_hits_3', data, layout);\n"])
	HTML.write(dec_text)
	dec_text="</script>\n</div>\n</div>\n</div>\n</div>\n</div>\n</div>\n\n"
	HTML.write(dec_text)
	###end
	dec_text = "\n\n</body></html>"
	HTML.write(dec_text)
	HTML.close()
	FQ.close()
	
if __name__ == "__main__":
        main(sys.argv[1:])
