#!/bin/bash

# change permissions
chmod 755 *.py
cd src/EMA/
chmod 755 ema
cd ../long_fragment/
chmod 755 construct_fragment
cd ../
chmod 755 correct_barcode_stlfr
cd ..

###check SV calling tools: Aquila
if ! [ -x "$(command -v Aquila_step1)" ];
then
    echo 'Error: Aquila is not installed ...'
    echo 'using conda to install ...'
    conda install -c bioconda aquila
else
    echo 'using existing aquila ...'
fi

###check the alignment tool: bwa
if ! [ -x "$(command -v bwa)" ];
then
    echo 'Error: bwa is not installed ...'
    echo 'using conda to install ...'
    conda install -c bioconda bwa
else
    echo 'using existing bwa ...'
fi

###check the SNV calling tool:freebayes
if ! [ -x "$(command -v freebayes)" ];
then
    echo 'Error: freebayes is not installed ...'
    echo 'using conda to install ...'
    conda install -c bioconda freebayes
else
    echo 'using existing freebayes ...'
fi

###check the SNV calling tool: GATK
if ! [ -x "$(command -v gatk3)" ];
then 
    echo 'Error: gatk3 is not installed ...'
    echo 'using conda to install ...'
    conda install -c bioconda gatk
else
    echo 'using existing gatk ...'
fi

###check the SNV calling tool: HapCUT2
if ! [ -x "$(command -v hapcut2)" ];
then
    echo 'Error: hapcut2 is not installed ...'
    echo 'using conda to install ...'
    conda install -c bioconda hapcut2
else
    echo 'using existing HapCUT2 ...'
fi

###check the SV calling tool:LinkedSV
if ! [ -x "$(command -v linkedsv.py)" ];
then
    echo 'Error: LinkedSV is not installed ...'
    echo 'install LinkedSV from source ...'
    git clone https://github.com/WGLab/LinkedSV.git
    cd LinkedSV && sh build.sh
else
    echo 'using existing LinkedSV ...'
fi

###check the tool:parallel
if ! [ -x "$(command -v parallel)" ];
then
    echo 'Error: parallel is not installed ...'
    echo 'using conda to install ...'
    conda install -c bioconda parallel
else
    echo 'using existing parallel ...'
fi

###check the tool:picard
if ! [ -x "$(command -v picard)" ];
then
    echo 'Error: picard is not installed ...'
    echo 'using conda to install ...'
    conda install -c bioconda picard
else
    echo 'using existing picard ...'
fi

###check the tool:samtools
if ! [ -x "$(command -v samtools)" ];
then
    echo 'Error: samtools is not installed ...'
    echo 'using conda to install ...'
    conda install -c bioconda samtools
else
    echo 'using existing samtools ...'
fi

###check the phasing tool: SpecHap
if ! [ -x "$(command -v SpecHap)" ];
then
    echo 'Error: SpecHap is not installed ...'
    echo 'install SpecHap from source ...'
    git clone https://github.com/deepomicslab/SpecHap.git
    cd SpecHap && mkdir build && cd build && cmake .. && make && make install
else
    echo 'using existing SpecHap ...'
fi

###check the phasing tool: whatshap
if ! [ -x "$(command -v whatshap)" ];
then
    echo 'Error: whatshap is not installed ...'
    echo 'using conda to install whatshap ...'
    conda install -c bioconda whatshap
else
    echo 'using existing whatshap ...'
fi

###check the phasing tool: valor
if ! [ -x "$(command -v valor)" ];
then
    echo 'Error: valor is not installed ...'
    echo 'install valor from source ...'
    git clone https://github.com/BilkentCompGen/valor.git --recursive
    cd valor && make libs && make && cp -r ../valor ymp
else
    echo 'using existing valor ...'
fi

###Example
if ! [ -x "$(command -v vcfsort)" ];
then
    echo 'Error: vcflib is not installed ...'
    echo 'using conda to install ...'
    conda install -c bioconda vcflib
else
    echo 'using existing vcflib ...'
fi
