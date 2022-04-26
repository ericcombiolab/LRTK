#!/bin/bash
# change permission

chmod 755 *.py
cd src/EMA/
chmod 755 ema
cd ../long_fragment/
chmod 755 construct_fragment
cd ../
chmod 755 correct_barcode_stlfr
cd ..

if ! [ -x "$(command -v samtools)" ];
then
    echo 'Error: samtools is not installed...'
    echo 'using Conda to install...'
    conda install -c bioconda samtools
else
    echo 'using existing samtools...'
fi

if ! [ -x "$(command -v picard)" ];
then
    echo 'Error: picard is not installed...'
    echo 'using Conda to picard...'
    conda install -c bioconda picard
else
    echo 'using existing picard...'
fi

if ! [ -x "$(command -v bwa)" ];
then
    echo 'Error: bwa is not installed...'
    echo 'using Conda to install...'
    conda install -c bioconda bwa
else
    echo 'using existing bwa...'
fi
if ! [ -x "$(command -v freebayes)" ];
then
    echo 'Error: freebayes is not installed...'
    echo 'using Conda to install...'
    conda install -c bioconda freebayes
else
    echo 'using existing freebayes...'
fi

if ! [ -x "$(command -v HAPCUT2)" ];
then
    echo 'Error: HAPCUT2 is not installed...'
    echo 'using Conda to install...'
    conda install -c bioconda hapcut2
else
    echo 'using existing HAPCUT2...'
fi

if ! [ -x "$(command -v Aquila_step1)" ];
then
    echo 'Error: Aquila is not installed...'
    echo 'using Conda to install...'
    conda install -c bioconda aquila
else
    echo 'using existing Aquila...'
fi


echo 'You have installed LRTK dependencies successfully!'
