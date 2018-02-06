CGI_PIPELINES=/lustre/scratch117/cellgen/cellgeni/pipelines
version=1.4

RNAseq_path=$CGI_PIPELINES/RNAseq$version
PICARD_HOME=$RNAseq_path/picard-tools-2.0.1

# RNAseq

# Install FastQC
curl -fsSL http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip -o $RNAseq_path/fastqc_v0.11.5.zip
unzip $RNAseq_path/fastqc_v0.11.5.zip -d $RNAseq_path
chmod 755 $RNAseq_path/FastQC/fastqc
ln -s $RNAseq_path/FastQC/fastqc $RNAseq_path/bin/fastqc
rm $RNAseq_path/fastqc_v0.11.5.zip

# Install bedops
mkdir $RNAseq_path/bedops
curl -fsSL https://github.com/bedops/bedops/releases/download/v2.4.20/bedops_linux_x86_64-v2.4.20.v2.tar.bz2 -o $RNAseq_path/bedops_linux_x86_64-v2.4.20.v2.tar.bz2
tar xvjf $RNAseq_path/bedops_linux_x86_64-v2.4.20.v2.tar.bz2 -C $RNAseq_path/bedops
ln -s $RNAseq_path/bedops/bin/* $RNAseq_path/bin/
rm $RNAseq_path/bedops_linux_x86_64-v2.4.20.v2.tar.bz2

# Install cutadapt
virtualenv $RNAseq_path/venv
$RNAseq_path/venv/bin/pip install --install-option="--install-scripts=$RNAseq_path/bin" cutadapt==1.9.1

# Install TrimGalore
curl -fsSL http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/trim_galore_v0.4.2.zip -o $RNAseq_path/trim_galore_v0.4.2.zip
unzip $RNAseq_path/trim_galore_v0.4.2.zip -d $RNAseq_path/TrimGalore
ln -s $RNAseq_path/TrimGalore/trim_galore $RNAseq_path/bin/trim_galore
rm $RNAseq_path/trim_galore_v0.4.2.zip

# Install STAR
git clone https://github.com/alexdobin/STAR.git $RNAseq_path/STAR
ln -s $RNAseq_path/STAR/bin/Linux_x86_64/STAR $RNAseq_path/bin/STAR
ln -s $RNAseq_path/STAR/bin/Linux_x86_64/STARlong $RNAseq_path/bin/STARlong

# Install RSeQC
$RNAseq_path/venv/bin/pip install --install-option="--install-scripts=$RNAseq_path/bin" RSeQC==2.6.4

# Install SAMTools
curl -fsSL https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 -o $RNAseq_path/samtools-1.3.1.tar.bz2
tar xvjf $RNAseq_path/samtools-1.3.1.tar.bz2 -C $RNAseq_path/
cd $RNAseq_path/samtools-1.3.1
./configure --prefix=$RNAseq_path
make
make install
cd -
rm $RNAseq_path/samtools-1.3.1.tar.bz2

# Install PreSeq
curl -fsSL http://smithlabresearch.org/downloads/preseq_linux_v2.0.tar.bz2 -o $RNAseq_path/preseq_linux_v2.0.tar.bz2
tar xvjf $RNAseq_path/preseq_linux_v2.0.tar.bz2 -C $RNAseq_path/
ln -s $RNAseq_path/preseq_v2.0/preseq $RNAseq_path/bin/preseq
ln -s $RNAseq_path/preseq_v2.0/bam2mr $RNAseq_path/bin/bam2mr
rm $RNAseq_path/preseq_linux_v2.0.tar.bz2

# Install PicardTools
curl -fsSL https://github.com/broadinstitute/picard/releases/download/2.0.1/picard-tools-2.0.1.zip -o $RNAseq_path/picard-tools-2.0.1.zip
unzip $RNAseq_path/picard-tools-2.0.1.zip -d $RNAseq_path/
rm $RNAseq_path/picard-tools-2.0.1.zip

# Install R
RUN curl -fsSL https://cran.r-project.org/src/base/R-3/R-3.4.2.tar.gz -o /opt/R-3.4.2.tar.gz && \
    tar xvzf /opt/R-3.4.2.tar.gz -C /opt/ && \
    cd /opt/R-3.4.2 && \
    ./configure && \
    make && \
    make install && \
    rm /opt/R-3.4.2.tar.gz

# Install R Packages v2
RUN echo 'source("https://bioconductor.org/biocLite.R")' > /opt/packages.r && \
    echo 'biocLite()' >> /opt/packages.r && \
    echo 'biocLite(c("Rsubread", "dupRadar", "limma", "lattice", "locfit", "edgeR", "chron", "data.table", "gtools", "gdata", "bitops", "caTools", "gplots", "markdown"))' >> /opt/packages.r && \
    Rscript /opt/packages.r && \
    mkdir -p  /usr/local/lib/R/site-library

# Install featureCounts
curl -fsSL http://downloads.sourceforge.net/project/subread/subread-1.5.1/subread-1.5.1-Linux-x86_64.tar.gz -o $RNAseq_path/subread-1.5.1-Linux-x86_64.tar.gz
tar xvzf $RNAseq_path/subread-1.5.1-Linux-x86_64.tar.gz -C $RNAseq_path/
ln -s $RNAseq_path/subread-1.5.1-Linux-x86_64/bin/featureCounts $RNAseq_path/bin/featureCounts
rm $RNAseq_path/subread-1.5.1-Linux-x86_64.tar.gz

# Install StringTie
curl -fsSL http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.3.Linux_x86_64.tar.gz -o $RNAseq_path/stringtie-1.3.3.Linux_x86_64.tar.gz
tar xvzf $RNAseq_path/stringtie-1.3.3.Linux_x86_64.tar.gz -C $RNAseq_path/
ln -s $RNAseq_path/stringtie-1.3.3.Linux_x86_64/stringtie $RNAseq_path/bin/stringtie
rm $RNAseq_path/stringtie-1.3.3.Linux_x86_64.tar.gz

# Install MultiQC
$RNAseq_path/venv/bin/pip install --install-option="--install-scripts=$RNAseq_path/bin" git+git://github.com/ewels/MultiQC.git
