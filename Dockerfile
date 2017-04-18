FROM ubuntu:trusty
RUN apt-get update && apt-get install -y wget
RUN apt-get install -y gcc
RUN apt-get install -y g++
RUN apt-get install -y make
RUN apt-get install -y libncurses5-dev
RUN apt-get install -y zlib1g-dev
RUN apt-get install -y python-minimal
RUN mkdir /tmp/sam-tools-1.3.1
RUN cd /tmp/sam-tools-1.3.1 && wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
RUN cd /tmp/sam-tools-1.3.1 && bzip2 -d samtools-1.3.1.tar.bz2 && tar -xf samtools-1.3.1.tar
RUN cd /tmp/sam-tools-1.3.1/samtools-1.3.1/ && ./configure && make && make install

RUN mkdir /tmp/bwa-0.7.12
RUN cd /tmp/bwa-0.7.12 && wget https://github.com/lh3/bwa/archive/0.7.12.tar.gz    
RUN cd /tmp/bwa-0.7.12 && gunzip 0.7.12.tar.gz && tar -xf 0.7.12.tar
RUN cd /tmp/bwa-0.7.12/bwa-0.7.12 && make
RUN cp /tmp/bwa-0.7.12/bwa-0.7.12/bwa /usr/local/bin


RUN mkdir /tmp/bedtools-2.26.0
RUN cd /tmp/bedtools-2.26.0 && wget https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz
RUN cd /tmp/bedtools-2.26.0 && gunzip bedtools-2.26.0.tar.gz && tar -xf bedtools-2.26.0.tar
RUN cd /tmp/bedtools-2.26.0/bedtools2 && make
RUN cd /tmp/bedtools-2.26.0/bedtools2/bin && cp bedtools /usr/local/bin


RUN mkdir /tmp/bcftools-1.3.1
RUN cd /tmp/bcftools-1.3.1 && wget https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2
RUN cd /tmp/bcftools-1.3.1 && bzip2 -d bcftools-1.3.1.tar.bz2 && tar -xf bcftools-1.3.1.tar
RUN cd /tmp/bcftools-1.3.1/bcftools-1.3.1 && make
RUN cd /tmp/bcftools-1.3.1/bcftools-1.3.1 && cp bcftools /usr/local/bin

RUN mkdir /usr/local/circos
RUN cd /usr/local/circos && wget http://circos.ca/distribution/circos-0.69-3.tgz
RUN cd /usr/local/circos && tar xvfz circos-0.69-3.tgz
RUN cd /usr/local/circos && ln -s circos-0.69-3 current
ENV PATH /usr/local/circos/current/bin:$PATH
ENV . ~/.bashrc

RUN perl -MCPAN -e 'install Clone'
RUN perl -MCPAN -e 'install Config::General'
RUN perl -MCPAN -e 'install Digest::MD5'
RUN perl -MCPAN -e 'install File::Spec::Functions'
RUN perl -MCPAN -e 'install Font::TTF::Font'
RUN perl -MCPAN -e 'install List::MoreUtils'
RUN perl -MCPAN -e 'install Math::Bezier'
RUN perl -MCPAN -e 'install Math::Round'
RUN perl -MCPAN -e 'install Math::VecStat'
RUN perl -MCPAN -e 'install Params::Validate'
RUN perl -MCPAN -e 'install Readonly'
RUN perl -MCPAN -e 'install Regexp::Common'
RUN perl -MCPAN -e 'install SVG'
RUN perl -MCPAN -e 'install Set::IntSpan'
RUN perl -MCPAN -e 'install Statistics::Basic'
RUN perl -MCPAN -e 'install Text::Balanced'
RUN perl -MCPAN -e 'install Text::Format'
RUN apt-get install -y libgd-gd2-perl
RUN perl -MCPAN -e 'install GD::Polyline'