#!/bin/sh

RED='\033[1;31m'
NC='\033[0m' # No Color

downloadDir=$(pwd)
requiredSamtoolsVersion="1.3.1"
requiredBWAVersion="0.7.12"
requiredBedtoolsVersion="2.26.0"
requiredCircosVersion="0.69-3"

# apt-get or yum
# 1 == no
useAptGet=1
useYum=1

# System Prerequisites
needInstallGcc=1
needInstallLibncurses=1
needInstallZLib=1
needInstallGPlusPlus=1
needInstallMake=1
needInstallCPan=1

# Bioinformatics Tool Requirements
needInstallSamtools=1
needInstallBwa=1
needInstallBedtools=1
needInstallCircos=1

# Function definitions

check_installed_tools() {
    test -x "/usr/bin/apt-get"
    useAptGet=$?

    test -x "/usr/bin/yum"
    useYum=$?

    test -x "/usr/bin/gcc"; 
    needInstallGcc=$?

    test -x "/usr/bin/g++"; 
    needInstallGPlusPlus=$?

    test -x "/usr/bin/make";
    needInstallMake=$?

    test -x "/usr/bin/cpan";
    needInstallCPan=$?

    test -f "/usr/include/ncurses.h"; 
    needInstallLibncurses=$?

    test -f "/usr/include/zlib.h"; 
    needInstallZLib=$?

    test -x "/usr/local/bin/samtools"; 
    needInstallSamtools=$?

    test -x "/usr/local/bin/bwa"; 
    needInstallBwa=$?

    test -x "/usr/local/bin/bedtools"; 
    needInstallBedtools=$?

    test -d "/usr/local/circos";
    needInstallCircos=$?

    if [ "$needInstallGcc" = 1 ] || [ "$needInstallGPlusPlus" = 1 ] || [ "$needInstallMake" = 1 ] || [ "$needInstallLibncurses" = 1 ] || [ "$needInstallZLib" = 1 ] || [ "$needInstallCPan" = 1 ]; then
        install_system_prerequisites
    fi

    # Warn if any software is already installed (ie, needInstall[NameOfTool] = 0)
    if [ "$needInstallSamtools" = 0 ] || [ "$needInstallBwa" = 0 ] || [ "$needInstallBedtools" = 0 ] || [ "$needInstallCircos" = 0 ]; then
        print_installed_software
    fi
}

install_system_prerequisites() {    
    if [ "$useAptGet" = 1 ] && [ "$useYum" = 1 ]; then
        printf "\nFNBToolkit requires the following system prerequisites to be installed:\n\n"

        if [ "$needInstallGcc" = 1 ]; then
            printf "    gcc\n"
        fi

        if [ "$needInstallGPlusPlus" = 1 ]; then
            printf "    g++\n"
        fi

        if [ "$needInstallMake" = 1 ]; then
            printf "    make\n"
        fi

        if [ "$needInstallCPan" = 1 ]; then
            printf "    cpan\n"
        fi

        if [ "$needInstallLibncurses" = 1 ]; then
            printf "    libncurses5-dev (on Debian or Ubuntu Linux) or ncurses-devel (on RPM-based Linux distributions)\n"
        fi

        if [ "$needInstallZLib" = 1 ]; then
            printf "    zlib1g-dev (on Debian or Ubuntu Linux) or zlib-devel (on RPM-based Linux distributions)\n"
        fi

        printf "\n${RED}Please install these prerequisites and then re-run this installer.${NC}\n\n"

        exit 1
    else
        if [ "$useAptGet" = 0 ]; then
            sudo apt-get update
        fi
                
        if [ "$needInstallGcc" = 1 ]; then
            if [ "$useAptGet" = 0 ]; then
                sudo apt-get -y install gcc
            fi

            if [ "$useYum" = 0 ]; then
                sudo yum -y install gcc
            fi
        fi

        if [ "$needInstallGPlusPlus" = 1 ]; then
            if [ "$useAptGet" = 0 ]; then
                sudo apt-get -y install g++
            fi

            if [ "$useYum" = 0 ]; then
                sudo yum -y install gcc-c++
            fi
        fi

        if [ "$needInstallMake" = 1 ]; then
            if [ "$useAptGet" = 0 ]; then
                sudo apt-get -y install make
            fi

            if [ "$useYum" = 0 ]; then
                sudo yum -y install make
            fi
        fi

        if [ "$needInstallLibncurses" = 1 ]; then
            if [ "$useAptGet" = 0 ]; then
                sudo apt-get -y install libncurses5-dev
                sudo apt-get -y install ncurses-devel
            fi

            if [ "$useYum" = 0 ]; then
                sudo yum -y install libncurses5-dev
                sudo yum -y install ncurses-devel
            fi
        fi

        if [ "$needInstallZLib" = 1 ]; then
            if [ "$useAptGet" = 0 ]; then
                sudo apt-get -y install zlib1g-dev
                sudo apt-get -y install zlib-devel
            fi

            if [ "$useYum" = 0 ]; then
                sudo yum -y install zlib1g-dev
                sudo yum -y install zlib-devel
            fi
        fi

        if [ "$needInstallCPan" = 1 ]; then
            if [ "$useAptGet" = 0 ]; then
                sudo apt-get -y install perl-CPAN
            fi

            if [ "$useYum" = 0 ]; then
                sudo yum -y install perl-CPAN
            fi
        fi
    fi
}

print_installed_software() {
    printf "\nFNBToolkit requires several Bioinformatics tools to be installed in order to run correctly.\nWe have detected that you have one or more of these programs installed already."
    printf "\n${RED}This installer will overwrite the following tools that you already have installed:${NC}\n\n"

    if [ "$needInstallSamtools" = 0 ] && [ "$(samtools --version | head -n1 | cut -d" " -f2)" != "$requiredSamtoolsVersion" ]; then
        printf "    Samtools\n"
    fi

    if [ "$needInstallBwa" = 0 ]; then
        printf "    BWA\n"
    fi

    if [ "$needInstallBedtools" = 0 ] && [ "$(bedtools --version | head -n1 | cut -d" " -f2)" != "v$requiredBedtoolsVersion" ]; then
        printf "    Bedtools\n"
    fi

    if [ "$needInstallCircos" = 0 ]; then
        printf "    Circos\n"
    fi

    printf "\n${RED}Do you want to continue and overwrite your current tools with the versions listed?${NC} (y/n):"

    while true; do
        read -p "" yn
        case $yn in
            [Yy]*) install_prerequisites; exit;;
            [Nn]*) printf "Unable to install FNB Toolkit.\n"; exit;;
            *) printf "Please answer yes or no.";;
        esac
    done
}

install_samtools() {
    printf "\nDownloading Samtools ${requiredSamtoolsVersion}\n\n"
    mkdir /tmp/sam-tools-${requiredSamtoolsVersion}
    cd /tmp/sam-tools-${requiredSamtoolsVersion} && wget https://github.com/samtools/samtools/releases/download/${requiredSamtoolsVersion}/samtools-${requiredSamtoolsVersion}.tar.bz2
    
    printf "\nInstalling Samtools ${requiredSamtoolsVersion}\n\n"
    sudo bzip2 -d samtools-${requiredSamtoolsVersion}.tar.bz2 && tar -xf samtools-${requiredSamtoolsVersion}.tar
    cd samtools-${requiredSamtoolsVersion}/ && ./configure

    if [ "$?" -ne 0 ]; then
        printf "${RED}There was a problem with samtools configuration.\n\nPlease review samtools-config.log for more information and then retry the installation.${NC}\n"
        cd $downloadDir
        sudo cp /tmp/sam-tools-${requiredSamtoolsVersion}/samtools-${requiredSamtoolsVersion}/config.log ./samtools-config.log
        clean
        exit -1
    fi

    sudo make
    sudo make install

    printf "\nSamtools successfully installed.\n"
}

install_bwa() {
    printf "\nDownloading BWA ${requiredBWAVersion}\n\n"
    mkdir /tmp/bwa-${requiredBWAVersion}
    cd /tmp/bwa-${requiredBWAVersion} && wget https://github.com/lh3/bwa/archive/0.7.12.tar.gz
    
    printf "\nInstalling BWA ${requiredBWAVersion}\n\n"
    sudo gunzip ${requiredBWAVersion}.tar.gz && tar -xf ${requiredBWAVersion}.tar
    cd bwa-${requiredBWAVersion}
    sudo make
    sudo cp bwa /usr/local/bin

    printf "\nBWA successfully installed.\n"
}

install_bedtools() {
    printf "\nDownloading Bedtools ${requiredBedtoolsVersion}\n\n"
    mkdir /tmp/bedtools-${requiredBedtoolsVersion}
    cd /tmp/bedtools-${requiredBedtoolsVersion} && wget https://github.com/arq5x/bedtools2/releases/download/v${requiredBedtoolsVersion}/bedtools-${requiredBedtoolsVersion}.tar.gz
    
    printf "\nInstalling Bedtools ${requiredBedtoolsVersion}\n\n"
    sudo gunzip bedtools-${requiredBedtoolsVersion}.tar.gz && tar -xf bedtools-${requiredBedtoolsVersion}.tar
    cd bedtools2
    sudo make
    sudo cp ./bin/bedtools /usr/local/bin

    printf "\nBedtools successfully installed.\n"
}

install_circos() {
    printf "\nDownloading Circos ${requiredCircosVersion}\n\n"
    mkdir /usr/local/circos
    cd /usr/local/circos && wget http://circos.ca/distribution/circos-${requiredCircosVersion}.tgz
    
    printf "\nInstalling Circos ${requiredCircosVersion}\n\n"
    sudo tar xvfz circos-${requiredCircosVersion}.tgz
    sudo ln -s circos-${requiredCircosVersion} current

    if [ "$useYum" = 0 ]; then
        sudo yum install -y perl-Module-Build
    fi

    cpan YAML

    perl -MCPAN -e 'install Clone'
    perl -MCPAN -e 'install Config::General'
    perl -MCPAN -e 'install Digest::MD5'
    perl -MCPAN -e 'install File::Spec::Functions'
    perl -MCPAN -e 'install Bundle::Compress::Zlib'
    perl -MCPAN -e 'install Font::TTF::Font'
    perl -MCPAN -e 'install List::MoreUtils'
    perl -MCPAN -e 'install Math::Bezier'
    perl -MCPAN -e 'install Math::Round'
    perl -MCPAN -e 'install Math::VecStat'
    perl -MCPAN -e 'install Params::Validate'
    perl -MCPAN -e 'install Readonly'
    perl -MCPAN -e 'install Regexp::Common'
    perl -MCPAN -e 'install SVG'
    perl -MCPAN -e 'install Set::IntSpan'
    perl -MCPAN -e 'install Statistics::Basic'
    perl -MCPAN -e 'install Text::Balanced'
    perl -MCPAN -e 'install Text::Format'

    if [ "$useAptGet" = 0 ]; then
        sudo apt-get install -y libgd-gd2-perl
    fi
    
    if [ "$useYum" = 0 ]; then
        sudo yum install -y libgd-gd2-perl
        sudo yum install -y gd gd-devel perl-GD
    fi
    perl -MCPAN -e 'install GD::Polyline'

    printf "\nCircos successfully installed.\n"
}

install_fnb_main() {
    printf "Copying FNBTools files to /usr/local/bin\n"
    cd $downloadDir
    sudo cp fnbalign.pl /usr/local/bin/fnbalign
    sudo chmod 755 /usr/local/bin/fnbalign
    printf "\nfnbalign copied!\n"

    cd $downloadDir
    sudo cp fnbscan.pl /usr/local/bin/fnbscan
    sudo chmod 755 /usr/local/bin/fnbscan
    printf "\nfnbscan copied!\n"

    cd $downloadDir
    sudo cp circosvis.py /usr/local/bin/circosvis
    sudo chmod 755 /usr/local/bin/circosvis
    printf "\ncircosvis copied!\n"

    cd $downloadDir
    if [ ! -d "/usr/local/fnbtools" ]; then
        sudo mkdir /usr/local/fnbtools
    fi

    sudo cp -r ./ /usr/local/fnbtools

    sudo rm /usr/local/fnbtools/Install.sh
    sudo rm /usr/local/fnbtools/fnbtools.tar.gz

    printf "\nFNBTools dependencies copied!\n"
}

clean() {
    sudo rm -rf /tmp/sam-tools-${requiredSamtoolsVersion}/
    sudo rm -rf /tmp/bwa-${requiredBWAVersion}
    sudo rm -rf /tmp/bedtools-${requiredBedtoolsVersion}
    sudo rm -rf /tmp/circos-${requiredCircosVersion}
}

install_prerequisites() {
    if [ "$(samtools --version | head -n1 | cut -d" " -f2)" != "$requiredSamtoolsVersion" ]; then
        install_samtools
    fi

    install_bwa

    if [ "$(bedtools --version | head -n1 | cut -d" " -f2)" != "v$requiredBedtoolsVersion" ]; then
        install_bedtools
    fi

    install_circos

    clean
}

# Main installation sequence
check_installed_tools
install_prerequisites
install_fnb_main
printf "\n\nFNBTools successfully installed -- Please continue with Step 5 of the README installation section to enable Circos visualizations of FNBTools outputs.\n\n"