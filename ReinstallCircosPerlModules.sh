#!/bin/sh

useAptGet=1
useYum=1

test -x "/usr/bin/apt-get"
useAptGet=$?

test -x "/usr/bin/yum"
useYum=$?

printf "\nInstalling Circos required Perl modules...\n\n"

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

printf "\nInstallation complete.  Please run circos -modules to verify that all modules were installed correctly.\n\nIf any are missing, re-run this script with sudo ./ReinstallCircosPerlModules.sh to re-try the installation."
