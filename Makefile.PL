use lib qw(inc);
use Devel::CheckLib;
use Cwd;
use Config;
use ExtUtils::PkgConfig;

my $libexif = 'libexif';

ExtUtils::PkgConfig->atleast_version($libexif,'0.6.18') or die "N/A: Need libexif >= 0.6.18\n";
my $include = ExtUtils::PkgConfig->cflags_only_I($libexif);
my $cflags  = ExtUtils::PkgConfig->cflags_only_other($libexif);
my $libs    = ExtUtils::PkgConfig->libs($libexif);

my @ll = do {
	my %uniq;
	grep length && !$uniq{$_}++,
	( split /\s+/, $Config::Config{libsdirs}),
	( split /\s+/, $Config::Config{loclibpth}),
};
my @ii = do {
	my %uniq;
	grep length && !$uniq{$_}++,
	( split /\s+/, $Config::Config{usrinc}),
	( split /\s+/, $Config::Config{locincpth}),
};


unshift @ll,'/usr/local/lib';
print "# LL = [@ll]\n";
print "# II = [@ii]\n";

check_lib_or_exit(
    debug   => 1,
    lib     => 'jpeg',
    libpath => [ @ll ],
    incpath => [ Cwd::cwd(), @ii ],
    header  => [['"jinclude.h"','<jpeglib.h>']],
    define0 => "#define JPEG_INTERNALS\n",
    define  => "#if JPEG_LIB_VERSION < 80\n#error jpeglib version >= 80 required\n#endif",
);

check_lib_or_exit(
    debug   => 1,
    LIBS    => $libs,
    INC     => ExtUtils::PkgConfig->cflags($libexif),
    header  => [['<libexif/exif-data.h>']],
    #header  => [['"jinclude.h"','<jpeglib.h>']],
    #define0 => "#define JPEG_INTERNALS\n",
    #define  => "#if JPEG_LIB_VERSION < 80\n#error jpeglib version >= 80 required\n#endif",
);

$include .= ' '.join(' ',map "-I$_", @ii);
$libs    .= ' '.join(' ',map "-L$_", @ll);

#perl -Muni::perl=:dumper -MConfig -E 'say dumper \%Config::Config'                

use 5.008008;
use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
	NAME              => 'Image::JpegTran',
	VERSION_FROM      => 'lib/Image/JpegTran.pm', # finds $VERSION
	PREREQ_PM         => {}, # e.g., Module::Name => 1.1
	ABSTRACT_FROM     => 'lib/Image/JpegTran.pm', # retrieve abstract from module
	AUTHOR            => 'Mons Anderson <mons@cpan.org>',
	LIBS              => ['-ljpeg '.$libs], # e.g., '-lm'
	DEFINE            => '', # e.g., '-DHAVE_SOMETHING'
	INC               => '-I. '.$include, # e.g., '-I. -I/usr/include/other'
	CCFLAGS           => $cflags,
	#CCFLAGS           => '-Wunused -std=c99',
	OBJECT            => 'JpegTran.o transupp.o', # link all the C files too
	META_ADD          => { license => 'perl' },
);
