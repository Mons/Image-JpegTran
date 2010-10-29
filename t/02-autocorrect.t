use strict;
use warnings;
use bytes();
use uni::perl ':dumper';

use Test::More tests => 97;

use Image::JpegTran ':all';
use Image::LibExif;
use File::Copy;

our $DIR = "t/fx";

sub diff($$) {
	use bytes;
	no utf8;
	my ($srcf,$dstf) = @_;
	my $src = do { open my $f, '<:raw',$srcf or die "$srcf: $!"; local $/; <$f> };
	my $dst = do { open my $f, '<:raw',$dstf or die "$dstf: $!"; local $/; <$f> };
	return 1 if $src eq $dst;
	my $sl = length $src;
	my $dl = length $dst;
	my $i = 0;
	while(substr($src,$i,1) eq substr($dst,$i,1)) {
		$i++;
	}
	print "# diff at byte $i from the beginning\n";
	printf "\t".('%02x 'x10)."\n", unpack('C10', substr($src,$i - 3,10));
	printf "\t".('%02x 'x10)."\n", unpack('C10', substr($dst,$i - 3,10));

	my $e = 1;
	while(substr($src,$sl-$e,1) eq substr($dst,$dl - $e,1)) {
		$e++;
	}
	print "# diff at byte $e from the end (src:@{[ $sl - $e ]} / dst: @{[ $dl - $e ]})\n";
	printf "\t".('%02x 'x10)."\n", unpack('C10', substr($src,$sl - $e - 2,10));
	printf "\t".('%02x 'x10)."\n", unpack('C10', substr($dst,$dl - $e - 2,10));
	if (( my $x = $sl - $e - $i ) < 200) {
		my $y = $dl - $e - $i;
		$x+=4;
		printf "\t".('%02x 'x$x)."\n", unpack('C'.$x, substr($src,$i - 2,$x));
		printf "\t".('%02x 'x$y)."\n", unpack('C'.$y, substr($dst,$i - 2,$y));
	}
	my $sexif = image_exif($srcf); my $dexif = image_exif($dstf);
	my %uniq;
	warn dumper $sexif,$dexif;
	for ( sort grep { !$uniq{$_}++ } keys %$sexif, keys %$dexif ) {
		if (exists $sexif->{$_} and exists $dexif->{$_}) {
			if (ref $sexif->{$_}) {
				if (${$sexif->{$_}} ne ${$dexif->{$_}}) {
					warn "$_: x <> x\n";
				}
			} else {
				if ($sexif->{$_} ne $dexif->{$_}) {
					warn "$_: $sexif->{$_} <> $dexif->{$_}\n";
				}
			}
		}
		elsif (!exists $sexif->{$_}) {
			warn "no in src: $_: $dexif->{$_}\n";
		}
		else {
			warn "no in dst: $_: $sexif->{$_}\n";
			
		}
	}
	return 0;
}
#my $orig = do { open my $f, '<:raw',"$DIR/X.jpg"; local $/; <$f> };

for my $file (<$DIR/F*.jpg>) {
	( my $dst = $file ) =~ s{$DIR/F}{$DIR/O};
	if( my $rc = jpegautotran($file,$dst,{ perfect => 1, copy => 'exif', discard_thumbnail => 1 }) ) {
		diag "copy:exif, nothumb: $file";
		my $exif = image_exif($dst);
		is $exif->{ImageWidth},  80, 'width';
		is $exif->{ImageLength}, 112, 'height';
		is $exif->{Orientation}+0, 1, 'orientation';
		ok !defined $exif->{ThumbnailImage}, 'no thumb';
		#diff("$DIR/X.jpg", $dst);
		#system("diff $DIR/X.jpg $dst");
		unlink $dst unless $dst =~ /O\.jpg$/;
	} else {
		say "$file - skip";
	}
}

for my $file (<$DIR/F*.jpg>) {
	( my $dst = $file ) =~ s{$DIR/F}{$DIR/O};
	if( my $rc = jpegautotran($file,$dst,{ perfect => 1, discard_thumbnail => 0 }) ) {
		diag "thumb: $file";
		my $exif = image_exif($dst);
		is $exif->{ImageWidth},  80, 'width';
		is $exif->{ImageLength}, 112, 'height';
		is $exif->{Orientation}+0, 1, 'orientation';
		ok ref $exif->{ThumbnailImage}, 'have thumb';
		ok diff("$DIR/F.jpg", $dst), 'no diff';
		unlink $dst;
	} else {
		diag "$file - skip";
	}
}

for my $file (<$DIR/F*.jpg>) {
	( my $dst = $file ) =~ s{$DIR/F}{$DIR/I};
	copy $file,$dst or die "$!";
}

for my $file (<$DIR/I*.jpg>) {
	if( my $rc = jpegautotran($file,undef,{ perfect => 1, discard_thumbnail => 1, copy => 'none' }) ) {
		diag "inplace: copy:none, nothumb: $file";
		my $exif = image_exif($file);
		ok !$exif, 'no exif';
		ok diff("$DIR/N.jpg", $file), 'no diff';
		unlink $file;
	} else {
		diag "$file - skip";
	}
}

{
	my $orig = "$DIR/X.jpg";
	my $file = "$DIR/Z.jpg";
	copy $orig, $file or die $!;

	my $exif = image_exif($file);
	is $exif->{Orientation}+0,1,"before";
	
	for my $tran (
		{ rotate => 90,    o => 8 },
		{ rotate => 180,   o => 6, },
		{ rotate => 270,   o => 3, },
		{ transverse => 1, o => 5, },
		{ flip => 'h',     o => 8, },
		{ transpose => 1,  o => 4, },
		{ flip => 'v'  ,   o => 1, },
	) {
		my $o = delete $tran->{o};
		jpegtran($file,undef,{ perfect => 1, %$tran });
		my $exif = image_exif($file);
		is $exif->{Orientation}+0,$o,"after @{[ %$tran ]}";
	}
	
	ok diff($orig, $file), 'no diff for all trans';
	unlink $file;
	
}