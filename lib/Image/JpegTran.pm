package Image::JpegTran;

use 5.008008;
use strict;
use warnings;
use base 'Exporter';
use Carp;

our @EXPORT_OK = our @EXPORT = qw( jpegtran );
our $VERSION = '0.01';

use XSLoader;
XSLoader::load('Image::JpegTran', $VERSION);

sub jpegtran($$;%) {
	 my $src = shift;
	-e( $src ) or croak "Can't find source file `$src'";
	my $dst = shift;
	my %args = (
		@_==1 && ref $_[0] ? %{$_[0]} : @_
	);
	_jpegtran($src,$dst,\%args);
}


1;
__END__
__END__
=head1 NAME

Image::JpegTran - XS wrapper around lossless JPEG transformation utility - jpegtran

=head1 SYNOPSIS

    use Image::JpegTran;
    
    jpegtran 'source.jpg','result.jpg', -rotate => 90, -trim, -perfect;

=head1 DESCRIPTION

=head1 AUTHOR

Mons Anderson, E<lt>mons@cpan.orgE<gt>

=head1 COPYRIGHT AND LICENSE

The main part of this module is copyright (C) 1991-2010

The Independent JPEG Group's JPEG software

Thomas G. Lane, Guido Vollbeding.

See README.IJG

=cut
