#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include "ppport.h"

#define ENTROPY_OPT_SUPPORTED  1

#include "jinclude.h"
#include <jpeglib.h>
#include "transupp.h"
#include <stdio.h>
#include <stdlib.h>

#include <libexif/exif-loader.h>
#include <libexif/exif-data.h>
#include <libexif/exif-log.h>
#include <libexif/exif-mem.h>

#define DES_DECO_SRC 0x0001
#define DES_COMP_DST 0x0002

#define my_croak(...) \
	STMT_START { \
		if (fp != NULL) { fclose(fp); } \
		if (SvIV(what_destroy) & DES_DECO_SRC) { \
			/*(void) jpeg_finish_decompress(&srcinfo);*/ \
			jpeg_destroy_decompress(&srcinfo); \
		}\
		if (SvIV(what_destroy) & DES_COMP_DST) { \
			/*jpeg_finish_compress(&dstinfo);*/ \
			jpeg_destroy_compress(&dstinfo); \
		}\
		croak(__VA_ARGS__);\
	} STMT_END

static unsigned char autocorrect[] = {
	 0                   // skip
	, 0                  //normal, do nothing
	,JXFORM_FLIP_H       //flipped h
	,JXFORM_ROT_180      //rotated 180
	,JXFORM_FLIP_V       //flipped v
	,JXFORM_TRANSPOSE    //transposed
	,JXFORM_ROT_90       //rotated 270
	,JXFORM_TRANSVERSE   //transversed
	,JXFORM_ROT_270      //rotated 90
};

typedef struct  {
	unsigned int         copy;
	jpeg_transform_info  trans;
	long                 max_memory_to_use;
	unsigned char        optimize_coding;
	unsigned char        arith_code;
	unsigned char        progressive;
} jpegtran_config;

			/*
			char *opt = SvPV_nolen(*key);
			char *begin;
			STRLEN len;
			begin = opt;
			while (*opt) {
				if (*opt == ' ' || *opt == ';' || *opt == ',') {
					if (begin == opt) { begin++; opt++; continue; }
					warn("option: %*s",)
				}
				opt++;
			}
			*/

static char * parse_config(HV* conf, jpegtran_config * config) {
	SV **key;
	memzero(config, sizeof(*config));
	
	config->trans.transform = JXFORM_NONE;
	
	if ((key = hv_fetch(conf, "perfect", 7, 0)) && SvTRUE(*key)) {
		config->trans.perfect = TRUE;
	} else {
		config->trans.perfect = FALSE;
	}
	
	if ((key = hv_fetch(conf, "trim", 4, 0)) && SvTRUE(*key)) {
		config->trans.trim    = TRUE;
	} else {
		config->trans.trim    = FALSE;
	}
	
	if ((key = hv_fetch(conf, "grayscale", 9, 0)) && SvTRUE(*key)) {
		config->trans.force_grayscale = TRUE;
	} else {
		config->trans.force_grayscale = FALSE;
	}

	if (key = hv_fetch(conf, "rotate", 6, 0)) {
		if( config->trans.transform != JXFORM_NONE) { croak("Can't apply several transforms at once"); }
		if (SvIOK( *key )) {
			config->trans.transform = 
				SvIV(*key) ==  90 ? JXFORM_ROT_90 :
				SvIV(*key) == 180 ? JXFORM_ROT_180 :
				SvIV(*key) == 270 ? JXFORM_ROT_270 :
				JXFORM_NONE;
			if( config->trans.transform == JXFORM_NONE ) {
				croak("Bad value for rotate");
			}
		} else {
			croak("Bad value for rotate");
		}
	}
	if ((key = hv_fetch(conf, "transpose", 9, 0)) && SvTRUE(*key)) {
		if( config->trans.transform != JXFORM_NONE){ croak("Can't apply several transforms at once"); }
		config->trans.transform = JXFORM_TRANSPOSE;
	}
	if ((key = hv_fetch(conf, "transverse", 10, 0)) && SvTRUE(*key)) {
		if( config->trans.transform != JXFORM_NONE){ croak("Can't apply several transforms at once"); }
		config->trans.transform = JXFORM_TRANSVERSE;
	}
	if (key = hv_fetch(conf, "flip", 4, 0)) {
		if( config->trans.transform != JXFORM_NONE){ croak("Can't apply several transforms at once"); }
		if (SvPOK( *key )) {
			if (strEQ(SvPV_nolen(*key),"horizontal") || strEQ(SvPV_nolen(*key),"horisontal") || strEQ(SvPV_nolen(*key),"h")) {
				config->trans.transform = JXFORM_FLIP_H;
			} else
			if (strEQ(SvPV_nolen(*key),"vertical") || strEQ(SvPV_nolen(*key),"v")) {
				config->trans.transform = JXFORM_FLIP_V;
			} else
			{
				croak("Bad value for flip. Could be [ horizontal | vertical | h | v ]");
			}
		} else {
			croak("Bad value for flip. Could be [ horizontal | vertical | h | v ]");
		}
	}
	if (key = hv_fetch(conf, "discard_thumbnail", 17, 0)) {
		if (SvTRUE(*key)) {
			config->trans.discard_thumbnail = TRUE;
		} else {
			config->trans.discard_thumbnail = FALSE;
		}
	} else {
		config->trans.discard_thumbnail = TRUE;
	}
	
	if (key = hv_fetch(conf, "copy", 4, 0)) {
		if (SvPOK(*key)) {
			char *copyopt = SvPV_nolen(*key);
			// none comments exif all
			if (!strcmp(copyopt,"none")) {
				config->copy = JCOPY_NONE;
			} else
			if (!strcmp(copyopt,"comments")) {
				config->copy = JCOPY_COM;
			} else 
			if (!strcmp(copyopt,"exif")) {
				config->copy = JCOPY_EXIF;
			} else 
			if (!strcmp(copyopt,"all")) {
				config->copy = JCOPY_ALL;
			} else
			{
				croak ("Bad value for copy. Available are: [ none | exif | comments | all ]");
			}
		} else {
			croak ("Bad value for copy. Available are: [ none | exif | comments | all ]");
		}
	} else {
		config->copy = JCOPY_ALL;
	}
	
	if ( key = hv_fetch(conf, "optimize", 8, 0) ) {
		// Enable entropy parm optimization.
		config->optimize_coding = SvTRUE(*key) ? TRUE : FALSE;
	} else {
		config->optimize_coding = TRUE;
	}
	
	if ((key = hv_fetch(conf, "arithmetic", 10, 0)) && SvTRUE(*key)) {
		config->arith_code = SvTRUE(*key) ? TRUE : FALSE;
	} else {
		config->arith_code = FALSE;
	}
	
	if ((key = hv_fetch(conf, "maxmemory", 9, 0))) {
		if (SvIOK(*key)) {
			config->max_memory_to_use = SvIV(*key);
		} else {
			croak("Bad value for maxmemory. Should be integer");
		}
	}
	if ((key = hv_fetch(conf, "progressive", 11, 0)) && SvTRUE(*key)) {
		config->progressive = TRUE;
	}
	//  TODO
	//  if (! jtransform_parse_crop_spec(&transformoption, argv[argn]));
	return 0;
}

typedef struct {
	char *src;
	char *dst;
	
} jpegtran_pair;


#define BY_REF 1
#define BY_NAME 2

struct my_error_mgr {
	struct jpeg_error_mgr pub;    /* "public" fields */
	jmp_buf setjmp_buffer;        /* for return to caller */
	int des_flag;
};

typedef struct my_error_mgr * my_error_ptr;

METHODDEF(void)
my_error_exit (j_common_ptr cinfo)
{
	my_error_ptr myerr = (my_error_ptr) cinfo->err;
	(*cinfo->err->output_message) (cinfo);
	longjmp(myerr->setjmp_buffer, 1);
}

METHODDEF(void)
my_output_message (j_common_ptr cinfo)
{
	char buffer[JMSG_LENGTH_MAX];
	(*cinfo->err->format_message) (cinfo, buffer);
	warn("JPEG Library Error: %s",buffer);
}

void jpegtran_execute (SV *src, SV *dst, char *src_name, char *dst_name, int src_arg_type, int dst_arg_type, jpegtran_config *config) {
		unsigned char *buf=NULL;
		unsigned long buf_size=0;
		
		struct jpeg_decompress_struct   srcinfo;
		struct jpeg_compress_struct     dstinfo;
		struct my_error_mgr             jsrcerr,jdsterr;
		jvirt_barray_ptr * src_coef_arrays;
		jvirt_barray_ptr * dst_coef_arrays;
		char buffer[JMSG_LENGTH_MAX];
		
		SV *what_destroy = sv_2mortal(newSViv(0));

		// We assume all-in-memory processing and can therefore use only a
		// single file pointer for sequential input and output operation.
		FILE * fp = NULL;

		srcinfo.err = jpeg_std_error(&jsrcerr.pub);
		jsrcerr.pub.error_exit = my_error_exit;
		jsrcerr.pub.output_message = my_output_message;

		sv_setiv(what_destroy,SvIV(what_destroy) | DES_DECO_SRC );
		if (setjmp(jsrcerr.setjmp_buffer)) {
			(*srcinfo.err->format_message) (&srcinfo, buffer);
			// warn("wd = %d; ", SvIV(what_destroy));
			my_croak("Real shit happens in src: %s",buffer);
		}
		jpeg_create_decompress(&srcinfo);

		dstinfo.err = jpeg_std_error(&jdsterr.pub);
		jdsterr.pub.error_exit = my_error_exit;
		jdsterr.pub.output_message = my_output_message;
		sv_setiv(what_destroy,SvIV(what_destroy) | DES_COMP_DST );
		if (setjmp(jdsterr.setjmp_buffer)) {
			(*dstinfo.err->format_message) (&dstinfo, buffer);
			my_croak("Real shit happens in dst: %s",buffer);
		}
		jpeg_create_compress(&dstinfo);

		dstinfo.optimize_coding = config->optimize_coding;
		dstinfo.arith_code      = config->arith_code;
		if ( config->max_memory_to_use )
			dstinfo.mem->max_memory_to_use = config->max_memory_to_use;

		if( src_arg_type == BY_NAME ) {
			if ((fp = fopen(src_name, READ_BINARY)) == NULL) {
				my_croak("can't open `%s' for reading",src_name);
			}

			jpeg_stdio_src(&srcinfo, fp);
		} else {
			// working with memory data

			SV * str_r;
			str_r = SvRV( src );
			STRLEN len;
			char *srcbin  = SvPVbyte( str_r,len );

			jpeg_mem_src(&srcinfo, srcbin, len);
		}
		// Enable saving of extra markers that we want to copy
		jcopy_markers_setup(&srcinfo, config->copy);
		
		// Read file header
		if (!jpeg_read_header(&srcinfo, TRUE)) {
			my_croak("bad source");
		}
		
		if (!jtransform_request_workspace(&srcinfo, &config->trans)) {
			my_croak("transformation is not perfect");
		}
		
		// Read source file as DCT coefficients
		src_coef_arrays = jpeg_read_coefficients(&srcinfo);
		
		// Initialize destination compression parameters from source values
		jpeg_copy_critical_parameters(&srcinfo, &dstinfo);
		
		// Adjust destination parameters if required by transform options;
		// also find out which set of coefficient arrays will hold the output.
		dst_coef_arrays = jtransform_adjust_parameters(&srcinfo, &dstinfo, src_coef_arrays, &config->trans);

		if( src_arg_type == BY_NAME ) {

			// Close input file, if we opened it.
			// Note: we assume that jpeg_read_coefficients consumed all input
			// until JPEG_REACHED_EOI, and that jpeg_finish_decompress will
			// only consume more while (! cinfo->inputctl->eoi_reached).
			// We cannot call jpeg_finish_decompress here since we still need the
			// virtual arrays allocated from the source object for processing.
			fclose(fp);
			fp = NULL;
		}

		if (config->progressive) {
			jpeg_simple_progression(&dstinfo);
		}
		
		if( dst_arg_type == BY_NAME ) {

			/* Open the output file. */
			if ((fp = fopen(dst_name, WRITE_BINARY)) == NULL) {
				my_croak("can't open `%s' for writing",dst);
			}
			// Specify data destination for compression
			jpeg_stdio_dest(&dstinfo, fp);

		} else {
			jpeg_mem_dest( &dstinfo, &buf, &buf_size );
		}

		//#ifdef C_MULTISCAN_FILES_SUPPORTED
		//if ((key = hv_fetch(conf, "scans", 5, 0))) {
		//	if (SvPOK(*key)) {
		//		char *scansarg = SvPV_nolen(*key);
		//		if (! read_scan_script(&dstinfo, scansarg)) {
		//			my_croak("Can't read scans script `%s'",scansarg);
		//		}
		//	}
		//}
		//#endif

		// Start compressor (note no image data is actually written here)
		jpeg_write_coefficients(&dstinfo, dst_coef_arrays);

		// Copy to the output file any extra markers that we want to preserve
		jcopy_markers_execute(&srcinfo, &dstinfo, config->copy);

		jtransform_execute_transformation(&srcinfo, &dstinfo, src_coef_arrays, &config->trans);

		// Finish compression and release memory
		jpeg_finish_compress(&dstinfo);
		jpeg_destroy_compress(&dstinfo);
		
		(void) jpeg_finish_decompress(&srcinfo);
		jpeg_destroy_decompress(&srcinfo);
		if( dst_arg_type == BY_NAME ) {	
			fclose(fp);
		} else {
			sv_setpvn( dst, buf, buf_size );
			free(buf);
		}
		if( jsrcerr.pub.num_warnings + jdsterr.pub.num_warnings ) {
			warn("Compression/decompression have warings");
		}

}

#define AUTOTRAN_USAGE "Usage: jpeg*tran(src, [ dst, [ conf ] ])\n"
#define JPEGTRAN_ARGS(src,dst,conf) \
		char * src_name; \
		char * dst_name; \
		int src_arg_type;\
		int dst_arg_type;\
		SV * dst; \
		HV * conf; \
		if (items > 1 && SvOK(ST(1))) { \
			if (SvROK(ST(1)) ) { \
				if( SvTYPE(SvRV(ST(1))) == SVt_PV ){ \
					dst = (SV *)SvRV(ST(1)); \
					dst_arg_type = BY_REF;\
				} else { \
					croak("Bad destination reference. " AUTOTRAN_USAGE); \
				} \
			} else { \
				dst_name = SvPV_nolen( ST(1) ); \
				dst_arg_type = BY_NAME;\
			} \
			src_arg_type = BY_REF;\
			if (! SvROK(src)) { \
				src_name = (char * ) SvPV_nolen( src ); \
				src_arg_type = BY_NAME;\
			} else { \
				if( SvTYPE(SvRV(src)) != SVt_PV ) { \
					croak("Bad source reference. " AUTOTRAN_USAGE); \
				} \
			} \
		} else { \
			if (SvROK(src)) { \
				croak("Second reference or filename is missing. " AUTOTRAN_USAGE); \
			} else { \
				dst = src; \
			} \
		} \
		if (items > 2) { \
			if (SvROK(ST(2)) && SvTYPE(SvRV(ST(2))) == SVt_PVHV ) { \
				conf = (HV *)SvRV(ST(2)); \
			} else { \
				croak ("Bad argumernts. " AUTOTRAN_USAGE); \
			} \
		} else { \
			conf = newHV(); \
			SV *keep = sv_2mortal(newRV_noinc((SV *)conf)); \
		} \
		if (items > 3) { \
			croak("Too many arguments. " AUTOTRAN_USAGE); \
		}
	
MODULE = Image::JpegTran		PACKAGE = Image::JpegTran

int
_jpegautotran(src,...)
		SV *src;
	PROTOTYPE: $;$$
	CODE:
		JPEGTRAN_ARGS(src,dst,conf);
		int k,l;
		jpegtran_config     config;
		ExifLoader  * loader;
		ExifData    * data;
		ExifContent * content;
		ExifEntry   * entry;
		
		parse_config( conf, &config );
		
		loader = exif_loader_new();
		if( src_arg_type == BY_REF ) { // this is binary by itself 
			SV * str_r;
			str_r = SvRV( src );
			STRLEN len;
			char *srcbin  = SvPVbyte( str_r,len );
			exif_loader_write( loader, srcbin, len );
		} else {
			exif_loader_write_file(loader, src_name);
		}
		data = exif_loader_get_data(loader);
		exif_loader_unref(loader);
		if (!data) XSRETURN_UNDEF;
		
		char val[255];
		char orientation = 0;
		ExifByteOrder o = exif_data_get_byte_order(data);
		for (k = 0; k < EXIF_IFD_COUNT; k++) {
			content = data->ifd[k];
			if (!content) continue;
			for (l = 0; l < content->count; l++) {
				entry = content->entries[l];
				if (entry->tag == 0x112) {
					orientation = 
						entry->format == EXIF_FORMAT_BYTE || entry->format == EXIF_FORMAT_SBYTE
							? (unsigned char) entry->data[0] :
						entry->format == EXIF_FORMAT_SHORT || entry->format == EXIF_FORMAT_SSHORT
							? exif_get_short(entry->data,o) :
						entry->format == EXIF_FORMAT_LONG || entry->format == EXIF_FORMAT_SLONG
							? exif_get_long(entry->data,o) :
							0;
					//warn("found orientation = %d (%s)\n", orientation, exif_entry_get_value(entry,val,255));
					k = EXIF_IFD_COUNT;
					break;
				}
			}
		}
		if (orientation > 1 && orientation < 9) { // don't touch orientation = 1 (normal)
			//warn("convert %d using %s\n", orientation, JXFORM_NAME[autocorrect[orientation]]);
			config.trans.transform = autocorrect[orientation];
			jpegtran_execute(src, dst, src_name, dst_name, src_arg_type, dst_arg_type, &config);
			ST(0) = sv_2mortal(newSViv(orientation));
			XSRETURN(1);
		} else {
			config.trans.transform = 0;
			jpegtran_execute(src, dst, src_name, dst_name, src_arg_type, dst_arg_type, &config);
			ST(0) = sv_2mortal(newSViv(orientation));
			XSRETURN(1);
			//XSRETURN_UNDEF;
		}
		

void
_jpegtran(src,...)
		SV *src;
	PROTOTYPE: $;$$
	CODE:
		JPEGTRAN_ARGS(src,dst,conf);
		
		jpegtran_config config;
		parse_config( conf, &config );
		jpegtran_execute(src, dst, src_name, dst_name, src_arg_type, dst_arg_type, &config);
		return;
