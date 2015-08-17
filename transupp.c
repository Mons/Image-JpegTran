/*
 * transupp.c
 *
 * Copyright (C) 1997-2009, Thomas G. Lane, Guido Vollbeding.
 * This file is part of the Independent JPEG Group's software.
 * For conditions of distribution and use, see the accompanying README file.
 *
 * This file contains image transformation routines and other utility code
 * used by the jpegtran sample application.  These are NOT part of the core
 * JPEG library.  But we keep these routines separate from jpegtran.c to
 * ease the task of maintaining jpegtran-like programs that have other user
 * interfaces.
 */

/* Although this file really shouldn't have access to the library internals,
 * it's helpful to let it call jround_up() and jcopy_block_row().
 */
#define JPEG_INTERNALS

#include "jinclude.h"
#include <jpeglib.h>
#include "transupp.h"		/* My own external interface */
#include <ctype.h>		/* to declare isdigit() */

#include <libexif/exif-loader.h>
#include <libexif/exif-data.h>
#include <libexif/exif-log.h>
#include <libexif/exif-mem.h>

#define match(x,y) (!memcmp(x,y,sizeof(y)))

/*
 * Lossless image transformation routines.  These routines work on DCT
 * coefficient arrays and thus do not require any lossy decompression
 * or recompression of the image.
 * Thanks to Guido Vollbeding for the initial design and code of this feature,
 * and to Ben Jackson for introducing the cropping feature.
 *
 * Horizontal flipping is done in-place, using a single top-to-bottom
 * pass through the virtual source array.  It will thus be much the
 * fastest option for images larger than main memory.
 *
 * The other routines require a set of destination virtual arrays, so they
 * need twice as much memory as jpegtran normally does.  The destination
 * arrays are always written in normal scan order (top to bottom) because
 * the virtual array manager expects this.  The source arrays will be scanned
 * in the corresponding order, which means multiple passes through the source
 * arrays for most of the transforms.  That could result in much thrashing
 * if the image is larger than main memory.
 *
 * If cropping or trimming is involved, the destination arrays may be smaller
 * than the source arrays.  Note it is not possible to do horizontal flip
 * in-place when a nonzero Y crop offset is specified, since we'd have to move
 * data from one block row to another but the virtual array manager doesn't
 * guarantee we can touch more than one row at a time.  So in that case,
 * we have to use a separate destination array.
 *
 * Some notes about the operating environment of the individual transform
 * routines:
 * 1. Both the source and destination virtual arrays are allocated from the
 *    source JPEG object, and therefore should be manipulated by calling the
 *    source's memory manager.
 * 2. The destination's component count should be used.  It may be smaller
 *    than the source's when forcing to grayscale.
 * 3. Likewise the destination's sampling factors should be used.  When
 *    forcing to grayscale the destination's sampling factors will be all 1,
 *    and we may as well take that as the effective iMCU size.
 * 4. When "trim" is in effect, the destination's dimensions will be the
 *    trimmed values but the source's will be untrimmed.
 * 5. When "crop" is in effect, the destination's dimensions will be the
 *    cropped values but the source's will be uncropped.  Each transform
 *    routine is responsible for picking up source data starting at the
 *    correct X and Y offset for the crop region.  (The X and Y offsets
 *    passed to the transform routines are measured in iMCU blocks of the
 *    destination.)
 * 6. All the routines assume that the source and destination buffers are
 *    padded out to a full iMCU boundary.  This is true, although for the
 *    source buffer it is an undocumented property of jdcoefct.c.
 */


LOCAL(void)
do_crop (j_decompress_ptr srcinfo, j_compress_ptr dstinfo,
	 JDIMENSION x_crop_offset, JDIMENSION y_crop_offset,
	 jvirt_barray_ptr *src_coef_arrays,
	 jvirt_barray_ptr *dst_coef_arrays)
/* Crop.  This is only used when no rotate/flip is requested with the crop. */
{
  JDIMENSION dst_blk_y, x_crop_blocks, y_crop_blocks;
  int ci, offset_y;
  JBLOCKARRAY src_buffer, dst_buffer;
  jpeg_component_info *compptr;

  /* We simply have to copy the right amount of data (the destination's
   * image size) starting at the given X and Y offsets in the source.
   */
  for (ci = 0; ci < dstinfo->num_components; ci++) {
    compptr = dstinfo->comp_info + ci;
    x_crop_blocks = x_crop_offset * compptr->h_samp_factor;
    y_crop_blocks = y_crop_offset * compptr->v_samp_factor;
    for (dst_blk_y = 0; dst_blk_y < compptr->height_in_blocks;
	 dst_blk_y += compptr->v_samp_factor) {
      dst_buffer = (*srcinfo->mem->access_virt_barray)
	((j_common_ptr) srcinfo, dst_coef_arrays[ci], dst_blk_y,
	 (JDIMENSION) compptr->v_samp_factor, TRUE);
      src_buffer = (*srcinfo->mem->access_virt_barray)
	((j_common_ptr) srcinfo, src_coef_arrays[ci],
	 dst_blk_y + y_crop_blocks,
	 (JDIMENSION) compptr->v_samp_factor, FALSE);
      for (offset_y = 0; offset_y < compptr->v_samp_factor; offset_y++) {
	jcopy_block_row(src_buffer[offset_y] + x_crop_blocks,
			dst_buffer[offset_y],
			compptr->width_in_blocks);
      }
    }
  }
}


LOCAL(void)
do_flip_h_no_crop (j_decompress_ptr srcinfo, j_compress_ptr dstinfo,
		   JDIMENSION x_crop_offset,
		   jvirt_barray_ptr *src_coef_arrays)
/* Horizontal flip; done in-place, so no separate dest array is required.
 * NB: this only works when y_crop_offset is zero.
 */
{
  JDIMENSION MCU_cols, comp_width, blk_x, blk_y, x_crop_blocks;
  int ci, k, offset_y;
  JBLOCKARRAY buffer;
  JCOEFPTR ptr1, ptr2;
  JCOEF temp1, temp2;
  jpeg_component_info *compptr;

  /* Horizontal mirroring of DCT blocks is accomplished by swapping
   * pairs of blocks in-place.  Within a DCT block, we perform horizontal
   * mirroring by changing the signs of odd-numbered columns.
   * Partial iMCUs at the right edge are left untouched.
   */
  MCU_cols = srcinfo->output_width /
    (dstinfo->max_h_samp_factor * dstinfo->min_DCT_h_scaled_size);

  for (ci = 0; ci < dstinfo->num_components; ci++) {
    compptr = dstinfo->comp_info + ci;
    comp_width = MCU_cols * compptr->h_samp_factor;
    x_crop_blocks = x_crop_offset * compptr->h_samp_factor;
    for (blk_y = 0; blk_y < compptr->height_in_blocks;
	 blk_y += compptr->v_samp_factor) {
      buffer = (*srcinfo->mem->access_virt_barray)
	((j_common_ptr) srcinfo, src_coef_arrays[ci], blk_y,
	 (JDIMENSION) compptr->v_samp_factor, TRUE);
      for (offset_y = 0; offset_y < compptr->v_samp_factor; offset_y++) {
	/* Do the mirroring */
	for (blk_x = 0; blk_x * 2 < comp_width; blk_x++) {
	  ptr1 = buffer[offset_y][blk_x];
	  ptr2 = buffer[offset_y][comp_width - blk_x - 1];
	  /* this unrolled loop doesn't need to know which row it's on... */
	  for (k = 0; k < DCTSIZE2; k += 2) {
	    temp1 = *ptr1;	/* swap even column */
	    temp2 = *ptr2;
	    *ptr1++ = temp2;
	    *ptr2++ = temp1;
	    temp1 = *ptr1;	/* swap odd column with sign change */
	    temp2 = *ptr2;
	    *ptr1++ = -temp2;
	    *ptr2++ = -temp1;
	  }
	}
	if (x_crop_blocks > 0) {
	  /* Now left-justify the portion of the data to be kept.
	   * We can't use a single jcopy_block_row() call because that routine
	   * depends on memcpy(), whose behavior is unspecified for overlapping
	   * source and destination areas.  Sigh.
	   */
	  for (blk_x = 0; blk_x < compptr->width_in_blocks; blk_x++) {
	    jcopy_block_row(buffer[offset_y] + blk_x + x_crop_blocks,
			    buffer[offset_y] + blk_x,
			    (JDIMENSION) 1);
	  }
	}
      }
    }
  }
}


LOCAL(void)
do_flip_h (j_decompress_ptr srcinfo, j_compress_ptr dstinfo,
	   JDIMENSION x_crop_offset, JDIMENSION y_crop_offset,
	   jvirt_barray_ptr *src_coef_arrays,
	   jvirt_barray_ptr *dst_coef_arrays)
/* Horizontal flip in general cropping case */
{
  JDIMENSION MCU_cols, comp_width, dst_blk_x, dst_blk_y;
  JDIMENSION x_crop_blocks, y_crop_blocks;
  int ci, k, offset_y;
  JBLOCKARRAY src_buffer, dst_buffer;
  JBLOCKROW src_row_ptr, dst_row_ptr;
  JCOEFPTR src_ptr, dst_ptr;
  jpeg_component_info *compptr;

  /* Here we must output into a separate array because we can't touch
   * different rows of a single virtual array simultaneously.  Otherwise,
   * this is essentially the same as the routine above.
   */
  MCU_cols = srcinfo->output_width /
    (dstinfo->max_h_samp_factor * dstinfo->min_DCT_h_scaled_size);

  for (ci = 0; ci < dstinfo->num_components; ci++) {
    compptr = dstinfo->comp_info + ci;
    comp_width = MCU_cols * compptr->h_samp_factor;
    x_crop_blocks = x_crop_offset * compptr->h_samp_factor;
    y_crop_blocks = y_crop_offset * compptr->v_samp_factor;
    for (dst_blk_y = 0; dst_blk_y < compptr->height_in_blocks;
	 dst_blk_y += compptr->v_samp_factor) {
      dst_buffer = (*srcinfo->mem->access_virt_barray)
	((j_common_ptr) srcinfo, dst_coef_arrays[ci], dst_blk_y,
	 (JDIMENSION) compptr->v_samp_factor, TRUE);
      src_buffer = (*srcinfo->mem->access_virt_barray)
	((j_common_ptr) srcinfo, src_coef_arrays[ci],
	 dst_blk_y + y_crop_blocks,
	 (JDIMENSION) compptr->v_samp_factor, FALSE);
      for (offset_y = 0; offset_y < compptr->v_samp_factor; offset_y++) {
	dst_row_ptr = dst_buffer[offset_y];
	src_row_ptr = src_buffer[offset_y];
	for (dst_blk_x = 0; dst_blk_x < compptr->width_in_blocks; dst_blk_x++) {
	  if (x_crop_blocks + dst_blk_x < comp_width) {
	    /* Do the mirrorable blocks */
	    dst_ptr = dst_row_ptr[dst_blk_x];
	    src_ptr = src_row_ptr[comp_width - x_crop_blocks - dst_blk_x - 1];
	    /* this unrolled loop doesn't need to know which row it's on... */
	    for (k = 0; k < DCTSIZE2; k += 2) {
	      *dst_ptr++ = *src_ptr++;	 /* copy even column */
	      *dst_ptr++ = - *src_ptr++; /* copy odd column with sign change */
	    }
	  } else {
	    /* Copy last partial block(s) verbatim */
	    jcopy_block_row(src_row_ptr + dst_blk_x + x_crop_blocks,
			    dst_row_ptr + dst_blk_x,
			    (JDIMENSION) 1);
	  }
	}
      }
    }
  }
}


LOCAL(void)
do_flip_v (j_decompress_ptr srcinfo, j_compress_ptr dstinfo,
	   JDIMENSION x_crop_offset, JDIMENSION y_crop_offset,
	   jvirt_barray_ptr *src_coef_arrays,
	   jvirt_barray_ptr *dst_coef_arrays)
/* Vertical flip */
{
  JDIMENSION MCU_rows, comp_height, dst_blk_x, dst_blk_y;
  JDIMENSION x_crop_blocks, y_crop_blocks;
  int ci, i, j, offset_y;
  JBLOCKARRAY src_buffer, dst_buffer;
  JBLOCKROW src_row_ptr, dst_row_ptr;
  JCOEFPTR src_ptr, dst_ptr;
  jpeg_component_info *compptr;

  /* We output into a separate array because we can't touch different
   * rows of the source virtual array simultaneously.  Otherwise, this
   * is a pretty straightforward analog of horizontal flip.
   * Within a DCT block, vertical mirroring is done by changing the signs
   * of odd-numbered rows.
   * Partial iMCUs at the bottom edge are copied verbatim.
   */
  MCU_rows = srcinfo->output_height /
    (dstinfo->max_v_samp_factor * dstinfo->min_DCT_v_scaled_size);

  for (ci = 0; ci < dstinfo->num_components; ci++) {
    compptr = dstinfo->comp_info + ci;
    comp_height = MCU_rows * compptr->v_samp_factor;
    x_crop_blocks = x_crop_offset * compptr->h_samp_factor;
    y_crop_blocks = y_crop_offset * compptr->v_samp_factor;
    for (dst_blk_y = 0; dst_blk_y < compptr->height_in_blocks;
	 dst_blk_y += compptr->v_samp_factor) {
      dst_buffer = (*srcinfo->mem->access_virt_barray)
	((j_common_ptr) srcinfo, dst_coef_arrays[ci], dst_blk_y,
	 (JDIMENSION) compptr->v_samp_factor, TRUE);
      if (y_crop_blocks + dst_blk_y < comp_height) {
	/* Row is within the mirrorable area. */
	src_buffer = (*srcinfo->mem->access_virt_barray)
	  ((j_common_ptr) srcinfo, src_coef_arrays[ci],
	   comp_height - y_crop_blocks - dst_blk_y -
	   (JDIMENSION) compptr->v_samp_factor,
	   (JDIMENSION) compptr->v_samp_factor, FALSE);
      } else {
	/* Bottom-edge blocks will be copied verbatim. */
	src_buffer = (*srcinfo->mem->access_virt_barray)
	  ((j_common_ptr) srcinfo, src_coef_arrays[ci],
	   dst_blk_y + y_crop_blocks,
	   (JDIMENSION) compptr->v_samp_factor, FALSE);
      }
      for (offset_y = 0; offset_y < compptr->v_samp_factor; offset_y++) {
	if (y_crop_blocks + dst_blk_y < comp_height) {
	  /* Row is within the mirrorable area. */
	  dst_row_ptr = dst_buffer[offset_y];
	  src_row_ptr = src_buffer[compptr->v_samp_factor - offset_y - 1];
	  src_row_ptr += x_crop_blocks;
	  for (dst_blk_x = 0; dst_blk_x < compptr->width_in_blocks;
	       dst_blk_x++) {
	    dst_ptr = dst_row_ptr[dst_blk_x];
	    src_ptr = src_row_ptr[dst_blk_x];
	    for (i = 0; i < DCTSIZE; i += 2) {
	      /* copy even row */
	      for (j = 0; j < DCTSIZE; j++)
		*dst_ptr++ = *src_ptr++;
	      /* copy odd row with sign change */
	      for (j = 0; j < DCTSIZE; j++)
		*dst_ptr++ = - *src_ptr++;
	    }
	  }
	} else {
	  /* Just copy row verbatim. */
	  jcopy_block_row(src_buffer[offset_y] + x_crop_blocks,
			  dst_buffer[offset_y],
			  compptr->width_in_blocks);
	}
      }
    }
  }
}


LOCAL(void)
do_transpose (j_decompress_ptr srcinfo, j_compress_ptr dstinfo,
	      JDIMENSION x_crop_offset, JDIMENSION y_crop_offset,
	      jvirt_barray_ptr *src_coef_arrays,
	      jvirt_barray_ptr *dst_coef_arrays)
/* Transpose source into destination */
{
  JDIMENSION dst_blk_x, dst_blk_y, x_crop_blocks, y_crop_blocks;
  int ci, i, j, offset_x, offset_y;
  JBLOCKARRAY src_buffer, dst_buffer;
  JCOEFPTR src_ptr, dst_ptr;
  jpeg_component_info *compptr;

  /* Transposing pixels within a block just requires transposing the
   * DCT coefficients.
   * Partial iMCUs at the edges require no special treatment; we simply
   * process all the available DCT blocks for every component.
   */
  for (ci = 0; ci < dstinfo->num_components; ci++) {
    compptr = dstinfo->comp_info + ci;
    x_crop_blocks = x_crop_offset * compptr->h_samp_factor;
    y_crop_blocks = y_crop_offset * compptr->v_samp_factor;
    for (dst_blk_y = 0; dst_blk_y < compptr->height_in_blocks;
	 dst_blk_y += compptr->v_samp_factor) {
      dst_buffer = (*srcinfo->mem->access_virt_barray)
	((j_common_ptr) srcinfo, dst_coef_arrays[ci], dst_blk_y,
	 (JDIMENSION) compptr->v_samp_factor, TRUE);
      for (offset_y = 0; offset_y < compptr->v_samp_factor; offset_y++) {
	for (dst_blk_x = 0; dst_blk_x < compptr->width_in_blocks;
	     dst_blk_x += compptr->h_samp_factor) {
	  src_buffer = (*srcinfo->mem->access_virt_barray)
	    ((j_common_ptr) srcinfo, src_coef_arrays[ci],
	     dst_blk_x + x_crop_blocks,
	     (JDIMENSION) compptr->h_samp_factor, FALSE);
	  for (offset_x = 0; offset_x < compptr->h_samp_factor; offset_x++) {
	    dst_ptr = dst_buffer[offset_y][dst_blk_x + offset_x];
	    src_ptr = src_buffer[offset_x][dst_blk_y + offset_y + y_crop_blocks];
	    for (i = 0; i < DCTSIZE; i++)
	      for (j = 0; j < DCTSIZE; j++)
		dst_ptr[j*DCTSIZE+i] = src_ptr[i*DCTSIZE+j];
	  }
	}
      }
    }
  }
}


LOCAL(void)
do_rot_90 (j_decompress_ptr srcinfo, j_compress_ptr dstinfo,
	   JDIMENSION x_crop_offset, JDIMENSION y_crop_offset,
	   jvirt_barray_ptr *src_coef_arrays,
	   jvirt_barray_ptr *dst_coef_arrays)
/* 90 degree rotation is equivalent to
 *   1. Transposing the image;
 *   2. Horizontal mirroring.
 * These two steps are merged into a single processing routine.
 */
{
  JDIMENSION MCU_cols, comp_width, dst_blk_x, dst_blk_y;
  JDIMENSION x_crop_blocks, y_crop_blocks;
  int ci, i, j, offset_x, offset_y;
  JBLOCKARRAY src_buffer, dst_buffer;
  JCOEFPTR src_ptr, dst_ptr;
  jpeg_component_info *compptr;

  /* Because of the horizontal mirror step, we can't process partial iMCUs
   * at the (output) right edge properly.  They just get transposed and
   * not mirrored.
   */
  MCU_cols = srcinfo->output_height /
    (dstinfo->max_h_samp_factor * dstinfo->min_DCT_h_scaled_size);

  for (ci = 0; ci < dstinfo->num_components; ci++) {
    compptr = dstinfo->comp_info + ci;
    comp_width = MCU_cols * compptr->h_samp_factor;
    x_crop_blocks = x_crop_offset * compptr->h_samp_factor;
    y_crop_blocks = y_crop_offset * compptr->v_samp_factor;
    for (dst_blk_y = 0; dst_blk_y < compptr->height_in_blocks;
	 dst_blk_y += compptr->v_samp_factor) {
      dst_buffer = (*srcinfo->mem->access_virt_barray)
	((j_common_ptr) srcinfo, dst_coef_arrays[ci], dst_blk_y,
	 (JDIMENSION) compptr->v_samp_factor, TRUE);
      for (offset_y = 0; offset_y < compptr->v_samp_factor; offset_y++) {
	for (dst_blk_x = 0; dst_blk_x < compptr->width_in_blocks;
	     dst_blk_x += compptr->h_samp_factor) {
	  if (x_crop_blocks + dst_blk_x < comp_width) {
	    /* Block is within the mirrorable area. */
	    src_buffer = (*srcinfo->mem->access_virt_barray)
	      ((j_common_ptr) srcinfo, src_coef_arrays[ci],
	       comp_width - x_crop_blocks - dst_blk_x -
	       (JDIMENSION) compptr->h_samp_factor,
	       (JDIMENSION) compptr->h_samp_factor, FALSE);
	  } else {
	    /* Edge blocks are transposed but not mirrored. */
	    src_buffer = (*srcinfo->mem->access_virt_barray)
	      ((j_common_ptr) srcinfo, src_coef_arrays[ci],
	       dst_blk_x + x_crop_blocks,
	       (JDIMENSION) compptr->h_samp_factor, FALSE);
	  }
	  for (offset_x = 0; offset_x < compptr->h_samp_factor; offset_x++) {
	    dst_ptr = dst_buffer[offset_y][dst_blk_x + offset_x];
	    if (x_crop_blocks + dst_blk_x < comp_width) {
	      /* Block is within the mirrorable area. */
	      src_ptr = src_buffer[compptr->h_samp_factor - offset_x - 1]
		[dst_blk_y + offset_y + y_crop_blocks];
	      for (i = 0; i < DCTSIZE; i++) {
		for (j = 0; j < DCTSIZE; j++)
		  dst_ptr[j*DCTSIZE+i] = src_ptr[i*DCTSIZE+j];
		i++;
		for (j = 0; j < DCTSIZE; j++)
		  dst_ptr[j*DCTSIZE+i] = -src_ptr[i*DCTSIZE+j];
	      }
	    } else {
	      /* Edge blocks are transposed but not mirrored. */
	      src_ptr = src_buffer[offset_x]
		[dst_blk_y + offset_y + y_crop_blocks];
	      for (i = 0; i < DCTSIZE; i++)
		for (j = 0; j < DCTSIZE; j++)
		  dst_ptr[j*DCTSIZE+i] = src_ptr[i*DCTSIZE+j];
	    }
	  }
	}
      }
    }
  }
}


LOCAL(void)
do_rot_270 (j_decompress_ptr srcinfo, j_compress_ptr dstinfo,
	    JDIMENSION x_crop_offset, JDIMENSION y_crop_offset,
	    jvirt_barray_ptr *src_coef_arrays,
	    jvirt_barray_ptr *dst_coef_arrays)
/* 270 degree rotation is equivalent to
 *   1. Horizontal mirroring;
 *   2. Transposing the image.
 * These two steps are merged into a single processing routine.
 */
{
  JDIMENSION MCU_rows, comp_height, dst_blk_x, dst_blk_y;
  JDIMENSION x_crop_blocks, y_crop_blocks;
  int ci, i, j, offset_x, offset_y;
  JBLOCKARRAY src_buffer, dst_buffer;
  JCOEFPTR src_ptr, dst_ptr;
  jpeg_component_info *compptr;

  /* Because of the horizontal mirror step, we can't process partial iMCUs
   * at the (output) bottom edge properly.  They just get transposed and
   * not mirrored.
   */
  MCU_rows = srcinfo->output_width /
    (dstinfo->max_v_samp_factor * dstinfo->min_DCT_v_scaled_size);

  for (ci = 0; ci < dstinfo->num_components; ci++) {
    compptr = dstinfo->comp_info + ci;
    comp_height = MCU_rows * compptr->v_samp_factor;
    x_crop_blocks = x_crop_offset * compptr->h_samp_factor;
    y_crop_blocks = y_crop_offset * compptr->v_samp_factor;
    for (dst_blk_y = 0; dst_blk_y < compptr->height_in_blocks;
	 dst_blk_y += compptr->v_samp_factor) {
      dst_buffer = (*srcinfo->mem->access_virt_barray)
	((j_common_ptr) srcinfo, dst_coef_arrays[ci], dst_blk_y,
	 (JDIMENSION) compptr->v_samp_factor, TRUE);
      for (offset_y = 0; offset_y < compptr->v_samp_factor; offset_y++) {
	for (dst_blk_x = 0; dst_blk_x < compptr->width_in_blocks;
	     dst_blk_x += compptr->h_samp_factor) {
	  src_buffer = (*srcinfo->mem->access_virt_barray)
	    ((j_common_ptr) srcinfo, src_coef_arrays[ci],
	     dst_blk_x + x_crop_blocks,
	     (JDIMENSION) compptr->h_samp_factor, FALSE);
	  for (offset_x = 0; offset_x < compptr->h_samp_factor; offset_x++) {
	    dst_ptr = dst_buffer[offset_y][dst_blk_x + offset_x];
	    if (y_crop_blocks + dst_blk_y < comp_height) {
	      /* Block is within the mirrorable area. */
	      src_ptr = src_buffer[offset_x]
		[comp_height - y_crop_blocks - dst_blk_y - offset_y - 1];
	      for (i = 0; i < DCTSIZE; i++) {
		for (j = 0; j < DCTSIZE; j++) {
		  dst_ptr[j*DCTSIZE+i] = src_ptr[i*DCTSIZE+j];
		  j++;
		  dst_ptr[j*DCTSIZE+i] = -src_ptr[i*DCTSIZE+j];
		}
	      }
	    } else {
	      /* Edge blocks are transposed but not mirrored. */
	      src_ptr = src_buffer[offset_x]
		[dst_blk_y + offset_y + y_crop_blocks];
	      for (i = 0; i < DCTSIZE; i++)
		for (j = 0; j < DCTSIZE; j++)
		  dst_ptr[j*DCTSIZE+i] = src_ptr[i*DCTSIZE+j];
	    }
	  }
	}
      }
    }
  }
}


LOCAL(void)
do_rot_180 (j_decompress_ptr srcinfo, j_compress_ptr dstinfo,
	    JDIMENSION x_crop_offset, JDIMENSION y_crop_offset,
	    jvirt_barray_ptr *src_coef_arrays,
	    jvirt_barray_ptr *dst_coef_arrays)
/* 180 degree rotation is equivalent to
 *   1. Vertical mirroring;
 *   2. Horizontal mirroring.
 * These two steps are merged into a single processing routine.
 */
{
  JDIMENSION MCU_cols, MCU_rows, comp_width, comp_height, dst_blk_x, dst_blk_y;
  JDIMENSION x_crop_blocks, y_crop_blocks;
  int ci, i, j, offset_y;
  JBLOCKARRAY src_buffer, dst_buffer;
  JBLOCKROW src_row_ptr, dst_row_ptr;
  JCOEFPTR src_ptr, dst_ptr;
  jpeg_component_info *compptr;

  MCU_cols = srcinfo->output_width /
    (dstinfo->max_h_samp_factor * dstinfo->min_DCT_h_scaled_size);
  MCU_rows = srcinfo->output_height /
    (dstinfo->max_v_samp_factor * dstinfo->min_DCT_v_scaled_size);

  for (ci = 0; ci < dstinfo->num_components; ci++) {
    compptr = dstinfo->comp_info + ci;
    comp_width = MCU_cols * compptr->h_samp_factor;
    comp_height = MCU_rows * compptr->v_samp_factor;
    x_crop_blocks = x_crop_offset * compptr->h_samp_factor;
    y_crop_blocks = y_crop_offset * compptr->v_samp_factor;
    for (dst_blk_y = 0; dst_blk_y < compptr->height_in_blocks;
	 dst_blk_y += compptr->v_samp_factor) {
      dst_buffer = (*srcinfo->mem->access_virt_barray)
	((j_common_ptr) srcinfo, dst_coef_arrays[ci], dst_blk_y,
	 (JDIMENSION) compptr->v_samp_factor, TRUE);
      if (y_crop_blocks + dst_blk_y < comp_height) {
	/* Row is within the vertically mirrorable area. */
	src_buffer = (*srcinfo->mem->access_virt_barray)
	  ((j_common_ptr) srcinfo, src_coef_arrays[ci],
	   comp_height - y_crop_blocks - dst_blk_y -
	   (JDIMENSION) compptr->v_samp_factor,
	   (JDIMENSION) compptr->v_samp_factor, FALSE);
      } else {
	/* Bottom-edge rows are only mirrored horizontally. */
	src_buffer = (*srcinfo->mem->access_virt_barray)
	  ((j_common_ptr) srcinfo, src_coef_arrays[ci],
	   dst_blk_y + y_crop_blocks,
	   (JDIMENSION) compptr->v_samp_factor, FALSE);
      }
      for (offset_y = 0; offset_y < compptr->v_samp_factor; offset_y++) {
	dst_row_ptr = dst_buffer[offset_y];
	if (y_crop_blocks + dst_blk_y < comp_height) {
	  /* Row is within the mirrorable area. */
	  src_row_ptr = src_buffer[compptr->v_samp_factor - offset_y - 1];
	  for (dst_blk_x = 0; dst_blk_x < compptr->width_in_blocks; dst_blk_x++) {
	    dst_ptr = dst_row_ptr[dst_blk_x];
	    if (x_crop_blocks + dst_blk_x < comp_width) {
	      /* Process the blocks that can be mirrored both ways. */
	      src_ptr = src_row_ptr[comp_width - x_crop_blocks - dst_blk_x - 1];
	      for (i = 0; i < DCTSIZE; i += 2) {
		/* For even row, negate every odd column. */
		for (j = 0; j < DCTSIZE; j += 2) {
		  *dst_ptr++ = *src_ptr++;
		  *dst_ptr++ = - *src_ptr++;
		}
		/* For odd row, negate every even column. */
		for (j = 0; j < DCTSIZE; j += 2) {
		  *dst_ptr++ = - *src_ptr++;
		  *dst_ptr++ = *src_ptr++;
		}
	      }
	    } else {
	      /* Any remaining right-edge blocks are only mirrored vertically. */
	      src_ptr = src_row_ptr[x_crop_blocks + dst_blk_x];
	      for (i = 0; i < DCTSIZE; i += 2) {
		for (j = 0; j < DCTSIZE; j++)
		  *dst_ptr++ = *src_ptr++;
		for (j = 0; j < DCTSIZE; j++)
		  *dst_ptr++ = - *src_ptr++;
	      }
	    }
	  }
	} else {
	  /* Remaining rows are just mirrored horizontally. */
	  src_row_ptr = src_buffer[offset_y];
	  for (dst_blk_x = 0; dst_blk_x < compptr->width_in_blocks; dst_blk_x++) {
	    if (x_crop_blocks + dst_blk_x < comp_width) {
	      /* Process the blocks that can be mirrored. */
	      dst_ptr = dst_row_ptr[dst_blk_x];
	      src_ptr = src_row_ptr[comp_width - x_crop_blocks - dst_blk_x - 1];
	      for (i = 0; i < DCTSIZE2; i += 2) {
		*dst_ptr++ = *src_ptr++;
		*dst_ptr++ = - *src_ptr++;
	      }
	    } else {
	      /* Any remaining right-edge blocks are only copied. */
	      jcopy_block_row(src_row_ptr + dst_blk_x + x_crop_blocks,
			      dst_row_ptr + dst_blk_x,
			      (JDIMENSION) 1);
	    }
	  }
	}
      }
    }
  }
}


LOCAL(void)
do_transverse (j_decompress_ptr srcinfo, j_compress_ptr dstinfo,
	       JDIMENSION x_crop_offset, JDIMENSION y_crop_offset,
	       jvirt_barray_ptr *src_coef_arrays,
	       jvirt_barray_ptr *dst_coef_arrays)
/* Transverse transpose is equivalent to
 *   1. 180 degree rotation;
 *   2. Transposition;
 * or
 *   1. Horizontal mirroring;
 *   2. Transposition;
 *   3. Horizontal mirroring.
 * These steps are merged into a single processing routine.
 */
{
  JDIMENSION MCU_cols, MCU_rows, comp_width, comp_height, dst_blk_x, dst_blk_y;
  JDIMENSION x_crop_blocks, y_crop_blocks;
  int ci, i, j, offset_x, offset_y;
  JBLOCKARRAY src_buffer, dst_buffer;
  JCOEFPTR src_ptr, dst_ptr;
  jpeg_component_info *compptr;

  MCU_cols = srcinfo->output_height /
    (dstinfo->max_h_samp_factor * dstinfo->min_DCT_h_scaled_size);
  MCU_rows = srcinfo->output_width /
    (dstinfo->max_v_samp_factor * dstinfo->min_DCT_v_scaled_size);

  for (ci = 0; ci < dstinfo->num_components; ci++) {
    compptr = dstinfo->comp_info + ci;
    comp_width = MCU_cols * compptr->h_samp_factor;
    comp_height = MCU_rows * compptr->v_samp_factor;
    x_crop_blocks = x_crop_offset * compptr->h_samp_factor;
    y_crop_blocks = y_crop_offset * compptr->v_samp_factor;
    for (dst_blk_y = 0; dst_blk_y < compptr->height_in_blocks;
	 dst_blk_y += compptr->v_samp_factor) {
      dst_buffer = (*srcinfo->mem->access_virt_barray)
	((j_common_ptr) srcinfo, dst_coef_arrays[ci], dst_blk_y,
	 (JDIMENSION) compptr->v_samp_factor, TRUE);
      for (offset_y = 0; offset_y < compptr->v_samp_factor; offset_y++) {
	for (dst_blk_x = 0; dst_blk_x < compptr->width_in_blocks;
	     dst_blk_x += compptr->h_samp_factor) {
	  if (x_crop_blocks + dst_blk_x < comp_width) {
	    /* Block is within the mirrorable area. */
	    src_buffer = (*srcinfo->mem->access_virt_barray)
	      ((j_common_ptr) srcinfo, src_coef_arrays[ci],
	       comp_width - x_crop_blocks - dst_blk_x -
	       (JDIMENSION) compptr->h_samp_factor,
	       (JDIMENSION) compptr->h_samp_factor, FALSE);
	  } else {
	    src_buffer = (*srcinfo->mem->access_virt_barray)
	      ((j_common_ptr) srcinfo, src_coef_arrays[ci],
	       dst_blk_x + x_crop_blocks,
	       (JDIMENSION) compptr->h_samp_factor, FALSE);
	  }
	  for (offset_x = 0; offset_x < compptr->h_samp_factor; offset_x++) {
	    dst_ptr = dst_buffer[offset_y][dst_blk_x + offset_x];
	    if (y_crop_blocks + dst_blk_y < comp_height) {
	      if (x_crop_blocks + dst_blk_x < comp_width) {
		/* Block is within the mirrorable area. */
		src_ptr = src_buffer[compptr->h_samp_factor - offset_x - 1]
		  [comp_height - y_crop_blocks - dst_blk_y - offset_y - 1];
		for (i = 0; i < DCTSIZE; i++) {
		  for (j = 0; j < DCTSIZE; j++) {
		    dst_ptr[j*DCTSIZE+i] = src_ptr[i*DCTSIZE+j];
		    j++;
		    dst_ptr[j*DCTSIZE+i] = -src_ptr[i*DCTSIZE+j];
		  }
		  i++;
		  for (j = 0; j < DCTSIZE; j++) {
		    dst_ptr[j*DCTSIZE+i] = -src_ptr[i*DCTSIZE+j];
		    j++;
		    dst_ptr[j*DCTSIZE+i] = src_ptr[i*DCTSIZE+j];
		  }
		}
	      } else {
		/* Right-edge blocks are mirrored in y only */
		src_ptr = src_buffer[offset_x]
		  [comp_height - y_crop_blocks - dst_blk_y - offset_y - 1];
		for (i = 0; i < DCTSIZE; i++) {
		  for (j = 0; j < DCTSIZE; j++) {
		    dst_ptr[j*DCTSIZE+i] = src_ptr[i*DCTSIZE+j];
		    j++;
		    dst_ptr[j*DCTSIZE+i] = -src_ptr[i*DCTSIZE+j];
		  }
		}
	      }
	    } else {
	      if (x_crop_blocks + dst_blk_x < comp_width) {
		/* Bottom-edge blocks are mirrored in x only */
		src_ptr = src_buffer[compptr->h_samp_factor - offset_x - 1]
		  [dst_blk_y + offset_y + y_crop_blocks];
		for (i = 0; i < DCTSIZE; i++) {
		  for (j = 0; j < DCTSIZE; j++)
		    dst_ptr[j*DCTSIZE+i] = src_ptr[i*DCTSIZE+j];
		  i++;
		  for (j = 0; j < DCTSIZE; j++)
		    dst_ptr[j*DCTSIZE+i] = -src_ptr[i*DCTSIZE+j];
		}
	      } else {
		/* At lower right corner, just transpose, no mirroring */
		src_ptr = src_buffer[offset_x]
		  [dst_blk_y + offset_y + y_crop_blocks];
		for (i = 0; i < DCTSIZE; i++)
		  for (j = 0; j < DCTSIZE; j++)
		    dst_ptr[j*DCTSIZE+i] = src_ptr[i*DCTSIZE+j];
	      }
	    }
	  }
	}
      }
    }
  }
}


/* Parse an unsigned integer: subroutine for jtransform_parse_crop_spec.
 * Returns TRUE if valid integer found, FALSE if not.
 * *strptr is advanced over the digit string, and *result is set to its value.
 */

LOCAL(boolean)
jt_read_integer (const char ** strptr, JDIMENSION * result)
{
  const char * ptr = *strptr;
  JDIMENSION val = 0;

  for (; isdigit(*ptr); ptr++) {
    val = val * 10 + (JDIMENSION) (*ptr - '0');
  }
  *result = val;
  if (ptr == *strptr)
    return FALSE;		/* oops, no digits */
  *strptr = ptr;
  return TRUE;
}


/* Parse a crop specification (written in X11 geometry style).
 * The routine returns TRUE if the spec string is valid, FALSE if not.
 *
 * The crop spec string should have the format
 *	<width>x<height>{+-}<xoffset>{+-}<yoffset>
 * where width, height, xoffset, and yoffset are unsigned integers.
 * Each of the elements can be omitted to indicate a default value.
 * (A weakness of this style is that it is not possible to omit xoffset
 * while specifying yoffset, since they look alike.)
 *
 * This code is loosely based on XParseGeometry from the X11 distribution.
 */

GLOBAL(boolean)
jtransform_parse_crop_spec (jpeg_transform_info *info, const char *spec)
{
  info->crop = FALSE;
  info->crop_width_set = JCROP_UNSET;
  info->crop_height_set = JCROP_UNSET;
  info->crop_xoffset_set = JCROP_UNSET;
  info->crop_yoffset_set = JCROP_UNSET;

  if (isdigit(*spec)) {
    /* fetch width */
    if (! jt_read_integer(&spec, &info->crop_width))
      return FALSE;
    info->crop_width_set = JCROP_POS;
  }
  if (*spec == 'x' || *spec == 'X') {	
    /* fetch height */
    spec++;
    if (! jt_read_integer(&spec, &info->crop_height))
      return FALSE;
    info->crop_height_set = JCROP_POS;
  }
  if (*spec == '+' || *spec == '-') {
    /* fetch xoffset */
    info->crop_xoffset_set = (*spec == '-') ? JCROP_NEG : JCROP_POS;
    spec++;
    if (! jt_read_integer(&spec, &info->crop_xoffset))
      return FALSE;
  }
  if (*spec == '+' || *spec == '-') {
    /* fetch yoffset */
    info->crop_yoffset_set = (*spec == '-') ? JCROP_NEG : JCROP_POS;
    spec++;
    if (! jt_read_integer(&spec, &info->crop_yoffset))
      return FALSE;
  }
  /* We had better have gotten to the end of the string. */
  if (*spec != '\0')
    return FALSE;
  info->crop = TRUE;
  return TRUE;
}


/* Trim off any partial iMCUs on the indicated destination edge */

LOCAL(void)
trim_right_edge (jpeg_transform_info *info, JDIMENSION full_width)
{
  JDIMENSION MCU_cols;

  MCU_cols = info->output_width / info->iMCU_sample_width;
  if (MCU_cols > 0 && info->x_crop_offset + MCU_cols ==
      full_width / info->iMCU_sample_width)
    info->output_width = MCU_cols * info->iMCU_sample_width;
}

LOCAL(void)
trim_bottom_edge (jpeg_transform_info *info, JDIMENSION full_height)
{
  JDIMENSION MCU_rows;

  MCU_rows = info->output_height / info->iMCU_sample_height;
  if (MCU_rows > 0 && info->y_crop_offset + MCU_rows ==
      full_height / info->iMCU_sample_height)
    info->output_height = MCU_rows * info->iMCU_sample_height;
}


/* Request any required workspace.
 *
 * This routine figures out the size that the output image will be
 * (which implies that all the transform parameters must be set before
 * it is called).
 *
 * We allocate the workspace virtual arrays from the source decompression
 * object, so that all the arrays (both the original data and the workspace)
 * will be taken into account while making memory management decisions.
 * Hence, this routine must be called after jpeg_read_header (which reads
 * the image dimensions) and before jpeg_read_coefficients (which realizes
 * the source's virtual arrays).
 *
 * This function returns FALSE right away if -perfect is given
 * and transformation is not perfect.  Otherwise returns TRUE.
 */

GLOBAL(boolean)
jtransform_request_workspace (j_decompress_ptr srcinfo,
			      jpeg_transform_info *info)
{
  jvirt_barray_ptr *coef_arrays;
  boolean need_workspace, transpose_it;
  jpeg_component_info *compptr;
  JDIMENSION xoffset, yoffset;
  JDIMENSION width_in_iMCUs, height_in_iMCUs;
  JDIMENSION width_in_blocks, height_in_blocks;
  int ci, h_samp_factor, v_samp_factor;

  /* Determine number of components in output image */
  if (info->force_grayscale &&
      srcinfo->jpeg_color_space == JCS_YCbCr &&
      srcinfo->num_components == 3)
    /* We'll only process the first component */
    info->num_components = 1;
  else
    /* Process all the components */
    info->num_components = srcinfo->num_components;

  /* Compute output image dimensions and related values. */
  jpeg_core_output_dimensions(srcinfo);

  /* Return right away if -perfect is given and transformation is not perfect.
   */
  if (info->perfect) {
    if (info->num_components == 1) {
      if (!jtransform_perfect_transform(srcinfo->output_width,
	  srcinfo->output_height,
	  srcinfo->min_DCT_h_scaled_size,
	  srcinfo->min_DCT_v_scaled_size,
	  info->transform))
	return FALSE;
    } else {
      if (!jtransform_perfect_transform(srcinfo->output_width,
	  srcinfo->output_height,
	  srcinfo->max_h_samp_factor * srcinfo->min_DCT_h_scaled_size,
	  srcinfo->max_v_samp_factor * srcinfo->min_DCT_v_scaled_size,
	  info->transform))
	return FALSE;
    }
  }

  /* If there is only one output component, force the iMCU size to be 1;
   * else use the source iMCU size.  (This allows us to do the right thing
   * when reducing color to grayscale, and also provides a handy way of
   * cleaning up "funny" grayscale images whose sampling factors are not 1x1.)
   */
  switch (info->transform) {
  case JXFORM_TRANSPOSE:
  case JXFORM_TRANSVERSE:
  case JXFORM_ROT_90:
  case JXFORM_ROT_270:
    info->output_width = srcinfo->output_height;
    info->output_height = srcinfo->output_width;
    if (info->num_components == 1) {
      info->iMCU_sample_width = srcinfo->min_DCT_v_scaled_size;
      info->iMCU_sample_height = srcinfo->min_DCT_h_scaled_size;
    } else {
      info->iMCU_sample_width =
	srcinfo->max_v_samp_factor * srcinfo->min_DCT_v_scaled_size;
      info->iMCU_sample_height =
	srcinfo->max_h_samp_factor * srcinfo->min_DCT_h_scaled_size;
    }
    break;
  default:
    info->output_width = srcinfo->output_width;
    info->output_height = srcinfo->output_height;
    if (info->num_components == 1) {
      info->iMCU_sample_width = srcinfo->min_DCT_h_scaled_size;
      info->iMCU_sample_height = srcinfo->min_DCT_v_scaled_size;
    } else {
      info->iMCU_sample_width =
	srcinfo->max_h_samp_factor * srcinfo->min_DCT_h_scaled_size;
      info->iMCU_sample_height =
	srcinfo->max_v_samp_factor * srcinfo->min_DCT_v_scaled_size;
    }
    break;
  }

  /* If cropping has been requested, compute the crop area's position and
   * dimensions, ensuring that its upper left corner falls at an iMCU boundary.
   */
  if (info->crop) {
    /* Insert default values for unset crop parameters */
    if (info->crop_xoffset_set == JCROP_UNSET)
      info->crop_xoffset = 0;	/* default to +0 */
    if (info->crop_yoffset_set == JCROP_UNSET)
      info->crop_yoffset = 0;	/* default to +0 */
    if (info->crop_xoffset >= info->output_width ||
	info->crop_yoffset >= info->output_height)
      ERREXIT(srcinfo, JERR_BAD_CROP_SPEC);
    if (info->crop_width_set == JCROP_UNSET)
      info->crop_width = info->output_width - info->crop_xoffset;
    if (info->crop_height_set == JCROP_UNSET)
      info->crop_height = info->output_height - info->crop_yoffset;
    /* Ensure parameters are valid */
    if (info->crop_width <= 0 || info->crop_width > info->output_width ||
	info->crop_height <= 0 || info->crop_height > info->output_height ||
	info->crop_xoffset > info->output_width - info->crop_width ||
	info->crop_yoffset > info->output_height - info->crop_height)
      ERREXIT(srcinfo, JERR_BAD_CROP_SPEC);
    /* Convert negative crop offsets into regular offsets */
    if (info->crop_xoffset_set == JCROP_NEG)
      xoffset = info->output_width - info->crop_width - info->crop_xoffset;
    else
      xoffset = info->crop_xoffset;
    if (info->crop_yoffset_set == JCROP_NEG)
      yoffset = info->output_height - info->crop_height - info->crop_yoffset;
    else
      yoffset = info->crop_yoffset;
    /* Now adjust so that upper left corner falls at an iMCU boundary */
    info->output_width =
      info->crop_width + (xoffset % info->iMCU_sample_width);
    info->output_height =
      info->crop_height + (yoffset % info->iMCU_sample_height);
    /* Save x/y offsets measured in iMCUs */
    info->x_crop_offset = xoffset / info->iMCU_sample_width;
    info->y_crop_offset = yoffset / info->iMCU_sample_height;
  } else {
    info->x_crop_offset = 0;
    info->y_crop_offset = 0;
  }

  /* Figure out whether we need workspace arrays,
   * and if so whether they are transposed relative to the source.
   */
  need_workspace = FALSE;
  transpose_it = FALSE;
  switch (info->transform) {
  case JXFORM_NONE:
    if (info->x_crop_offset != 0 || info->y_crop_offset != 0)
      need_workspace = TRUE;
    /* No workspace needed if neither cropping nor transforming */
    break;
  case JXFORM_FLIP_H:
    if (info->trim)
      trim_right_edge(info, srcinfo->output_width);
    if (info->y_crop_offset != 0)
      need_workspace = TRUE;
    /* do_flip_h_no_crop doesn't need a workspace array */
    break;
  case JXFORM_FLIP_V:
    if (info->trim)
      trim_bottom_edge(info, srcinfo->output_height);
    /* Need workspace arrays having same dimensions as source image. */
    need_workspace = TRUE;
    break;
  case JXFORM_TRANSPOSE:
    /* transpose does NOT have to trim anything */
    /* Need workspace arrays having transposed dimensions. */
    need_workspace = TRUE;
    transpose_it = TRUE;
    break;
  case JXFORM_TRANSVERSE:
    if (info->trim) {
      trim_right_edge(info, srcinfo->output_height);
      trim_bottom_edge(info, srcinfo->output_width);
    }
    /* Need workspace arrays having transposed dimensions. */
    need_workspace = TRUE;
    transpose_it = TRUE;
    break;
  case JXFORM_ROT_90:
    if (info->trim)
      trim_right_edge(info, srcinfo->output_height);
    /* Need workspace arrays having transposed dimensions. */
    need_workspace = TRUE;
    transpose_it = TRUE;
    break;
  case JXFORM_ROT_180:
    if (info->trim) {
      trim_right_edge(info, srcinfo->output_width);
      trim_bottom_edge(info, srcinfo->output_height);
    }
    /* Need workspace arrays having same dimensions as source image. */
    need_workspace = TRUE;
    break;
  case JXFORM_ROT_270:
    if (info->trim)
      trim_bottom_edge(info, srcinfo->output_width);
    /* Need workspace arrays having transposed dimensions. */
    need_workspace = TRUE;
    transpose_it = TRUE;
    break;
  }

  /* Allocate workspace if needed.
   * Note that we allocate arrays padded out to the next iMCU boundary,
   * so that transform routines need not worry about missing edge blocks.
   */
  if (need_workspace) {
    coef_arrays = (jvirt_barray_ptr *)
      (*srcinfo->mem->alloc_small) ((j_common_ptr) srcinfo, JPOOL_IMAGE,
		SIZEOF(jvirt_barray_ptr) * info->num_components);
    width_in_iMCUs = (JDIMENSION)
      jdiv_round_up((long) info->output_width,
		    (long) info->iMCU_sample_width);
    height_in_iMCUs = (JDIMENSION)
      jdiv_round_up((long) info->output_height,
		    (long) info->iMCU_sample_height);
    for (ci = 0; ci < info->num_components; ci++) {
      compptr = srcinfo->comp_info + ci;
      if (info->num_components == 1) {
	/* we're going to force samp factors to 1x1 in this case */
	h_samp_factor = v_samp_factor = 1;
      } else if (transpose_it) {
	h_samp_factor = compptr->v_samp_factor;
	v_samp_factor = compptr->h_samp_factor;
      } else {
	h_samp_factor = compptr->h_samp_factor;
	v_samp_factor = compptr->v_samp_factor;
      }
      width_in_blocks = width_in_iMCUs * h_samp_factor;
      height_in_blocks = height_in_iMCUs * v_samp_factor;
      coef_arrays[ci] = (*srcinfo->mem->request_virt_barray)
	((j_common_ptr) srcinfo, JPOOL_IMAGE, FALSE,
	 width_in_blocks, height_in_blocks, (JDIMENSION) v_samp_factor);
    }
    info->workspace_coef_arrays = coef_arrays;
  } else
    info->workspace_coef_arrays = NULL;

  return TRUE;
}


/* Transpose destination image parameters */

LOCAL(void)
transpose_critical_parameters (j_compress_ptr dstinfo)
{
  int tblno, i, j, ci, itemp;
  jpeg_component_info *compptr;
  JQUANT_TBL *qtblptr;
  JDIMENSION jtemp;
  UINT16 qtemp;

  /* Transpose image dimensions */
  jtemp = dstinfo->image_width;
  dstinfo->image_width = dstinfo->image_height;
  dstinfo->image_height = jtemp;
  itemp = dstinfo->min_DCT_h_scaled_size;
  dstinfo->min_DCT_h_scaled_size = dstinfo->min_DCT_v_scaled_size;
  dstinfo->min_DCT_v_scaled_size = itemp;

  /* Transpose sampling factors */
  for (ci = 0; ci < dstinfo->num_components; ci++) {
    compptr = dstinfo->comp_info + ci;
    itemp = compptr->h_samp_factor;
    compptr->h_samp_factor = compptr->v_samp_factor;
    compptr->v_samp_factor = itemp;
  }

  /* Transpose quantization tables */
  for (tblno = 0; tblno < NUM_QUANT_TBLS; tblno++) {
    qtblptr = dstinfo->quant_tbl_ptrs[tblno];
    if (qtblptr != NULL) {
      for (i = 0; i < DCTSIZE; i++) {
	for (j = 0; j < i; j++) {
	  qtemp = qtblptr->quantval[i*DCTSIZE+j];
	  qtblptr->quantval[i*DCTSIZE+j] = qtblptr->quantval[j*DCTSIZE+i];
	  qtblptr->quantval[j*DCTSIZE+i] = qtemp;
	}
      }
    }
  }
}

#define EXIF_FOUND_WIDTH  1
#define EXIF_FOUND_HEIGHT 2
#define EXIF_FOUND_ORIENT 4
#define EXIF_FOUND_ALL (1|2|4)

static ExifEntry * exif_content_add_tag_intvalue( ExifContent * c, ExifTag tag, int value ) {
	if (!c || !c->parent) return NULL;
	ExifEntry *e;
	ExifByteOrder o = exif_data_get_byte_order(c->parent);
	e = exif_entry_new();
	exif_content_add_entry(c, e);
	exif_entry_unref(e);
	exif_entry_initialize(e, tag);
	switch (e->format) {
		case EXIF_FORMAT_LONG:
			exif_set_long (e->data, o, value);
			break;
		case EXIF_FORMAT_SLONG:
			exif_set_slong (e->data, o, value);
			break;
		case EXIF_FORMAT_BYTE:
			e->data[0] = (unsigned char) value;
			break;
		case EXIF_FORMAT_SBYTE:
			e->data[0] = (char) value;
			break;
		case EXIF_FORMAT_SSHORT:
			exif_set_sshort (e->data, o, value);
			break;
		case EXIF_FORMAT_SHORT:
			exif_set_short (e->data, o, value);
			break;
		default:
			fprintf(stderr, "Bad format for `exif_content_add_tag_intvalue': %d\n",e->format);
			return NULL;
	}
	exif_entry_fix(e);
	return e;
}

static ExifEntry * exif_entry_set_intvalue( ExifEntry * e, ExifByteOrder o, int value ) {
	if (!e || !e->parent || !e->parent->parent) return NULL;
	switch (e->format) {
		case EXIF_FORMAT_LONG:
			exif_set_long (e->data, o, value);
			break;
		case EXIF_FORMAT_SLONG:
			exif_set_slong (e->data, o, value);
			break;
		case EXIF_FORMAT_BYTE:
			e->data[0] = (unsigned char) value;
			break;
		case EXIF_FORMAT_SBYTE:
			e->data[0] = (char) value;
			break;
		case EXIF_FORMAT_SSHORT:
			exif_set_sshort (e->data, o, value);
			break;
		case EXIF_FORMAT_SHORT:
			exif_set_short (e->data, o, value);
			break;
		default:
			fprintf(stderr, "Bad format for `exif_content_add_tag_intvalue': %d\n",e->format);
			return NULL;
	}
	exif_entry_fix(e);
	return e;
}


//JXFORM_FLIP_H	JXFORM_FLIP_V	JXFORM_TRANSPOSE	JXFORM_TRANSVERSE	JXFORM_ROT_90	JXFORM_ROT_180	JXFORM_ROT_270
//normal
//flipped h
//rot 180
//flipped v
//transposed
//rotated 270
//transversed
//rotated 90

static char transorient[9][9] = {
	{0,	0,	0,	0,	0,	0,	0,	0},
	{0,	2,	4,	5,	7,	8,	3,	6},
	{0,	1,	3,	6,	8,	7,	4,	5},
	{0,	4,	2,	7,	5,	6,	1,	8},
	{0,	3,	1,	8,	6,	5,	2,	7},
	{0,	8,	6,	1,	3,	2,	7,	4},
	{0,	7,	5,	2,	4,	1,	8,	3},
	{0,	6,	8,	3,	1,	4,	5,	2},
	{0,	5,	7,	4,	2,	3,	6,	1}
};

LOCAL(void)
adjust_exif (
	jpeg_saved_marker_ptr   exif,
	j_compress_ptr          dstinfo,
	jpeg_transform_info    *info
) {
	
	//fprintf(stderr,"[%02x] %02x %02x %02x %02x .(%d). %02x %02x\n", data[-1], data[0],data[1],data[2],data[3],length,data[length-2],data[length-1]);
	ExifData    * edata;
	ExifContent * content;
	ExifEntry   * entry, *width, *height;
	
	int k,l, orientation, found = 0;
	unsigned char val[255];
	ExifByteOrder byteorder;
	
	edata = exif_data_new();
	exif_data_load_data( edata, (const unsigned char *)exif->data, exif->data_length );
	exif_data_fix(edata);

	if (edata) {
		//exif_data_dump(edata);
		byteorder = exif_data_get_byte_order(edata);
		//fprintf(stderr,"adjust exif with byteorder %s\n",exif_byte_order_get_name(byteorder));
		
		for (k = 0; k < EXIF_IFD_COUNT; k++) {
			content = edata->ifd[k];
			if (!content) break;
			for (l = 0; l < content->count; l++) {
				entry = content->entries[l];
				//warn("entry %08x -> %s = %s",entry->tag, exif_tag_get_name(entry->tag), exif_entry_get_value(entry,val,255));
				if (entry->tag == 0x112 ) {
					switch(entry->format) {
						case EXIF_FORMAT_SSHORT:
						case EXIF_FORMAT_BYTE:
						case EXIF_FORMAT_SBYTE:
						case EXIF_FORMAT_SHORT:
							//fprintf(stderr,"Get short\n");
							orientation = exif_get_short(entry->data,byteorder);
							break;
						case EXIF_FORMAT_LONG:
						case EXIF_FORMAT_SLONG:
							//fprintf(stderr,"Get long\n");
							orientation = exif_get_long(entry->data,byteorder);
							break;
						default:
							fprintf(stderr,"Bad format for exif orientation: %s. Reset orientation to normal\n", exif_format_get_name(entry->format));
							orientation = 0;
							exif_content_remove_entry(content,entry);
							break;
					}
					if (orientation < 1 || orientation > 8) {
						fprintf(stderr,"Bad value for orientation: %d. Reset orientation to normal\n", orientation);
						exif_content_remove_entry(content,entry);
						continue;
					}
					found |= EXIF_FOUND_ORIENT;
					//exif_entry_dump(entry,1);
					char neworient = transorient[orientation][info->transform];
					if (!neworient) neworient = 1;
					exif_entry_set_intvalue(entry, byteorder, neworient);
					//fprintf(stderr,"(%s) %d -> (%d) -> %d\n",exif_format_get_name(entry->format),
					//	orientation, info->transform, neworient);
				}
				else
				if (entry->tag == 0x100) {
					found |= EXIF_FOUND_WIDTH;
					exif_entry_set_intvalue(entry, byteorder, dstinfo->jpeg_width);
				}
				else
				if(entry->tag == 0x101) {
					found |= EXIF_FOUND_HEIGHT;
					exif_entry_set_intvalue(entry, byteorder, dstinfo->jpeg_height);
				}
				if(found == EXIF_FOUND_ALL) {
					k = EXIF_IFD_COUNT;
					break;
				}
			}
		}
		if (found != EXIF_FOUND_ALL) {
			//fprintf(stderr,"Found not all entries\n");
			ExifEntry *e;
			if (! ( found & EXIF_FOUND_ORIENT )) {
				fprintf(stderr,"Create orientation = 1\n");
				e = exif_content_add_tag_intvalue( edata->ifd[0], EXIF_TAG_ORIENTATION, transorient[1][info->transform] );
				exif_entry_dump(e,1);
			}
			if (! ( found & EXIF_FOUND_WIDTH )) {
				e = exif_content_add_tag_intvalue( edata->ifd[0], EXIF_TAG_IMAGE_WIDTH,  dstinfo->jpeg_width );
				//exif_entry_dump(e,1);
			}
			if (! ( found & EXIF_FOUND_HEIGHT )) {
				e = exif_content_add_tag_intvalue( edata->ifd[0], EXIF_TAG_IMAGE_LENGTH, dstinfo->jpeg_height );
				//exif_entry_dump(e,1);
			}
		}
		
		if (info->discard_thumbnail) {
			ExifData    * newdata;
			newdata = exif_data_new();
			exif_data_set_byte_order(newdata,byteorder);
			int i = 0;
			for(i=0;i < EXIF_IFD_COUNT; i++) {
				if( newdata->ifd[i] = edata->ifd[i] ) {
					newdata->ifd[i]->parent = newdata;
					exif_content_ref(newdata->ifd[i]);
				}
			}
			exif_data_unref(edata);
			edata = newdata;
		}
		//exif_data_dump(newdata);
		
		JOCTET FAR * newexif;
		unsigned int newlength;
		exif_data_save_data(edata, &newexif, &newlength);
		//fprintf(stderr,"created exif: %d (was %d) -> %p (was %p)\n",newlength, exif->data_length, newexif, exif->data);
		exif->data = newexif;
		exif->data_length = newlength;
		return;
	} else {
		fprintf(stderr,"bad edata\n");
	}
	return;
}



/* Adjust output image parameters as needed.
 *
 * This must be called after jpeg_copy_critical_parameters()
 * and before jpeg_write_coefficients().
 *
 * The return value is the set of virtual coefficient arrays to be written
 * (either the ones allocated by jtransform_request_workspace, or the
 * original source data arrays).  The caller will need to pass this value
 * to jpeg_write_coefficients().
 */

GLOBAL(jvirt_barray_ptr *)
jtransform_adjust_parameters (
	j_decompress_ptr      srcinfo,
	j_compress_ptr        dstinfo,
	jvirt_barray_ptr     *src_coef_arrays,
	jpeg_transform_info  *info
) {
  /* If force-to-grayscale is requested, adjust destination parameters */
  if (info->force_grayscale) {
    /* First, ensure we have YCbCr or grayscale data, and that the source's
     * Y channel is full resolution.  (No reasonable person would make Y
     * be less than full resolution, so actually coping with that case
     * isn't worth extra code space.  But we check it to avoid crashing.)
     */
    if (((dstinfo->jpeg_color_space == JCS_YCbCr &&
	  dstinfo->num_components == 3) ||
	 (dstinfo->jpeg_color_space == JCS_GRAYSCALE &&
	  dstinfo->num_components == 1)) &&
	srcinfo->comp_info[0].h_samp_factor == srcinfo->max_h_samp_factor &&
	srcinfo->comp_info[0].v_samp_factor == srcinfo->max_v_samp_factor) {
      /* We use jpeg_set_colorspace to make sure subsidiary settings get fixed
       * properly.  Among other things, it sets the target h_samp_factor &
       * v_samp_factor to 1, which typically won't match the source.
       * We have to preserve the source's quantization table number, however.
       */
      int sv_quant_tbl_no = dstinfo->comp_info[0].quant_tbl_no;
      jpeg_set_colorspace(dstinfo, JCS_GRAYSCALE);
      dstinfo->comp_info[0].quant_tbl_no = sv_quant_tbl_no;
    } else {
      /* Sorry, can't do it */
      ERREXIT(dstinfo, JERR_CONVERSION_NOTIMPL);
    }
  } else if (info->num_components == 1) {
    /* For a single-component source, we force the destination sampling factors
     * to 1x1, with or without force_grayscale.  This is useful because some
     * decoders choke on grayscale images with other sampling factors.
     */
    dstinfo->comp_info[0].h_samp_factor = 1;
    dstinfo->comp_info[0].v_samp_factor = 1;
  }

  /* Correct the destination's image dimensions as necessary
   * for rotate/flip, resize, and crop operations.
   */
  dstinfo->jpeg_width = info->output_width;
  dstinfo->jpeg_height = info->output_height;

  /* Transpose destination image parameters */
  switch (info->transform) {
  case JXFORM_TRANSPOSE:
  case JXFORM_TRANSVERSE:
  case JXFORM_ROT_90:
  case JXFORM_ROT_270:
    transpose_critical_parameters(dstinfo);
    break;
  default:
    break;
  }

	if (srcinfo->marker_list != NULL) {
		jpeg_saved_marker_ptr marker;
		marker = srcinfo->marker_list;
		while (marker) {
			//fprintf(stderr,"have marker %d (length = %d)\n",marker->marker,marker->data_length );
			if (
				marker->marker == JPEG_APP0+1 &&
				marker->data_length >= 6 && match(marker->data,JExifHeader)
			) {
				// Adjust Exif properties
				//fprintf(stderr,"found exif marker\n");
				adjust_exif(marker, dstinfo, info);
				// Suppress output of JFIF marker
				dstinfo->write_JFIF_header = FALSE;
				break;
			}
			marker = marker->next;
		}
	}
	
	/* Return the appropriate output data set */
	if (info->workspace_coef_arrays != NULL) return info->workspace_coef_arrays;
	return src_coef_arrays;
}


/* Execute the actual transformation, if any.
 *
 * This must be called *after* jpeg_write_coefficients, because it depends
 * on jpeg_write_coefficients to have computed subsidiary values such as
 * the per-component width and height fields in the destination object.
 *
 * Note that some transformations will modify the source data arrays!
 */

GLOBAL(void)
jtransform_execute_transform (j_decompress_ptr srcinfo,
			      j_compress_ptr dstinfo,
			      jvirt_barray_ptr *src_coef_arrays,
			      jpeg_transform_info *info)
{
  jvirt_barray_ptr *dst_coef_arrays = info->workspace_coef_arrays;

  /* Note: conditions tested here should match those in switch statement
   * in jtransform_request_workspace()
   */
  switch (info->transform) {
  case JXFORM_NONE:
    if (info->x_crop_offset != 0 || info->y_crop_offset != 0)
      do_crop(srcinfo, dstinfo, info->x_crop_offset, info->y_crop_offset,
	      src_coef_arrays, dst_coef_arrays);
    break;
  case JXFORM_FLIP_H:
    if (info->y_crop_offset != 0)
      do_flip_h(srcinfo, dstinfo, info->x_crop_offset, info->y_crop_offset,
		src_coef_arrays, dst_coef_arrays);
    else
      do_flip_h_no_crop(srcinfo, dstinfo, info->x_crop_offset,
			src_coef_arrays);
    break;
  case JXFORM_FLIP_V:
    do_flip_v(srcinfo, dstinfo, info->x_crop_offset, info->y_crop_offset,
	      src_coef_arrays, dst_coef_arrays);
    break;
  case JXFORM_TRANSPOSE:
    do_transpose(srcinfo, dstinfo, info->x_crop_offset, info->y_crop_offset,
		 src_coef_arrays, dst_coef_arrays);
    break;
  case JXFORM_TRANSVERSE:
    do_transverse(srcinfo, dstinfo, info->x_crop_offset, info->y_crop_offset,
		  src_coef_arrays, dst_coef_arrays);
    break;
  case JXFORM_ROT_90:
    do_rot_90(srcinfo, dstinfo, info->x_crop_offset, info->y_crop_offset,
	      src_coef_arrays, dst_coef_arrays);
    break;
  case JXFORM_ROT_180:
    do_rot_180(srcinfo, dstinfo, info->x_crop_offset, info->y_crop_offset,
	       src_coef_arrays, dst_coef_arrays);
    break;
  case JXFORM_ROT_270:
    do_rot_270(srcinfo, dstinfo, info->x_crop_offset, info->y_crop_offset,
	       src_coef_arrays, dst_coef_arrays);
    break;
  }
}

/* jtransform_perfect_transform
 *
 * Determine whether lossless transformation is perfectly
 * possible for a specified image and transformation.
 *
 * Inputs:
 *   image_width, image_height: source image dimensions.
 *   MCU_width, MCU_height: pixel dimensions of MCU.
 *   transform: transformation identifier.
 * Parameter sources from initialized jpeg_struct
 * (after reading source header):
 *   image_width = cinfo.image_width
 *   image_height = cinfo.image_height
 *   MCU_width = cinfo.max_h_samp_factor * cinfo.block_size
 *   MCU_height = cinfo.max_v_samp_factor * cinfo.block_size
 * Result:
 *   TRUE = perfect transformation possible
 *   FALSE = perfect transformation not possible
 *           (may use custom action then)
 */

GLOBAL(boolean)
jtransform_perfect_transform(
	JDIMENSION image_width,
	JDIMENSION image_height,
	int MCU_width,
	int MCU_height,
	JXFORM_CODE transform
) {
  boolean result = TRUE; /* initialize TRUE */

  switch (transform) {
  case JXFORM_FLIP_H:
  case JXFORM_ROT_270:
    if (image_width % (JDIMENSION) MCU_width)
      result = FALSE;
    break;
  case JXFORM_FLIP_V:
  case JXFORM_ROT_90:
    if (image_height % (JDIMENSION) MCU_height)
      result = FALSE;
    break;
  case JXFORM_TRANSVERSE:
  case JXFORM_ROT_180:
    if (image_width % (JDIMENSION) MCU_width)
      result = FALSE;
    if (image_height % (JDIMENSION) MCU_height)
      result = FALSE;
    break;
  default:
    break;
  }

  return result;
}

/* Setup decompression object to save desired markers in memory.
 * This must be called before jpeg_read_header() to have the desired effect.
 */

GLOBAL(void)
jcopy_markers_setup (
	j_decompress_ptr srcinfo,
	JCOPY_MASK mask // see transsupp.h -> JCOPY_
) {
	mask &= JCOPY_ALL;
	if (!mask) return;
	unsigned int i, x = 1;
	jpeg_save_markers(srcinfo, JPEG_APP0, 0xFFFF); // always copy APP0;
	for (i=1; i < 33; i++) {
		if ( mask & x ) {
			//fprintf(stderr,"save marker APP0+%d\n",i);
			jpeg_save_markers(srcinfo, JPEG_APP0 + i, 0xFFFF);
		}
		x <<= 1;
	}
}

/* Copy markers saved in the given source object to the destination object.
 * This should be called just after jpeg_start_compress() or
 * jpeg_write_coefficients().
 * Note that those routines will have written the SOI, and also the
 * JFIF APP0 or Adobe APP14 markers if selected.
 */

GLOBAL(void)
jcopy_markers_execute (
	j_decompress_ptr srcinfo,
	j_compress_ptr dstinfo,
	JCOPY_MASK copymask
) {
	jpeg_saved_marker_ptr marker;

	/* In the current implementation, we don't actually need to examine the
	 * option flag here; we just copy everything that got saved.
	 * But to avoid confusion, we do not output JFIF and Adobe APP14 markers
	 * if the encoder library already wrote one.
	 */
	for (marker = srcinfo->marker_list; marker != NULL; marker = marker->next) {
		if (marker->marker == JPEG_APP0 + 1) {
			if ( marker->data_length > 6 && match(marker->data,JExifHeader) ) {
			} else {
				if (copymask != JCOPY_ALL) {
					//fprintf(stderr,"skip exif of type %02x %02x\n", marker->data[0],marker->data[1]);
					continue;
				}
			}
		}
		
		if (
			dstinfo->write_JFIF_header && marker->marker == JPEG_APP0 &&
			marker->data_length >= 5 && match(marker->data,JFIFHeader)
		) continue;			// reject duplicate JFIF
		if (
			dstinfo->write_Adobe_marker && marker->marker == JPEG_APP0+14 &&
			marker->data_length >= 5 && match(marker->data,JAdobeHeader)
		) continue;			// reject duplicate Adobe
		//fprintf(stderr,"write marker %04x %p (length = %d)\n",marker->marker,marker->data,marker->data_length);
		jpeg_write_marker(dstinfo, marker->marker, marker->data, marker->data_length);
	}
}
