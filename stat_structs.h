/* stat_structs.h -- Structs required for and shared between stats.c and bamcheck.c

    Copyright (C) 2015 Genome Research Ltd.

    Author: Sam Nicholls <sam@samnicholls.net>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#ifndef __SAMTOOLS_STRUCTS_H__
#define __SAMTOOLS_STRUCTS_H__

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include "stats_isize.h"

// The GC-depth graph works as follows: split the reference sequence into
// segments and calculate GC content and depth in each bin. Then sort
// these segments by their GC and plot the depth distribution by means
// of 10th, 25th, etc. depth percentiles.
typedef struct
{
    float gc;
    uint32_t depth;
}
gc_depth_t;

// For coverage distribution, a simple pileup
typedef struct
{
    int64_t pos;
    int size, start;
    int *buffer;
}
round_buffer_t;

typedef struct { uint32_t from, to; } pos_t;
typedef struct
{
    int npos,mpos,cpos;
    pos_t *pos;
}
regions_t;

typedef struct
{
    uint64_t a;
    uint64_t c;
    uint64_t g;
    uint64_t t;
    uint64_t n;
    uint64_t other;
}
acgtno_count_t;

typedef struct
{
    // Auxiliary data
    int flag_require, flag_filter;
    faidx_t *fai;                   // Reference sequence for GC-depth graph
    int argc;                       // Command line arguments to be printed on the output
    char **argv;
    int gcd_bin_size;           // The size of GC-depth bin
    int nisize;         // The maximum insert size that the allocated array can hold - 0 indicates no limit
    int trim_qual;      // bwa trim quality
    float isize_main_bulk;  // There are always some unrealistically big insert sizes, report only the main part
    int cov_min,cov_max,cov_step;   // Minimum, maximum coverage and size of the coverage bins
    samFile* sam;
    bam_hdr_t* sam_header;

    // Filters
    int filter_readlen;

    // Misc
    char *split_tag;      // Tag on which to perform stats splitting
    char *split_prefix;   // Path or string prefix for filenames created when splitting
}
stats_info_t;

typedef struct {
    double mean_above_baseline;
    double mean_below_baseline;
    double max_above_baseline;
    double max_below_baseline;
    double max_baseline_deviation;
    double total_mean_deviation;
} bamcheck_bcd_t;

typedef struct {
    uint64_t iqr_inc_contig_length;
    uint64_t iqr_inc_contig_start;

    uint64_t runmed_mean_dec_contig_length;
    uint64_t runmed_mean_dec_contig_start;

    float runmed_mean_dec_high;
    float runmed_mean_dec_low;
    float runmed_mean_dec_range;

} bamcheck_quality_dropoff_t ;

typedef struct {
    double fwd_ins_above_baseline_pct;
    double fwd_ins_below_baseline_pct;

    double fwd_del_above_baseline_pct;
    double fwd_del_below_baseline_pct;

    double rev_ins_above_baseline_pct;
    double rev_ins_below_baseline_pct;

    double rev_del_above_baseline_pct;
    double rev_del_below_baseline_pct;

    bamcheck_bcd_t *bcd_a;
    bamcheck_bcd_t *bcd_c;
    bamcheck_bcd_t *bcd_g;
    bamcheck_bcd_t *bcd_t;

    bamcheck_quality_dropoff_t *fwd_dropoff;
    bamcheck_quality_dropoff_t *rev_dropoff;

} bamcheck_stats_t;


typedef struct {
    uint64_t n;
    double *delta;

    double above_total;
    uint64_t above_n;
    double above_max;
    double above_min;

    double below_total;
    uint64_t below_n;
    double below_max;
    double below_min;

    uint64_t total_count;

    double max_baseline_deviation;
    double total_mean_deviation;
} bamcheck_baseline_delta;

typedef struct {
    double pct_above_baseline;
    double pct_below_baseline;
    uint64_t total_count;
} bamcheck_cycles_t;

typedef struct
{
    // Dimensions of the quality histogram holder (quals_1st,quals_2nd), GC content holder (gc_1st,gc_2nd),
    //  insert size histogram holder
    int nquals;         // The number of quality bins
    int nbases;         // The maximum sequence length the allocated array can hold
    int ngc;            // The size of gc_1st and gc_2nd
    int nindels;        // The maximum indel length for indel distribution

    // Arrays for the histogram data
    uint64_t *quals_1st, *quals_2nd;
    uint64_t *gc_1st, *gc_2nd;
    acgtno_count_t *acgtno_cycles;
    uint64_t *read_lengths;
    uint64_t *insertions, *deletions;
    uint64_t *ins_cycles_1st, *ins_cycles_2nd, *del_cycles_1st, *del_cycles_2nd;
    isize_t *isize;

    // The extremes encountered
    int max_len;            // Maximum read length
    int max_qual;           // Maximum quality
    int is_sorted;

    // Summary numbers
    uint64_t total_len;
    uint64_t total_len_dup;
    uint64_t nreads_1st;
    uint64_t nreads_2nd;
    uint64_t nreads_filtered;
    uint64_t nreads_dup;
    uint64_t nreads_unmapped;
    uint64_t nreads_single_mapped;
    uint64_t nreads_paired_and_mapped;
    uint64_t nreads_properly_paired;
    uint64_t nreads_paired_tech;
    uint64_t nreads_anomalous;
    uint64_t nreads_mq0;
    uint64_t nbases_mapped;
    uint64_t nbases_mapped_cigar;
    uint64_t nbases_trimmed;  // bwa trimmed bases
    uint64_t nmismatches;
    uint64_t nreads_QCfailed, nreads_secondary;
    struct {
        uint32_t names, reads, quals;
    } checksum;

    // GC-depth related data
    uint32_t ngcd, igcd;        // The maximum number of GC depth bins and index of the current bin
    gc_depth_t *gcd;            // The GC-depth bins holder
    int32_t tid, gcd_pos;       // Position of the current bin
    int32_t pos;                // Position of the last read

    // Coverage distribution related data
    int ncov;                       // The number of coverage bins
    uint64_t *cov;                  // The coverage frequencies
    round_buffer_t cov_rbuf;        // Pileup round buffer

    // Mismatches by read cycle
    uint8_t *rseq_buf;              // A buffer for reference sequence to check the mismatches against
    int mrseq_buf;                  // The size of the buffer
    int32_t rseq_pos;               // The coordinate of the first base in the buffer
    int32_t nrseq_buf;              // The used part of the buffer
    uint64_t *mpc_buf;              // Mismatches per cycle

    // Target regions
    int nregions, reg_from,reg_to;
    regions_t *regions;

    // Auxiliary data
    double sum_qual;                // For calculating average quality value
    void *rg_hash;                  // Read groups to include, the array is null-terminated

    // Split
    char* split_name;

    stats_info_t* info;             // Pointer to options and settings struct
    bamcheck_stats_t *bamcheck;     // Pointer to bamcheck struct

}
stats_t;

#endif
