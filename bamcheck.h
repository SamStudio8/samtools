/* bamcheck.c -- A port of bamcheckR to samtools

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

#ifndef __SAMTOOLS_BAMCHECK_H__
#define __SAMTOOLS_BAMCHECK_H__

#include <stdint.h>
#include "stat_structs.h"

bamcheck_stats_t* bamcheck_stats_init();
void calculate_bamcheck(stats_t *curr_stats);

int64_t med_triplet(int64_t a, int64_t b, int64_t c);
int cmpfunc (const void * a, const void * b);
uint64_t* copy_arr(uint64_t *source, int n, int new_n, int filter);

void bamcheck_stats_destroy(bamcheck_stats_t *bamcheck);
void bamcheck_base_content_baseline(bamcheck_bcd_t *result, double *base_prop, int n);
void bamcheck_cycles(bamcheck_cycles_t *result, uint64_t *cycles_arr, int n, int k);
void bamcheck_baseline_d(bamcheck_baseline_delta *result, uint64_t *baseline, uint64_t *count, int baseline_n, int count_n, double scalar_baseline);
void bamcheck_baseline_d_double(bamcheck_baseline_delta *result, double *count, int count_n, double baseline);
void init_bamcheck_baseline_delta(bamcheck_baseline_delta *result, int n);
uint64_t* runmed(uint64_t *cycles_arr, int n, int k);

void bamcheck_quality_dropoff(stats_t *curr_stats);
void bamcheck_base_content_deviation(stats_t *curr_stats);
void bamcheck_indel_peaks(stats_t *curr_stats);
double percentile(double *values, uint64_t n, int tile);
double* copy_arr_double(double *source, int n, int new_n, int filter);
double* runmed_double(double *cycles_arr, int n, int k);
uint64_t wtd_percentile(double *cumsums, uint64_t cumsums_n, uint64_t n, int tile);
int cmpfunc_f (const void * a, const void * b);
double med_triplet_double(double a, double b, double c);
void bamcheck_quality_dropoff(stats_t *curr_stats);
void bamcheck_quality_dropoff_executor(bamcheck_quality_dropoff_t *result, uint64_t *cycle_counts, uint64_t cycles_n, uint64_t n_bin_quals, uint64_t max_quals_n);
void summarise_cycles(bamcheck_cycles_summary_t *summary, uint64_t *cycle_counts, uint64_t cycles_n, uint64_t n_bin_quals, uint64_t max_quals_n);
void init_bamcheck_cycles_summary(bamcheck_cycles_summary_t *summary, uint64_t cycles_n);
void destroy_bamcheck_cycles_summary(bamcheck_cycles_summary_t *summary);

#endif
