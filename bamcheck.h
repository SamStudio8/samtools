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

bamcheck_bcd_t* bamcheck_base_content_baseline(double *base_prop, int n);
bamcheck_cycles_t* bamcheck_cycles(uint64_t *cycles_arr, int n, int k);
bamcheck_baseline_delta* bamcheck_baseline_d(uint64_t *baseline, uint64_t *count, int baseline_n, int count_n, double scalar_baseline);
bamcheck_baseline_delta* bamcheck_baseline_d_double(double *count, int count_n, double baseline);
bamcheck_baseline_delta* init_bamcheck_baseline_delta(int n);
uint64_t* runmed(uint64_t *cycles_arr, int n, int k);

void bamcheck_quality_dropoff(stats_t *curr_stats);
void bamcheck_base_content_deviation(stats_t *curr_stats);
void bamcheck_indel_peaks(stats_t *curr_stats);
double percentile(double *values, uint64_t n, int tile);

#endif
