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

#include "bamcheck.h"

#include <stdint.h>
#include <stdio.h>

int cmpfunc (const void * a, const void * b)
{
       return ( *(uint64_t*)a - *(uint64_t*)b );
}

uint64_t* copy_arr(uint64_t *source, int n, int new_n, int filter){
    uint64_t *dest;
    dest = calloc(new_n,sizeof(uint64_t));

    size_t i = 0;
    size_t j = 0;
    for (; i < n; i++) {
        if ( source[i] == 0 && filter > 0 ) {
            continue;
        }
        dest[j] = source[i];
        j++;
    }
    return dest;

}

int64_t med_triplet(int64_t a, int64_t b, int64_t c){
    if (a < b) {
        if (c < b){
            // C[A]B or A[C]B
            if (a >= c) return a;
            return c;
        }
    }
    else {
        if (c > b){
            // B[A]C or B[C]A
            if (a <= c) return a;
            return c;
        }
    }
    return b;
}

//TODO Super function to also return mean and median (seeing as they are trivial)
uint64_t* runmed(uint64_t *cycles_arr, int n, int k){
    uint64_t *count;
    uint64_t *baseline = calloc(n, sizeof(uint64_t));

    size_t i;
    for(i = 0; i < n; i++){
        baseline[i] = 0;
    }

    uint64_t smooth_start[k/2];
    uint64_t smooth_end[k/2];
    for(i = 0; i < k/2; i++){
        smooth_start[i] = 0;
        smooth_end[i] = 0;
    }

    // Copy counts for the known number of valid IC lines (n)
    count = copy_arr(cycles_arr, n, n, 0);

    // Medianize sliding windows of k, beginning at i=k/2 (the first valid full window)
    uint64_t *k_window;
    for(i = (k/2); i < n-(k/2); i++){
        k_window = copy_arr(&count[i-(k/2)], k, k, 0);
        qsort(k_window, k, sizeof(uint64_t), cmpfunc);
        baseline[i] = k_window[k/2];
        free(k_window);
    }

    // Keep ends ( k/2 elements on each end )
    memcpy(baseline, &count[0], sizeof(uint64_t) * (k/2));
    memcpy(&baseline[n-(k/2)], &count[n-(k/2)], sizeof(uint64_t) * (k/2));

    // Smooth 1st and n-1th element
    smooth_start[1] = med_triplet(baseline[0], baseline[1], baseline[2]);
    smooth_end[(k/2)-2] = med_triplet(baseline[n-1], baseline[n-2], baseline[n-3]);

    // Smooth remaining elements to k/2

    // Smooth values for first and last k/2 elements (where the windows were too
    //  large to medianize initially). Work inwards starting with element pair (2, n-2).
    for(i = 2; i < (k/2)+1; i ++){
        int j = 2 * i - 1;

        k_window = copy_arr(&baseline[0], j, j, 0);
        qsort(k_window, j, sizeof(uint64_t), cmpfunc);
        smooth_start[i-1] = k_window[j/2];
        free(k_window);

        k_window = copy_arr(&baseline[n-j], j, j, 0);
        qsort(k_window, j, sizeof(uint64_t), cmpfunc);
        smooth_end[(k/2)-i] = k_window[j/2];
        free(k_window);
    }

    // First and last element with Tukey Rule
    smooth_start[0] = med_triplet(baseline[0], smooth_start[1], 3 * (int)smooth_start[1] - 2 * (int)smooth_start[2]);
    smooth_end[(k/2)-1] = med_triplet(baseline[n-1], smooth_end[(k/2)-2], 3 * (int)smooth_end[(k/2)-2] - 2 * (int)smooth_end[(k/2)-3]);

    // Move smoothed values over baseline
    memcpy(baseline, &smooth_start[0], sizeof(uint64_t) * (k/2));
    memcpy(&baseline[n-(k/2)], &smooth_end[0], sizeof(uint64_t) * (k/2));

    free(count);
    return baseline;
}

bamcheck_baseline_delta* init_bamcheck_baseline_delta(int n){
    bamcheck_baseline_delta *delta = calloc(1, sizeof(bamcheck_baseline_delta));

    delta->delta = calloc(n, sizeof(double));
    delta->n = n;

    delta->above_total = 0;
    delta->below_total = 0;

    delta->above_n = 0;
    delta->below_n = 0;

    delta->above_min = 0;
    delta->below_min = 0;

    delta->above_max = 0;
    delta->below_max = 0;

    delta->total_count = 0;

    delta->max_baseline_deviation = 0;
    delta->total_mean_deviation = 0;
    return delta;
}

bamcheck_baseline_delta* bamcheck_baseline_d(uint64_t *baseline, uint64_t *count, int baseline_n, int count_n, double scalar_baseline){

    /*
    int scalar_baseline = 0;
    if ( baseline_n != count_n ){
        // If baseline_n and count_n are not the same dimension,
        // and baseline_n is not a scalar (ie. a mean or median)
        // then we can't know what to do with the baseline.
        if ( baseline_n != 1 ){
            //TODO Throw a proper samtools error
            exit(1);
        }
        else{
            scalar_baseline = 1;
        }
    }
    */

    bamcheck_baseline_delta *result;
    result = init_bamcheck_baseline_delta(count_n);

    int i;
    double curr_delta, above_min, above_max, below_min, below_max;
    above_min = above_max = below_min = below_max = -1;

    for(i = 0; i < count_n; i++){
        if ( baseline_n == 0 ){
            curr_delta = (double)count[i] - scalar_baseline;
        }
        else{
            curr_delta = (double)count[i] - baseline[i];
        }
        result->delta[i] = curr_delta;
        result->total_count += count[i];

        if (curr_delta > 0){
            result->above_n++;
            result->above_total += curr_delta;

            if( above_min == -1 ) above_min = curr_delta;
            else {
                if ( curr_delta < above_min ) above_min = curr_delta;
            }

            if( above_max == -1 ) above_max = curr_delta;
            else {
                if ( curr_delta > above_max ) above_max = curr_delta;
            }
        }
        else if (curr_delta < 0){
            curr_delta = curr_delta*-1;
            result->below_n++;
            result->below_total += curr_delta;

            if( below_min == -1 ) below_min = curr_delta;
            else {
                if ( curr_delta < below_min ) below_min = curr_delta;
            }

            if( below_max == -1 ) below_max = curr_delta;
            else {
                if ( curr_delta > below_max ) below_max = curr_delta;
            }
        }
    }

    result->below_min = below_min;
    result->below_max = below_max;

    result->above_min = above_min;
    result->above_max = above_max;

    return result;
}

bamcheck_cycles_t* bamcheck_cycles(uint64_t *cycles_arr, int n, int k){
    bamcheck_cycles_t *result = calloc(1, sizeof(bamcheck_cycles_t));

    uint64_t *count;
    count = copy_arr(cycles_arr, n, n, 0);
    uint64_t *baseline;
    baseline = runmed(cycles_arr, n, k);

    // Calculate baseline and counts above/below
    bamcheck_baseline_delta* deviation;
    deviation = bamcheck_baseline_d(baseline, count, n, n, 0);

    result->pct_above_baseline = (double)deviation->above_total / (double)deviation->total_count * 100;
    result->pct_below_baseline = (double)deviation->below_total / (double)deviation->total_count * 100;
    result->total_count = deviation->total_count;

    free(deviation->delta);
    free(deviation);
    free(count);
    free(baseline);

    return result;
}

void bamcheck_indel_peaks(stats_t *curr_stats){

    int k = 25;
    int ilen;
    int ic_lines = 0;
    for (ilen=0; ilen<=curr_stats->nbases; ilen++){
        if ( curr_stats->ins_cycles_1st[ilen]>0 || curr_stats->ins_cycles_2nd[ilen]>0 || curr_stats->del_cycles_1st[ilen]>0 || curr_stats->del_cycles_2nd[ilen]>0 ){
            ic_lines++;
        }
    }

    if (ic_lines > 0){
        bamcheck_cycles_t *result;
        result = bamcheck_cycles(curr_stats->ins_cycles_1st, ic_lines, k);
        curr_stats->bamcheck->fwd_ins_above_baseline_pct = result->pct_above_baseline;
        curr_stats->bamcheck->fwd_ins_below_baseline_pct = result->pct_below_baseline;
        free(result);

        result = bamcheck_cycles(curr_stats->ins_cycles_2nd, ic_lines, k);
        curr_stats->bamcheck->rev_ins_above_baseline_pct = result->pct_above_baseline;
        curr_stats->bamcheck->rev_ins_below_baseline_pct = result->pct_below_baseline;
        free(result);

        result = bamcheck_cycles(curr_stats->del_cycles_1st, ic_lines, k);
        curr_stats->bamcheck->fwd_del_above_baseline_pct = result->pct_above_baseline;
        curr_stats->bamcheck->fwd_del_below_baseline_pct = result->pct_below_baseline;
        free(result);

        result = bamcheck_cycles(curr_stats->del_cycles_2nd, ic_lines, k);
        curr_stats->bamcheck->rev_del_above_baseline_pct = result->pct_above_baseline;
        curr_stats->bamcheck->rev_del_below_baseline_pct = result->pct_below_baseline;
        free(result);
    }

}

bamcheck_bcd_t* bamcheck_base_content_baseline(uint64_t *base_prop, int n) {
    bamcheck_bcd_t *result = calloc(1, sizeof(bamcheck_bcd_t));

    bamcheck_baseline_delta* deviation;

    double mean = 0.0;
    int64_t i;
    double total = 0;
    for( i = 0; i < n; i++ ) {
        total += base_prop[i];
    }
    mean = total / (double)n;
    deviation = bamcheck_baseline_d(NULL, base_prop, 0, n, mean);

    result->mean_above_baseline = deviation->above_total / deviation->n;
    result->mean_below_baseline = deviation->below_total / deviation->n;
    result->max_above_baseline = deviation->above_max;
    result->max_below_baseline = deviation->below_max;

    result->total_mean_deviation = result->mean_above_baseline + result->mean_below_baseline;

    if ( result->max_above_baseline > result->max_below_baseline ) {
        result->max_baseline_deviation = result->max_above_baseline;
    }
    else {
        result->max_baseline_deviation = result->max_below_baseline;
    }

    free(deviation->delta);
    free(deviation);
    return result;
}

void bamcheck_base_content_deviation(stats_t *curr_stats){

    uint64_t gcc_a[curr_stats->max_len];
    uint64_t gcc_c[curr_stats->max_len];
    uint64_t gcc_g[curr_stats->max_len];
    uint64_t gcc_t[curr_stats->max_len];
    //uint64_t gcc_n[curr_stats->max_len];
    //uint64_t gcc_o[curr_stats->max_len];

    int ibase;
    for (ibase=0; ibase<curr_stats->max_len; ibase++) {
        acgtno_count_t *acgtno_count = &(curr_stats->acgtno_cycles[ibase]);
        uint64_t acgt_sum = acgtno_count->a + acgtno_count->c + acgtno_count->g + acgtno_count->t;

        gcc_a[ibase] = 100.*acgtno_count->a/acgt_sum;
        gcc_c[ibase] = 100.*acgtno_count->c/acgt_sum;
        gcc_g[ibase] = 100.*acgtno_count->g/acgt_sum;
        gcc_t[ibase] = 100.*acgtno_count->t/acgt_sum;
        //gcc_n[ibase] = 100.*acgtno_count->n/acgt_sum;
        //gcc_o[ibase] = 100.*acgtno_count->other/acgt_sum;
    }

    curr_stats->bamcheck->bcd_a = bamcheck_base_content_baseline(gcc_a, curr_stats->max_len);
    curr_stats->bamcheck->bcd_c = bamcheck_base_content_baseline(gcc_c, curr_stats->max_len);
    curr_stats->bamcheck->bcd_g = bamcheck_base_content_baseline(gcc_g, curr_stats->max_len);
    curr_stats->bamcheck->bcd_t = bamcheck_base_content_baseline(gcc_t, curr_stats->max_len);
}

void calculate_bamcheck(stats_t *curr_stats){
    bamcheck_indel_peaks(curr_stats);
    bamcheck_base_content_deviation(curr_stats);
    bamcheck_quality_dropoff(curr_stats);
}

bamcheck_stats_t* bamcheck_stats_init()
{
    bamcheck_stats_t* bamcheck = calloc(1, sizeof(bamcheck_stats_t));

    bamcheck->fwd_ins_above_baseline_pct = 0;
    bamcheck->fwd_ins_below_baseline_pct = 0;

    bamcheck->fwd_del_above_baseline_pct = 0;
    bamcheck->fwd_del_below_baseline_pct = 0;

    bamcheck->rev_ins_above_baseline_pct = 0;
    bamcheck->rev_ins_below_baseline_pct = 0;

    bamcheck->rev_del_above_baseline_pct = 0;
    bamcheck->rev_del_below_baseline_pct = 0;

    return bamcheck;
}

void bamcheck_quality_dropoff(stats_t *curr_stats){
}

