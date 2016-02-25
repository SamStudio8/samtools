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

//TODO(samstudio8) Should probably read doubles over uint64_t

#include "bamcheck.h"

#include <stdint.h>
#include <stdio.h>
#include <math.h>

int cmpfunc (const void * a, const void * b){
       return ( *(uint64_t*)a - *(uint64_t*)b );
}

int cmpfunc_f (const void * a, const void * b){
    if (*(double*)a > *(double*)b) return 1;
    else if (*(double*)a < *(double*)b) return -1;
}

uint64_t* copy_arr(uint64_t *source, int n, int new_n, int filter){
    uint64_t *dest;
    dest = calloc(new_n,sizeof(uint64_t));

    //TODO
    if (dest == NULL){
        printf("Memory not allocated.");
    }

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

double* copy_arr_double(double *source, int n, int new_n, int filter){
    double *dest;
    dest = calloc(new_n,sizeof(double));

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

double med_triplet_double(double a, double b, double c){
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

double* runmed_double(double *cycles_arr, int n, int k){
    double *count;
    double *baseline = calloc(n, sizeof(double));

    size_t i;
    for(i = 0; i < n; i++){
        baseline[i] = 0;
    }

    double smooth_start[k/2];
    double smooth_end[k/2];
    for(i = 0; i < k/2; i++){
        smooth_start[i] = 0;
        smooth_end[i] = 0;
    }

    // Copy counts for the known number of valid IC lines (n)
    count = copy_arr_double(cycles_arr, n, n, 0);

    // Medianize sliding windows of k, beginning at i=k/2 (the first valid full window)
    double *k_window;
    for(i = (k/2); i < n-(k/2); i++){
        k_window = copy_arr_double(&count[i-(k/2)], k, k, 0);
        qsort(k_window, k, sizeof(double), cmpfunc_f);
        baseline[i] = k_window[k/2];
        free(k_window);
    }

    // Keep ends ( k/2 elements on each end )
    memcpy(baseline, &count[0], sizeof(double) * (k/2));
    memcpy(&baseline[n-(k/2)], &count[n-(k/2)], sizeof(double) * (k/2));

    // Smooth 1st and n-1th element
    smooth_start[1] = med_triplet_double(baseline[0], baseline[1], baseline[2]);
    smooth_end[(k/2)-2] = med_triplet_double(baseline[n-1], baseline[n-2], baseline[n-3]);

    // Smooth values for first and last k/2 elements (where the windows were too
    //  large to medianize initially). Work inwards starting with element pair (2, n-2).
    for(i = 2; i < (k/2)+1; i ++){
        int j = 2 * i - 1;

        k_window = copy_arr_double(&baseline[0], j, j, 0);
        qsort(k_window, j, sizeof(double), cmpfunc_f);
        smooth_start[i-1] = k_window[j/2];
        free(k_window);

        k_window = copy_arr_double(&baseline[n-j], j, j, 0);
        qsort(k_window, j, sizeof(double), cmpfunc_f);
        smooth_end[(k/2)-i] = k_window[j/2];
        free(k_window);
    }

    // First and last element with Tukey Rule
    smooth_start[0] = med_triplet_double(baseline[0], smooth_start[1], 3 * smooth_start[1] - 2 * smooth_start[2]);
    smooth_end[(k/2)-1] = med_triplet_double(baseline[n-1], smooth_end[(k/2)-2], 3 * smooth_end[(k/2)-2] - 2 * smooth_end[(k/2)-3]);

    // Move smoothed values over baseline
    memcpy(baseline, &smooth_start[0], sizeof(double) * (k/2));
    memcpy(&baseline[n-(k/2)], &smooth_end[0], sizeof(double) * (k/2));

    free(count);
    return baseline;
}
//TODO Super function to also return mean and median (seeing as they are trivial)
//TODO(samstudio8) Implement a skip list? Speed doesn't seem too poor.
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

void init_bamcheck_baseline_delta(bamcheck_baseline_delta *delta, int n){
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
}

void bamcheck_baseline_d_double(bamcheck_baseline_delta *result, double *count, int count_n, double baseline){

    int i;
    double curr_delta, above_min, above_max, below_min, below_max;
    above_min = above_max = below_min = below_max = -1;

    for(i = 0; i < count_n; i++){
        curr_delta = (double)count[i] - baseline;
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
}

void bamcheck_baseline_d(bamcheck_baseline_delta *result, uint64_t *baseline, uint64_t *count, int baseline_n, int count_n, double scalar_baseline){

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
}

void bamcheck_cycles(bamcheck_cycles_t *result, uint64_t *cycles_arr, int n, int k){
    uint64_t *count;
    count = copy_arr(cycles_arr, n, n, 0);
    uint64_t *baseline;
    baseline = runmed(cycles_arr, n, k);

    // Calculate baseline and counts above/below
    bamcheck_baseline_delta *deviation = calloc(1, sizeof(bamcheck_baseline_delta));
    init_bamcheck_baseline_delta(result, n);
    bamcheck_baseline_d(deviation, baseline, count, n, n, 0);

    result->pct_above_baseline = (double)deviation->above_total / (double)deviation->total_count * 100;
    result->pct_below_baseline = (double)deviation->below_total / (double)deviation->total_count * 100;
    result->total_count = deviation->total_count;

    free(deviation->delta);
    free(deviation);
    free(count);
    free(baseline);
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

    bamcheck_cycles_t *result = calloc(1, sizeof(bamcheck_cycles_t));
    if (ic_lines > (k/2)){
        bamcheck_cycles(result, curr_stats->ins_cycles_1st, ic_lines, k);
        curr_stats->bamcheck->fwd_ins_above_baseline_pct = result->pct_above_baseline;
        curr_stats->bamcheck->fwd_ins_below_baseline_pct = result->pct_below_baseline;

        bamcheck_cycles(result, curr_stats->ins_cycles_2nd, ic_lines, k);
        curr_stats->bamcheck->rev_ins_above_baseline_pct = result->pct_above_baseline;
        curr_stats->bamcheck->rev_ins_below_baseline_pct = result->pct_below_baseline;

        bamcheck_cycles(result, curr_stats->del_cycles_1st, ic_lines, k);
        curr_stats->bamcheck->fwd_del_above_baseline_pct = result->pct_above_baseline;
        curr_stats->bamcheck->fwd_del_below_baseline_pct = result->pct_below_baseline;

        bamcheck_cycles(result, curr_stats->del_cycles_2nd, ic_lines, k);
        curr_stats->bamcheck->rev_del_above_baseline_pct = result->pct_above_baseline;
        curr_stats->bamcheck->rev_del_below_baseline_pct = result->pct_below_baseline;
    }
    else {
        // Not enough fragment qualities to look at baseline windows
        curr_stats->bamcheck->fwd_ins_above_baseline_pct = -1;
        curr_stats->bamcheck->fwd_ins_below_baseline_pct = -1;
        curr_stats->bamcheck->rev_ins_above_baseline_pct = -1;
        curr_stats->bamcheck->rev_ins_below_baseline_pct = -1;
        curr_stats->bamcheck->fwd_del_above_baseline_pct = -1;
        curr_stats->bamcheck->fwd_del_below_baseline_pct = -1;
        curr_stats->bamcheck->rev_del_above_baseline_pct = -1;
        curr_stats->bamcheck->rev_del_below_baseline_pct = -1;
    }
    free(result);

}

void bamcheck_base_content_baseline(bamcheck_bcd_t *result, double *base_prop, int n) {
    bamcheck_baseline_delta *deviation = calloc(1, sizeof(bamcheck_baseline_delta));
    init_bamcheck_baseline_delta(deviation, n);

    double mean = 0.0;
    int64_t i;
    double total = 0;
    for( i = 0; i < n; i++ ) {
        total += base_prop[i];
    }
    mean = total / (double)n;
    bamcheck_baseline_d_double(deviation, base_prop, n, mean);

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
}

void bamcheck_base_content_deviation(stats_t *curr_stats){

    double gcc_a[curr_stats->max_len];
    double gcc_c[curr_stats->max_len];
    double gcc_g[curr_stats->max_len];
    double gcc_t[curr_stats->max_len];

    int ibase;
    for (ibase=0; ibase < curr_stats->max_len; ibase++) {
        acgtno_count_t *acgtno_count = &(curr_stats->acgtno_cycles[ibase]);
        uint64_t acgt_sum = acgtno_count->a + acgtno_count->c + acgtno_count->g + acgtno_count->t;

        gcc_a[ibase] = 100.*acgtno_count->a/acgt_sum;
        gcc_c[ibase] = 100.*acgtno_count->c/acgt_sum;
        gcc_g[ibase] = 100.*acgtno_count->g/acgt_sum;
        gcc_t[ibase] = 100.*acgtno_count->t/acgt_sum;
    }

    bamcheck_base_content_baseline(curr_stats->bamcheck->bcd_a, gcc_a, curr_stats->max_len);
    bamcheck_base_content_baseline(curr_stats->bamcheck->bcd_c, gcc_c, curr_stats->max_len);
    bamcheck_base_content_baseline(curr_stats->bamcheck->bcd_g, gcc_g, curr_stats->max_len);
    bamcheck_base_content_baseline(curr_stats->bamcheck->bcd_t, gcc_t, curr_stats->max_len);
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

    bamcheck->bcd_a = calloc(1, sizeof(bamcheck_bcd_t));
    bamcheck->bcd_c = calloc(1, sizeof(bamcheck_bcd_t));
    bamcheck->bcd_g = calloc(1, sizeof(bamcheck_bcd_t));
    bamcheck->bcd_t = calloc(1, sizeof(bamcheck_bcd_t));

    bamcheck->fwd_dropoff = calloc(1, sizeof(bamcheck_quality_dropoff_t));
    bamcheck->rev_dropoff = calloc(1, sizeof(bamcheck_quality_dropoff_t));

    return bamcheck;
}

void bamcheck_stats_destroy(bamcheck_stats_t *bamcheck){
    free(bamcheck->bcd_a);
    free(bamcheck->bcd_c);
    free(bamcheck->bcd_g);
    free(bamcheck->bcd_t);

    free(bamcheck->fwd_dropoff);
    free(bamcheck->rev_dropoff);
}

typedef struct {
    uint64_t start_longest_trend_pos;
    uint64_t start_longest_trend_neg;
    uint64_t start_longest_trend_above_thres;
    uint64_t start_longest_trend_below_thres;

    uint64_t length_longest_trend_pos;
    uint64_t length_longest_trend_neg;
    uint64_t length_longest_trend_above_thres;
    uint64_t length_longest_trend_below_thres;
} bamcheck_cycles_trends_t;

bamcheck_cycles_trends_t* cycle_trends(double *cycles, uint64_t n, double cutoff){
    bamcheck_cycles_trends_t *result = calloc(1, sizeof(bamcheck_cycles_trends_t));

    uint64_t i;

    int64_t inc_best_start, inc_best_length, inc_curr_start, inc_curr_length;
    inc_best_start = inc_best_length = inc_curr_start = inc_curr_length = -1;

    int64_t abv_best_start, abv_best_length, abv_curr_start, abv_curr_length;
    abv_best_start = abv_best_length = abv_curr_start = abv_curr_length = -1;

    int64_t blw_best_start, blw_best_length, blw_curr_start, blw_curr_length;
    blw_best_start = blw_best_length = blw_curr_start = blw_curr_length = -1;

    int64_t dec_best_start, dec_best_length, dec_curr_start, dec_curr_length;
    dec_best_start = dec_best_length = dec_curr_start = dec_curr_length = -1;

    for (i = 0; i < n; i++) {
        if (cycles[i] >= cutoff) {
            if (abv_curr_start < 0){
                abv_curr_start = i;
                abv_curr_length = 1;
            }
            else {
                abv_curr_length++;
            }
        }
        else {
            if (abv_curr_length > abv_best_length){
                abv_best_start = abv_curr_start;
                abv_best_length = abv_curr_length;
            }
            abv_curr_length = -1;
            abv_curr_start = -1;
        }
        if (cycles[i] <= cutoff) {
            if (blw_curr_start < 0){
                blw_curr_start = i;
                blw_curr_length = 1;
            }
            else {
                blw_curr_length++;
            }
        }
        else {
            if (blw_curr_length > blw_best_length){
                blw_best_start = blw_curr_start;
                blw_best_length = blw_curr_length;
            }
            blw_curr_length = -1;
            blw_curr_start = -1;
        }

        if (i > 0) {
            if(cycles[i] >= cycles[i-1]){
                if(inc_curr_start < 0){
                    inc_curr_start = i-1;
                    inc_curr_length = 1;
                }
                else{
                    inc_curr_length++;
                }
            }
            else{
                if (inc_curr_length > inc_best_length){
                    inc_best_length = inc_curr_length;
                    inc_best_start = inc_curr_start;
                }
                inc_curr_length = -1;
                inc_curr_start = -1;
            }

            if(cycles[i] <= cycles[i-1]){
                if(dec_curr_start < 0){
                    dec_curr_start = i-1;
                    dec_curr_length = 1;
                }
                else{
                    dec_curr_length++;
                }
            }
            else{
                if (dec_curr_length > dec_best_length){
                    dec_best_length = dec_curr_length;
                    dec_best_start = dec_curr_start;
                }
                dec_curr_length = -1;
                dec_curr_start = -1;
            }
        }
    }

    //Contiguous cycle may have continued to final cycle
    if (abv_curr_length > abv_best_length) {
        abv_best_start = abv_curr_start;
        abv_best_length = abv_curr_length;
    }
    if (inc_curr_length > inc_best_length) {
        inc_best_start = inc_curr_start;
        inc_best_length = inc_curr_length;
    }
    if (blw_curr_length > blw_best_length) {
        blw_best_start = blw_curr_start;
        blw_best_length = blw_curr_length;
    }
    if (dec_curr_length > dec_best_length) {
        dec_best_start = dec_curr_start;
        dec_best_length = dec_curr_length;
    }

    //TODO Handle situations where best is never found...
    if (abv_best_length == -1){
        abv_best_start = 0;
        abv_best_length = 0;
    }
    else { abv_best_start++; }

    if (inc_best_length == -1){
        inc_best_start = 0;
        inc_best_length = 0;
    }
    else { inc_best_start++; }

    if (blw_best_length == -1){
        blw_best_start = 0;
        blw_best_length = 0;
    }
    else{ blw_best_start++; }

    if (dec_best_length == -1){
        dec_best_start = 0;
        dec_best_length = 0;
    }
    else { dec_best_start++; }

    result->start_longest_trend_pos = inc_best_start;
    result->start_longest_trend_neg = dec_best_start;
    result->start_longest_trend_above_thres = abv_best_start;
    result->start_longest_trend_below_thres = blw_best_start;

    result->length_longest_trend_pos = inc_best_length;
    result->length_longest_trend_neg = dec_best_length;
    result->length_longest_trend_above_thres = abv_best_length;
    result->length_longest_trend_below_thres = blw_best_length;

    return result;
}

void bamcheck_quality_dropoff(stats_t *curr_stats){
    bamcheck_quality_dropoff_executor(curr_stats->bamcheck->fwd_dropoff, curr_stats->quals_1st, curr_stats->max_len - 1, curr_stats->nquals, curr_stats->max_qual);
    bamcheck_quality_dropoff_executor(curr_stats->bamcheck->rev_dropoff, curr_stats->quals_2nd, curr_stats->max_len - 1, curr_stats->nquals, curr_stats->max_qual);
}

void summarise_cycles(bamcheck_cycles_summary_t *summary, uint64_t *cycle_counts, uint64_t cycles_n, uint64_t n_bin_quals, uint64_t max_quals_n){

    double cycle_quality[max_quals_n];
    int64_t ibase;
    int iqual;

    // Fragments by cycle...
    for (ibase = 0; ibase <= cycles_n; ibase++) {
        // ...cycle fragments by quality
        double curr;
        double sum = 0;
        uint64_t zeros = 0;
        for (iqual = 0; iqual < max_quals_n; iqual++) {
            curr = (double)cycle_counts[ ibase*n_bin_quals + iqual];
            cycle_quality[iqual] = curr;
            sum += curr;
            if (curr == 0){
                zeros++;
            }
        }

        //qsort(cycle_quality, curr_stats->max_qual, sizeof(double), cmpfunc_f);
        double sum_weight = 0;
        double sum_weighted_qual = 0;
        double cum_sum_weight[max_quals_n - zeros];
        double weighted_quals[max_quals_n - zeros];
        uint64_t cum_sum_index[max_quals_n - zeros];
        int j = 0;
        for(iqual = 0; iqual < max_quals_n; iqual++){
            sum_weighted_qual += cycle_quality[iqual] * iqual;
            sum_weight += cycle_quality[iqual];
            if (cycle_quality[iqual] != 0){
                cum_sum_index[j] = iqual;
                cum_sum_weight[j] = sum_weight;
                j++;
            }
        }
        j=0;
        for(iqual = 0; iqual < max_quals_n; iqual++){
            if (cycle_quality[iqual] != 0){
                //weighted_quals[j] = (cycle_quality[iqual] * iqual) / sum_weight;
                //weighted_quals[j] = (100.0/sum_weight) * (cum_sum_weight[iqual] - (cycle_quality[iqual] / 2.0));
                weighted_quals[j] = cum_sum_weight[iqual];
                j++;
            }
        }

        qsort(weighted_quals, max_quals_n - zeros, sizeof(double), cmpfunc_f);

        summary->cycle_means[ibase] = sum_weighted_qual / sum_weight;
        summary->cycle_iqrs[ibase] = cum_sum_index[wtd_percentile(cum_sum_weight, max_quals_n - zeros, sum_weight, 75)] - cum_sum_index[wtd_percentile(cum_sum_weight, max_quals_n - zeros, sum_weight, 25)];
        summary->cycle_medians[ibase] = cum_sum_index[wtd_percentile(cum_sum_weight, max_quals_n - zeros, sum_weight, 50)];
    }
}

void init_bamcheck_cycles_summary(bamcheck_cycles_summary_t *summary, uint64_t cycles_n){
    summary->cycle_means = calloc(cycles_n, sizeof(double));
    summary->cycle_medians = calloc(cycles_n, sizeof(double));
    summary->cycle_iqrs = calloc(cycles_n, sizeof(double));
}
void destroy_bamcheck_cycles_summary(bamcheck_cycles_summary_t *summary){
    free(summary->cycle_means);
    free(summary->cycle_medians);
    free(summary->cycle_iqrs);
}

void bamcheck_quality_dropoff_executor(bamcheck_quality_dropoff_t *result, uint64_t *cycle_counts, uint64_t cycles_n, uint64_t n_bin_quals, uint64_t max_quals_n){

    bamcheck_cycles_summary_t *summary = calloc(1, sizeof(bamcheck_cycles_summary_t));
    init_bamcheck_cycles_summary(summary, cycles_n);
    summarise_cycles(summary, cycle_counts, cycles_n, n_bin_quals, max_quals_n);

    int drop = 3;
    int k = 25;
    double iqr_cutoff = 10.0;

    bamcheck_cycles_trends_t *trends;

    trends = cycle_trends(&summary->cycle_iqrs[drop], cycles_n - (drop+1), iqr_cutoff);
    result->iqr_inc_contig_start = trends->start_longest_trend_above_thres;
    result->iqr_inc_contig_length = trends->length_longest_trend_above_thres;
    free(trends);

    double *mean_baseline;
    mean_baseline = runmed_double(summary->cycle_means, cycles_n, k);
    trends = cycle_trends(&mean_baseline[drop], cycles_n - (drop+1), 0);
    result->runmed_mean_dec_contig_start = trends->start_longest_trend_neg;
    result->runmed_mean_dec_contig_length = trends->length_longest_trend_neg;
    result->runmed_mean_dec_high = mean_baseline[trends->start_longest_trend_neg+drop-1];
    result->runmed_mean_dec_low = mean_baseline[trends->start_longest_trend_neg+(drop-2)+ trends->length_longest_trend_neg];
    result->runmed_mean_dec_range = result->runmed_mean_dec_high - result->runmed_mean_dec_low;
    free(trends);
    free(mean_baseline);

    destroy_bamcheck_cycles_summary(summary);
    free(summary);
}

double percentile(double *values, uint64_t n, int tile){
    //TODO(samstudio8) Potential optimisation by using partial sorting
    float k_index = (tile/(float)100) * n;
    if ( ceil(k_index) != k_index ) {
        k_index = ceil(k_index);
        return values[(int)k_index-1];
    }
    else {
        return ( (values[(int)k_index-1] + values[(int)k_index]) / (float) 2 );
    }
}

//TODO
uint64_t wtd_percentile(double *cumsums, uint64_t cumsums_n, uint64_t n, int tile){
    //TODO(samstudio8) Potential optimisation by using partial sorting
    float k_index = (tile/(float)100) * n;
    uint64_t target;
    if ( ceil(k_index) != k_index ) {
        k_index = ceil(k_index);
        target = (int)k_index-1;
    }
    else {
        target = (int)k_index;
    }

    //TODO Avoid cases where n is too small
    if( cumsums_n == 0 ){
        return -1;
    }

    size_t i;
    for (i = cumsums_n-1; i >= 0; i--){
        if (target > cumsums[i]){
            //TODO Getting there... Need to interpolate...
            return i+1;
        }
    }
}

