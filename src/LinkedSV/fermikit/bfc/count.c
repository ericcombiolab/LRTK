#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bfc.h"
#include "bbf.h"
#include "htab.h"
#include "bseq.h"

/* A note on multi-threading

   The bloom filter is always the same regardless of how many threads in use.
   However, the k-mer inserted to the hash table may be slightly different.
   Suppose k-mers A and B are both singletons and that if A is inserted first,
   B is a false positive and gets inserted to the hash table. In the
   multi-threading mode, nonetheless, B may be inserted before A. In this case,
   B is not a false positive any more. This is not a bug. The k-mers put into
   the hash table depends on the order of input.
*/

#define CNT_BUF_SIZE 256

typedef struct { // cache to reduce locking
	uint64_t y[2];
	int is_high;
} insbuf_t;

typedef struct {
	const bfc_opt_t *opt;
	bseq_file_t *ks;
	bfc_bf_t *bf, *bf_high;
	bfc_ch_t *ch;
	int *n_buf;
	insbuf_t **buf;
} cnt_shared_t;

typedef struct {
	int n_seqs;
	bseq1_t *seqs;
	cnt_shared_t *cs;
} cnt_step_t;

static int bfc_kmer_bufclear(cnt_shared_t *cs, int forced, int tid)
{
	int i, k, r;
	if (cs->ch == 0) return 0;
	for (i = k = 0; i < cs->n_buf[tid]; ++i) {
		r = bfc_ch_insert(cs->ch, cs->buf[tid][i].y, cs->buf[tid][i].is_high, forced);
		if (r < 0) cs->buf[tid][k++] = cs->buf[tid][i];
	}
	cs->n_buf[tid] = k;
	return k;
}

static void bfc_kmer_insert(cnt_shared_t *cs, const bfc_kmer_t *x, int is_high, int tid)
{
	int k = cs->opt->k, ret;
	uint64_t y[2], hash;
	hash = bfc_kmer_hash(k, x->x, y);
	ret = bfc_bf_insert(cs->bf, hash);
	if (ret == cs->opt->n_hashes) {
		if (cs->ch && bfc_ch_insert(cs->ch, y, is_high, 0) < 0) { // counting with a hash table
			insbuf_t *p;
			if (bfc_kmer_bufclear(cs, 0, tid) == CNT_BUF_SIZE)
				bfc_kmer_bufclear(cs, 1, tid);
			p = &cs->buf[tid][cs->n_buf[tid]++];
			p->y[0] = y[0], p->y[1] = y[1], p->is_high = is_high;
		} else if (cs->bf_high) // keep high-occurrence k-mers
			bfc_bf_insert(cs->bf_high, hash);
	}
}

static void worker_count(void *_data, long k, int tid)
{
	cnt_step_t *data = (cnt_step_t*)_data;
	cnt_shared_t *cs = data->cs;
	bseq1_t *s = &data->seqs[k];
	const bfc_opt_t *o = cs->opt;
	int i, l;
	bfc_kmer_t x = bfc_kmer_null;
	uint64_t qmer = 0, mask = (1ULL<<o->k) - 1;
	for (i = l = 0; i < s->l_seq; ++i) {
		int c = seq_nt6_table[(uint8_t)s->seq[i]] - 1;
		if (c < 4) {
			bfc_kmer_append(o->k, x.x, c);
			qmer = (qmer<<1 | (s->qual == 0 || s->qual[i] - 33 >= o->q)) & mask;
			if (++l >= o->k) bfc_kmer_insert(cs, &x, (qmer == mask), tid);
		} else l = 0, qmer = 0, x = bfc_kmer_null;
	}
}

static void *bfc_count_cb(void *shared, int step, void *_data)
{
	cnt_shared_t *cs = (cnt_shared_t*)shared;
	if (step == 0) {
		cnt_step_t *ret;
		ret = calloc(1, sizeof(cnt_step_t));
		ret->seqs = bseq_read(cs->ks, cs->opt->chunk_size, 0, &ret->n_seqs);
		ret->cs = cs;
		fprintf(stderr, "[M::%s] read %d sequences\n", __func__, ret->n_seqs);
		if (ret->seqs) return ret;
		else free(ret);
	} else if (step == 1) {
		int i;
		double rt, eff;
		cnt_step_t *data = (cnt_step_t*)_data;
		kt_for(cs->opt->n_threads, worker_count, data, data->n_seqs);
		rt = realtime() - bfc_real_time;
		eff = 100. * cputime() / (rt + 1e-6);
		if (cs->ch) {
			fprintf(stderr, "[M::%s @%.1f*%.1f%%] processed %d sequences; # distinct k-mers: %ld\n",
					__func__, rt, eff, data->n_seqs, (long)bfc_ch_count(cs->ch));
		} else {
			fprintf(stderr, "[M::%s @%.1f*%.1f%%] processed %d sequences\n",
					__func__, rt, eff, data->n_seqs);
		}
		for (i = 0; i < cs->opt->n_threads; ++i)
			bfc_kmer_bufclear(cs, 1, i);
		for (i = 0; i < data->n_seqs; ++i) {
			bseq1_t *s = &data->seqs[i];
			free(s->seq); free(s->qual); free(s->comment); free(s->name);
		}
		free(data->seqs); free(data);
	}
	return 0;
}

void *bfc_count(const char *fn, const bfc_opt_t *opt)
{
	cnt_shared_t cs;
	void *ret;
	int i;

	memset(&cs, 0, sizeof(cnt_shared_t));
	cs.opt = opt;
	cs.bf = bfc_bf_init(opt->bf_shift, opt->n_hashes);
	if (!opt->filter_mode) {
		cs.ch = bfc_ch_init(opt->k, opt->l_pre);
		cs.n_buf = calloc(opt->n_threads, sizeof(int));
		cs.buf = calloc(opt->n_threads, sizeof(void*));
		for (i = 0; i < opt->n_threads; ++i)
			cs.buf[i] = malloc(CNT_BUF_SIZE * sizeof(insbuf_t));
		cs.ks = bseq_open(fn);
		kt_pipeline(opt->no_mt_io? 1 : 2, bfc_count_cb, &cs, 2);
		bseq_close(cs.ks);
		for (i = 0; i < opt->n_threads; ++i) free(cs.buf[i]);
		free(cs.buf); free(cs.n_buf);
		ret = (void*)cs.ch;
	} else {
		cs.bf_high = bfc_bf_init(opt->bf_shift, opt->n_hashes);
		cs.ks = bseq_open(fn);
		kt_pipeline(opt->no_mt_io? 1 : 2, bfc_count_cb, &cs, 2);
		bseq_close(cs.ks);
		ret = (void*)cs.bf_high;
	}
	bfc_bf_destroy(cs.bf);
	return ret;
}
