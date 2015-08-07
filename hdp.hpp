#pragma once

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "vec.hpp"
#include "pct.hpp"

#include "qlog.hpp"

// Proportional (for numerical stability) Exponential
template <typename T>
void prop_exp(T* x, uint32_t len)
{
	T mean = 0;
	for(uint32_t i=0;i<len;i++)
		if(likely(std::isfinite(x[i])))
			mean += x[i];
	mean /= len;
	for(uint32_t i=0;i<len;i++)
		x[i] = exp(x[i]-mean);
}

/**
 * HDP in CRF representation
 */
class HDP
{
public:
	const uint32_t n_doc;
	const uint32_t n_word;

	double alpha;
	double beta;
	double gamma;

	Vec<uint32_t>* dat;

	Vec<uint32_t>* table_stat;
	Vec<uint32_t> menu_stat;
	uint32_t menu_stat_sum;

	Vec<uint32_t*> word_stat;
	Vec<uint32_t> word_stat_sum;

	Vec<uint32_t>* table;
	Vec<uint32_t>* menu;

	PCT pct_log;
	PCT pct_log_b;
	PCT pct_log_nb;
	PCT pct_log_g;

	HDP(uint32_t _n_doc, uint32_t _n_word):
		n_doc(_n_doc), n_word(_n_word)
	{
		dat = (Vec<uint32_t>*)malloc(n_doc*sizeof(Vec<uint32_t>));
		memset(dat,0,n_doc*sizeof(Vec<uint32_t>));
		table_stat = (Vec<uint32_t>*)malloc(n_doc*sizeof(Vec<uint32_t>));
		memset(table_stat,0,n_doc*sizeof(Vec<uint32_t>));
		table = (Vec<uint32_t>*)malloc(n_doc*sizeof(Vec<uint32_t>));
		memset(table,0,n_doc*sizeof(Vec<uint32_t>));
		menu = (Vec<uint32_t>*)malloc(n_doc*sizeof(Vec<uint32_t>));
		memset(menu,0,n_doc*sizeof(Vec<uint32_t>));
		menu_stat_sum = 0;
	}

	~HDP() { dtor(); }

	void dtor()
	{
		for(uint32_t d=0;d<n_doc;d++)
			dat[d].dtor();
		free(dat);
		for(uint32_t d=0;d<n_doc;d++)
			table_stat[d].dtor();
		free(table_stat);
		for(uint32_t k=0;k<word_stat.len;k++)
			free(word_stat[k]);
		word_stat.clear();
		for(uint32_t d=0;d<n_doc;d++)
			table[d].dtor();
		free(table);	
		for(uint32_t d=0;d<n_doc;d++)
			menu[d].dtor();
		free(menu);
	}

	void config(double a, double b, double g, uint32_t buffer_size=65536*128)
	{
		alpha = a;
		beta = b;
		gamma = g;

		pct_log.make(buffer_size,log,0);
		pct_log_b.make(buffer_size,log,beta);
		pct_log_nb.make(buffer_size,log,beta*n_word);
		pct_log_g.make(buffer_size,log,gamma);
	}

	void add_entry(uint32_t doc, uint32_t word)
	{
		qassert(doc<n_doc);
		qassert(word<n_word);
		dat[doc].push_back(word);
	}

	void read_data(const char* filename) // in lda-c format
	{
		FILE * f = fopen(filename,"r");
		if(!f)
			error("Cannot open %s for reading.\n",filename);
		for(uint32_t d=0;d<n_doc;d++)
		{
			uint32_t n,w,m;
			qassert(1==fscanf(f,"%u ",&n));
			for(uint32_t i=0;i<n;i++)
			{
				qassert(2==fscanf(f,"%u:%u",&w,&m));
				for(uint32_t j=0;j<m;j++)
					add_entry(d,w);
			}
		}
	}

	void init0()
	{
		// 1 table for a doc, 1 menu for the franchise.
		word_stat.push_back((uint32_t*)malloc(n_word*sizeof(uint32_t)));
		memset(word_stat[0],0,n_word*sizeof(uint32_t));
		word_stat_sum.push_back(0);
		menu_stat.push_back(0);
		for(uint32_t d=0;d<n_doc;d++)
		{
			table_stat[d].push_back(0);
			menu[d].push_back(0); //menu[d][table[d][i]];
			menu_stat[0]++;
			menu_stat_sum++;
			for(uint32_t i=0;i<dat[d].len;i++)
			{
				table[d].push_back(0); // table[d][i]
				word_stat[0][dat[d][i]]++; //dat[d][i].n;
				word_stat_sum[0]++; //dat[d][i].n;
				table_stat[d][0]++; //dat[d][i].n;
			}
		}
	}

	void init()
	{
		static Vec<uint32_t> x(n_doc);
		if(x.len==0)
			for(uint32_t d=0;d<n_doc;d++)
				x.push_back(d);
		shuffle(&(x[0]),n_doc);
		for(uint32_t di=0;di<n_doc;di++)
			for(uint32_t i=0;i<dat[x[di]].len;i++)
			{
				uint32_t d = x[di];
				// assign_user(d,i)
				reassign_user(d,i,true);
			}
	}

	void gibbs_table() 
	{
		static Vec<uint32_t> x(n_doc);
		if(x.len==0)
			for(uint32_t d=0;d<n_doc;d++)
				x.push_back(d);
		shuffle(&(x[0]),n_doc);
		for(uint32_t d=0;d<n_doc;d++)
			for(uint32_t i=0;i<dat[x[d]].len;i++)
				reassign_user(x[d],i);
	}

	void reassign_user(uint32_t d, uint32_t i, bool firstrun=false)
	{
		static Vec<double> p; // prob without normalization
		static Vec<double> q; // cum prob to be filled by rmult()
		double lprob_t, lprob_k;
		uint32_t w = dat[d][i];
		// Remove statistics
		if(not firstrun) {
			uint32_t t = table[d][i];
			uint32_t k = menu[d][t];
			word_stat[k][w]--; //dat[d][i].n;
			word_stat_sum[k]--; //dat[d][i].n;
			table_stat[d][t]--; //dat[d][i].n;
		}
		// 1. Take an existing table.
		for(uint32_t t=0;t<table_stat[d].len;t++)
		{
			lprob_t = pct_log(table_stat[d][t]);
			p.push_back(lprob_t);
			q.push_back(0);
			// g_k = P(w_di | t_di)
			uint32_t k = menu[d][t];
			p[t] += pct_log_b(word_stat[k][w]) - pct_log_nb(word_stat_sum[k]);
		}
		// 2. Take a new table with an existing menu.
		lprob_t = log(alpha);
		for(uint32_t k=0;k<menu_stat.len;k++)
		{
			lprob_k = pct_log(menu_stat[k]) - pct_log_g(menu_stat_sum);
			p.push_back(lprob_t+lprob_k);
			q.push_back(0);
			p[p.len-1] += pct_log_b(word_stat[k][w]) - pct_log_nb(word_stat_sum[k]);
		}
		// 3. Take a new table with a new menu.
		lprob_k = pct_log_g(0) - pct_log_g(menu_stat_sum);
		p.push_back(lprob_t + lprob_k);
		q.push_back(0);
		p[p.len-1] += -pct_log(n_word);
		// Exponetial
		prop_exp(&p[0],p.len);
		// Draw random number
		uint32_t res = rmultinorm(&p[0],&q[0],p.len);
		p.clear();
		q.clear();
		uint32_t res_t = res<table_stat[d].len?res:table_stat[d].len;
		if(firstrun)
			table[d].push_back(0);
		table[d][i] = res_t;
		if(res_t==table_stat[d].len) // New table in d
		{
			uint32_t res_k = res - table_stat[d].len;
			table_stat[d].push_back(0);
			menu[d].push_back(res_k);
			if(res_k==menu_stat.len) // New menu!
			{
				word_stat.push_back((uint32_t*)malloc(n_word*sizeof(uint32_t)));
				memset(word_stat[word_stat.len-1],0,n_word*sizeof(uint32_t));
				word_stat_sum.push_back(0);
				menu_stat.push_back(0);
			}
			// Update for New table
			menu_stat[res_k]++;
			menu_stat_sum++;
		}
		// Update statistics
		if(true) {
			uint32_t t = table[d][i];
			uint32_t k = menu[d][t];
			word_stat[k][w]++;//dat[d][i].n;
			word_stat_sum[k]++;//dat[d][i].n;
			table_stat[d][t]++;//dat[d][i].n;
		}
	}

	void gibbs_menu()
	{
		static Vec<uint32_t> x(n_doc);
		if(x.len==0)
			for(uint32_t d=0;d<n_doc;d++)
				x.push_back(d);
		shuffle(&(x[0]),n_doc);
		for(uint32_t d=0;d<n_doc;d++)
			for(uint32_t t=0;t<menu[x[d]].len;t++)
				reassign_table(x[d],t);
	}
	
	void reassign_table(uint32_t d, uint32_t t)
	{
		static Vec<double> p;
		static Vec<double> q;
		if(true) //Remove statistics
		{
			uint32_t k = menu[d][t];
			for(uint32_t i=0;i<dat[d].len;i++)
				if(table[d][i]==t)
					word_stat[k][dat[d][i]]--;//dat[d][i].n;
			word_stat_sum[k] -= table_stat[d][t];
			menu_stat[k]--;
			menu_stat_sum--;
		}
		// 1. Take an existing menu
		for(uint32_t k=0;k<menu_stat.len;k++)
		{
			p.push_back(pct_log(menu_stat[k])); // possibly -inf
			q.push_back(0);
		}
		// 2. Take a new menu
		p.push_back(pct_log_g(0));
		q.push_back(0);
		// calculate g_k
		static Vec<uint32_t> local_stat(n_word);
		local_stat.len = n_word;
		memset(&(local_stat[0]),0,n_word*sizeof(uint32_t));
		uint32_t j = 0;
		for(uint32_t i=0;i<dat[d].len;i++) 
			if(table[d][i]==t)
			{
				for(uint32_t k=0;k<menu_stat.len;k++)
					p[k] += pct_log_b(word_stat[k][dat[d][i]]+local_stat[dat[d][i]]) - pct_log_nb(word_stat_sum[k]+j);
				p[menu_stat.len] += pct_log_b(local_stat[dat[d][i]]) - pct_log_nb(j);
				local_stat[dat[d][i]]++; // for duplicated word
				j++;
			}
		prop_exp(&(p[0]),p.len);
		// Draw random number
		uint32_t res = rmultinorm(&(p[0]),&(q[0]),p.len);
		p.clear();
		q.clear();
		menu[d][t] = res;
		if(res==menu_stat.len)
		{
			word_stat.push_back((uint32_t*)malloc(n_word*sizeof(uint32_t)));
			memset(&(word_stat[word_stat.len-1][0]),0,n_word*sizeof(uint32_t));
			word_stat_sum.push_back(0);
			menu_stat.push_back(0);
		}
		if(true) // Update statistics
		{
			menu_stat[res]++;
			menu_stat_sum++;
			for(uint32_t i=0;i<dat[d].len;i++)
				if(table[d][i]==t)
					word_stat[res][dat[d][i]]++; //=dat[d][i].n;
			word_stat_sum[res] += table_stat[d][t];
		}
	}

	void remove_empty()
	{
		// Remove empty tables
		for(uint32_t d=0;d<n_doc;d++)
		{
			for(int t=table_stat[d].len-1;t>=0;t--)
			{
				if(table_stat[d][t]>0)
					continue;
				// Remove table t's statistics
				menu_stat[menu[d][t]]--;
				menu_stat_sum--;
				// Move last table to table t
				uint32_t last = table_stat[d].len-1;
				for(uint32_t i=0;i<dat[d].len;i++)
					if(table[d][i]==last)
						table[d][i]=t;
				table_stat[d][t] = table_stat[d][last];
				table_stat[d].len--;
				menu[d][t] = menu[d][last];
				menu[d].len--;
			}
		}
		// Remove empty menu
		for(int k=menu_stat.len-1;k>=0;k--)
		{
			if(menu_stat[k]>0)
				continue;
			// Move last menu to menu k
			uint32_t last = menu_stat.len-1;
			for(uint32_t d=0;d<n_doc;d++)
				for(uint32_t t=0;t<menu[d].len;t++)
					if(menu[d][t]==last) // more efficient?
						menu[d][t]=k;
			menu_stat[k] = menu_stat[last];
			menu_stat.len--;
			word_stat_sum[k] = word_stat_sum[last];
			word_stat_sum.len--;
			free(word_stat[k]);
			word_stat[k] = word_stat[last];
			word_stat.len--;
		}
	}

	void summary(uint32_t verbosity)
	{
		if(verbosity>0)
			printf("#menu: %5u	#table: %8u\n",menu_stat.len,menu_stat_sum);
		if(verbosity>1)
		{
			for(uint32_t k=0;k<menu_stat.len;k++)
			{
				printf("menu %3u: ",k);
				printf("#table: %5u  ",menu_stat[k]);
				printf("#word: %7u\n",word_stat_sum[k]);
			}
		}
	}

	void output_topics(FILE* fo)
	{
		for(uint32_t k=0;k<word_stat.len;k++)
		{
			for(uint32_t w=0;w<n_word-1;w++)
				fprintf(fo,"%u\t",word_stat[k][w]);
			fprintf(fo,"%u\n",word_stat[k][n_word-1]);
		}
	}

	void output_assignments(FILE* fo)
	{
		uint32_t * cnt;
		qassert((cnt=(uint32_t*)malloc(menu_stat.len*sizeof(uint32_t))));
		for(uint32_t d=0;d<n_doc;d++)
		{
			memset(cnt,0,menu_stat.len*sizeof(uint32_t));
			for(uint32_t t=0;t<table_stat[d].len;t++)
				cnt[menu[d][t]] += table_stat[d][t];
			for(uint32_t k=0;k<menu_stat.len-1;k++)
				fprintf(fo,"%u\t",cnt[k]);
			fprintf(fo,"%u\n",cnt[menu_stat.len-1]);
		}
		free(cnt);
	}

	void check()
	{
		uint32_t * tmp;
		// Check num of menu
		qassert(word_stat.len==word_stat_sum.len);
		qassert(word_stat.len==menu_stat.len);
		// Check num of table
		for(uint32_t d=0;d<n_doc;d++)
			qassert(menu[d].len==table_stat[d].len);	
		debug("[PASS] N menu, N table\n");
		// Check word_stat_sum
		for(uint32_t k=0;k<word_stat.len;k++)
		{
			uint32_t s = 0;
			for(uint32_t w=0;w<n_word;w++)
				s += word_stat[k][w];
			qassert(word_stat_sum[k]==s);
		}
		for(uint32_t k=0;k<word_stat.len;k++)
		{
			uint32_t s=0;
			for(uint32_t d=0;d<n_doc;d++)
				for(uint32_t t=0;t<menu[d].len;t++)
					s += (menu[d][t]==k?table_stat[d][t]:0);
			qassert(s==word_stat_sum[k]);
		}
		debug("[PASS] word_stat_sum\n");
		// Check menu_stat_sum
		uint32_t s=0;
		for(uint32_t k=0;k<menu_stat.len;k++)
			s += menu_stat[k];
		qassert(s==menu_stat_sum);
		// Check table_stat
		for(uint32_t d=0;d<n_doc;d++)
		{
			tmp = (uint32_t*)malloc(table_stat[d].len*sizeof(uint32_t));
			memset(tmp,0,table_stat[d].len*sizeof(uint32_t));
			for(uint32_t i=0;i<dat[d].len;i++)
				tmp[table[d][i]]++; //=dat[d][i].n;
			for(uint32_t t=0;t<table_stat[d].len;t++)
				qassert(tmp[t]==table_stat[d][t]);
			free(tmp);
		}
		debug("[PASS] table_stat\n");
		// Check menu_stat
		tmp = (uint32_t*)malloc(menu_stat.len*sizeof(uint32_t));
		memset(tmp,0,menu_stat.len*sizeof(uint32_t));
		for(uint32_t d=0;d<n_doc;d++)
			for(uint32_t t=0;t<menu[d].len;t++)
				tmp[menu[d][t]]++;
		for(uint32_t k=0;k<menu_stat.len;k++)
			qassert(tmp[k]==menu_stat[k]);
		debug("[PASS] menu_stat\n");
		// Check word_stat
		tmp = (uint32_t*)malloc(word_stat.len*n_word*sizeof(uint32_t));
		memset(tmp,0,word_stat.len*n_word*sizeof(uint32_t));
		for(uint32_t d=0;d<n_doc;d++)
			for(uint32_t i=0;i<dat[d].len;i++)
			{
				uint32_t k=menu[d][table[d][i]];
				tmp[k*n_word+dat[d][i]]++; //=dat[d][i].n;
			}
		for(uint32_t k=0;k<word_stat.len;k++)
			for(uint32_t w=0;w<n_word;w++)
				qassert(tmp[k*n_word+w]==word_stat[k][w]);
		free(tmp);
		debug("[PASS] word_stat\n");
	}

};

