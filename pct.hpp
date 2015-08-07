#pragma once

#include <cstdint>

/*
 * Pre-Computed Table for lgamma(i+delta)
 */
class PCT
{
public:
	uint32_t len;
	double (*func)(double);
	double delta;
	double * x;

	PCT() : len(0), func(NULL), delta(0), x(NULL) {}

	PCT(uint32_t _len, double (*_func)(double), double _delta)
		: len(0), func(NULL), delta(0), x(NULL)
	{
		make(_len,_func,_delta);
	}

	void dtor()
	{
		if(NULL!=x)
		{
			free(x);
			x = NULL;
		}
	}

	void make(uint32_t _len, double (*_func)(double), double _delta)
	{
		len = _len;
		func = _func;
		delta = _delta;
		if(NULL!=x)
			free(x);
		x = (double*)malloc(len*sizeof(double));
		for(uint32_t i=0;i<len;i++)
			x[i] = i + delta;
		for(uint32_t i=0;i<len;i++)
			x[i] = func(x[i]);
	}

	double operator()(uint32_t i)
	{
		if(i<len)
			return x[i];
		return func(i+delta);
	}
};

