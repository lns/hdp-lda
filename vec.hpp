#pragma once

#include <cstdint>
#include <cstdlib>
#include "qlog.hpp"

template <typename T>
class Vec
{
public:
	uint32_t max_len;
	uint32_t len;
	T* head;

	Vec(): max_len(0), len(0), head(NULL) {};

	Vec(uint32_t _max_len): max_len(_max_len), len(0), head(NULL)
	{
		qassert(max_len>0);
		qassert((head=(T*)malloc(max_len*sizeof(T))));
	}

	~Vec() { dtor(); }

	void dtor()
	{
		if(NULL==head)
			return;
		free(head);
		head = NULL;
	}

	void push_back(const T& t)
	{
		if(len==max_len)
		{
			max_len = 4>2*max_len?4:2*max_len; // overflow
			head = (T*)realloc(head,max_len*sizeof(T));
		}
		head[len++] = t;
	}

	void clear()
	{
		len = 0;
	}

	T& operator[](uint32_t i)
	{
		qassert(i<len);
		return head[i];
	}
};

