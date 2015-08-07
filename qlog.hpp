/**
 * @file qlog.hpp
 * @brief Thread-safe (and colorful!) functions to generate log.
 * @author Qing Wang
 * @date 2013-10-03, 2015-04-04
 */

#ifndef _QLOG_HPP_
#define _QLOG_HPP_

// Would you like to have colorful output?
#define _COLORFUL
// Do we need parallel writing to stderr?
#define _PARALLEL

#ifdef __cplusplus
#include <cstdio>
#include <cstdlib>
#include <cerrno>
#else
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#endif

#include <pthread.h>

#define likely(x) __builtin_expect(!!(x), 1) 
#define unlikely(x) __builtin_expect(!!(x), 0)

#ifdef _PARALLEL
static pthread_mutex_t _mutex_stderr = PTHREAD_MUTEX_INITIALIZER;
#ifdef __cplusplus
/**
 * An object to auto remove the mutex on exit.
 */
static class _class_auto_remover
{public:
	~_class_auto_remover() { pthread_mutex_destroy(&_mutex_stderr); }
} _auto_remover ;
#endif
#define LOCK() pthread_mutex_lock(&_mutex_stderr)
#define UNLOCK() pthread_mutex_unlock(&_mutex_stderr)
#else
#define LOCK()
#define UNLOCK()
#endif

// Foreground color
// [31m red
// [32m green
// [33m yellow
// [34m blue
// [35m magenta
// [36m cyan
#ifdef _COLORFUL
#define S_ERROR "\033[31m[ERROR]\033[m "
#define S_WARNING "\033[33m[WARNING]\033[m "
#define S_DEBUG "\033[36m[DEBUG]\033[m "
#define S_INFO "\033[32m[INFO]\033[m "
#define S_ASSERT "\033[35m[ASSERT]\033[m "
#else
#define S_ERROR "[ERROR] "
#define S_WARNING "[WARNING] "
#define S_DEBUG "[DEBUG] "
#define S_INFO "[INFO] "
#define S_ASSERT "[ASSERT] "
#endif

/**
 * Print where it goes wrong and exit.
 */
#define error(...) \
	do { LOCK(); \
		fprintf(stderr,"In %s(), %s:%d\n",__func__, __FILE__, __LINE__); \
		if(errno) fprintf(stderr,"[ERRNO:%02d] '%s'\n",errno,strerror(errno)); \
		fprintf(stderr, S_ERROR); \
		fprintf(stderr, __VA_ARGS__); \
		UNLOCK(); \
		exit(-1); \
	}while(0)

/**
 * Print where it goes wrong and continue.
 */
#define warning(...) \
	do { LOCK(); \
		fprintf(stderr,"In %s(), %s:%d\n",__func__, __FILE__, __LINE__); \
		fprintf(stderr, S_WARNING); \
		fprintf(stderr, __VA_ARGS__); \
		UNLOCK(); \
	}while(0)

/**
 * Print where it goes wrong and continue.
 */
#define debug(...) \
	do { LOCK(); \
		fprintf(stderr,"In %s(), %s:%d\n",__func__, __FILE__, __LINE__); \
		fprintf(stderr, S_DEBUG); \
		fprintf(stderr, __VA_ARGS__); \
		UNLOCK(); \
	}while(0)

/**
 * Print something.
 */
#define info(...) \
	do { LOCK(); \
		fprintf(stderr, S_INFO); \
		fprintf(stderr, __VA_ARGS__); \
		UNLOCK(); \
	}while(0)

/**
 * Assert, print where and exit on failure.
 */
#define qassert(x) \
	do { if(unlikely(!(x))) { LOCK(); \
		fprintf(stderr,"In %s(), %s:%d\n",__func__, __FILE__, __LINE__); \
		fprintf(stderr, S_ASSERT); \
		fprintf(stderr, "'%s' failed.\n",#x); \
		UNLOCK(); \
		exit(-1); \
	}}while(0)

#endif // _QLOG_HPP_
