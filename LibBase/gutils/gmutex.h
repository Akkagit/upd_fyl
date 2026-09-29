
#ifndef MUTEX_H
#define MUTEX_H

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  Purpose: mutual exclusion
  Version: $ Rev: $

  2013-08-14 /JD: created

-*/

#if defined _WIN32 || defined _WIN64
#include <windows.h>
#else
#include <pthread.h>
#endif

#include "gutils/ExportBase.h"
#include <omp.h>
#include <thread>
#include <mutex>

namespace gnut
{

	class BASE_LIBRARY_EXPORT t_gmutex
	{
	public:
		t_gmutex();
		t_gmutex(const t_gmutex& Other);
		~t_gmutex();

		t_gmutex operator=(const t_gmutex& Other);

		void lock();
		void unlock();


	protected:
#ifdef USE_OPENMP
		omp_lock_t _mutex;
#else
		mutex _mutex;

#endif // USE_OPENMP


	};

	/*
	template< class mutex_type >
	class t_glock
	{
	  public:
	   t_glock( mutex_type &mobj );
	  ~t_glock();

	  private:
		mutex_type *m_pMutex;
	};
	*/

} // namespace

#endif // MUTEX_H
