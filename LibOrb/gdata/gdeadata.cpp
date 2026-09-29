

#include"gdeadata.h"

using namespace std;

namespace great
{
	/** @brief constructor. */
	t_gdeadata::t_gdeadata() :t_gdata()
	{
		id_type(t_gdata::DESAISCOPOLECOEF);
	}


	/** @brief whether the ocean tide data is empty.
	* @return  bool
	*	@retval   true   ocean tide data is empty
	*   @retval   false  ocean tide data is existent
	*/
	bool t_gdeadata::is_empty()
	{
		if (_Anm.size() == 0 || _Bnm.size() == 0 || _AInm.size() == 0 || _BInm.size() == 0)
		{
			return true;
		}
		else
		{
			return false;
		}
	}


	/**
	* @brief get all oceantide coefficient.
	* @return   all oceantide coefficient
	*/
	map <int, map< int, double >>& t_gdeadata:: Anm()
	{ 
		return _Anm; 
	}

	/**
	* @brief get all oceantide coefficient.
	* @return   all oceantide coefficient
	*/
	map <int, map< int, double >>& t_gdeadata::Bnm()
	{
		return _Bnm;
	}

	/**
	* @brief get all oceantide coefficient.
	* @return   all oceantide coefficient
	*/
	map <int, map< int, double >>& t_gdeadata::AInm()
	{
		return _AInm;
	}

	/**
	* @brief get all oceantide coefficient.
	* @return   all oceantide coefficient
	*/
	map <int, map< int, double >>& t_gdeadata::BInm()
	{
		return _BInm;
	}

	/**
	* @brief get all load coefficient.
	* @return   all load coefficient
	*/
	vector<double>& t_gdeadata::load()
	{ 
		return _load; 
	}

	/**
		* @brief add one oceanide record to _oct.
		* @param[in]   m		index m
		* @param[in]   n		index n
		* @param[in]   A B AI BI	coefficient
	*/
	
	void t_gdeadata::set_rec(const int& n, const int& m, double A, double B, double AI, double BI)
	{
    #ifdef BMUTEX   
		boost::mutex::scoped_lock lock(_mutex);
    #endif 
		_gmutex.lock();
		
		_Anm[n][m] = A;
		_Bnm[n][m] = B;
		_AInm[n][m] = AI;
		_BInm[n][m] = BI;
		
		_gmutex.unlock();
		return;
	}

	/**
	* @brief add a load coefficient value to vector.
	* @param[in]   load		load coefficient value
	*/
	void t_gdeadata::set_load( double load)
	{
    #ifdef BMUTEX   
		boost::mutex::scoped_lock lock(_mutex);
    #endif 
		_gmutex.lock();
		_load.push_back(load);
		_gmutex.unlock();
		return;
	}

	
}

