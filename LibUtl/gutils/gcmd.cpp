#include "gutils/gcmd.h"
#include "gutils/toml.hpp"

namespace great
{
	great::t_gcmd::t_gcmd()
	{

	}


	t_gcmd::t_gcmd(string app_name, string app_path)
	{
		_app_name = app_name;
		_app_path = app_path;
	}

	t_gcmd::~t_gcmd()
	{

	}

	bool t_gcmd::operator-(const pair<string, string>& par)
	{
		if (_app_pars.count(par.first) == 0) return true;

		_app_pars.erase(par.first);
		return true;
	}
	bool t_gcmd::operator+(const pair<string, string>& par)
	{
		if (_app_pars.count(par.first) > 0) return false;

		_app_pars[par.first] = par.second;
	}

	t_gcmd t_gcmd::operator=(const t_gcmd& cmd)
	{
		return t_gcmd();
	}

	void t_gcmd::set_log_path(string cmd_log_path, string app_log_path)
	{
		_cmd_log_path = cmd_log_path;
		_app_log_path = app_log_path;
	}

	bool t_gcmd::add_cmd_pars(string opt, string pars)
	{
		if (_app_pars.count(opt) > 0) return false;

		_app_pars[opt] = pars;
	}

	string t_gcmd::run_cmd()
	{
		_cmd.clear();
		for (const auto& it : _app_pars)
		{


		}

		return string("void");
	}
}

