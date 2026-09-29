#ifndef GCMD_H
#define GCMD_H

#include "gutils/ExportUtl.h"

namespace great
{
	class t_gcmd
	{
	public:
		t_gcmd();
		t_gcmd(string app_name, string app_path);
		virtual ~t_gcmd();

		bool operator-(const pair<string, string>& par);     // t_gcmd
		bool operator+(const pair<string, string>& par);     // t_gcmd
		t_gcmd operator=(const t_gcmd& cmd);

	protected:

		string _app_name;
		string _app_path;
		string _cmd_log;
		string _app_log;
		string _cmd_log_path;
		string _app_log_path;

		map<string, string> _app_pars;
		string _cmd;
	public:
		void set_log_path(string cmd_log_path, string app_log_path);
		bool add_cmd_pars(string opt, string pars);
		string run_cmd();
	private:
		
	};

}
#endif // !GCMD_H
