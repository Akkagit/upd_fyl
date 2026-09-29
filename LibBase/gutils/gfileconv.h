#ifndef GFILECONV_H
#define GFILECONV_H

/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)

  Purpose: file conversion utilities
  Version: $ Rev: $

  2016-06-26 /JD: created

-*/

#include <string>
#include "gutils/ExportBase.h"

#if defined _WIN32 || defined _WIN32
#define PATH_SEPARATOR "\\"
#else
#define PATH_SEPARATOR "/"
#endif

#define GFILE_PREFIX "file://"

#define MAX_PATH_LEN 256

#ifdef WIN32 || _WIN32
#include <io.h>
#define ACCESS(fileName, accessMode) _access(fileName, accessMode)
#define MKDIR(path) _mkdir(path)
#else
#define ACCESS(fileName, accessMode) access(fileName, accessMode)
#define MKDIR(path) mkdir(path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)
#endif

using namespace std;

namespace gnut
{
	string base_name(const string& path);               // extract file base name
	string  dir_name(const string& path);               // extract file dir  name
	string file_name(const string& path);
	bool  dir_exists(const string& path);               // check existance of path
	int    make_path(const string& path);               // create path recursively
	int    make_dir(const string& path);                // create single directory
} // namespace

#endif