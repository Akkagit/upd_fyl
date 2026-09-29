
#ifndef GVIEWER_H
#define GVIEWER_H

/**
* @file gviewer.h
* @brief trajectory display core algorithm
* @verbatim
		History
		-1.0	LiShengyu	       2020-05-20 create the file.
  @endverbatim

* @author LiShengyu
* @email
* @version 1.0
* @date 2020-05-20
* @license
*/

#include <Eigen/Core>
#include <GLFW/glfw3.h>
#include <vector>
#include <utility>
#include <thread>
#include <mutex>
#include <math.h>
#include <map>
#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif
#include "gutils/ExportVIS.h"

using namespace std;
typedef vector<pair<Eigen::Matrix3d, Eigen::Vector3d>> Frames;
typedef vector<Eigen::Vector3d> VPointCloud;
typedef vector<Eigen::Vector3d> Trajectory;

namespace gvis
{
	class VIS_LIBRARY_EXPORT gviewer
	{
	public:
		gviewer();
		~gviewer();
	public:
		void Show();
		void ClearView();
		void SetFrames(const Frames &t);
		void SetFrames(const vector<Frames> &vt);
		void AddNewFrame(const pair<Eigen::Matrix3d, Eigen::Vector3d> &f);

		void SetPointCloud(const VPointCloud &pc);
		void SetPointCloud(const vector<VPointCloud> &vpc);
		void AddNewPoint(const Eigen::Vector3d &p);

		void SetTrajectory(const Trajectory &t);
		void SetTrajectory(const vector<Trajectory> &vt);
		void AddNewPos(const Eigen::Vector3d &p);
		void Hide();
	private:
		void Run();

	private:
		static vector<gviewer*> vptr;
		GLFWwindow* window;
		thread *t;
		mutex m_mutex;
		vector<Trajectory> mv_trajectory;
		vector<VPointCloud> mv_pointCloud;
		vector<Frames> mv_frames;
	};
}
#endif
