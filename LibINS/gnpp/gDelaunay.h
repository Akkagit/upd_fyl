/**
*
* @verbatim
History
-1.0  2020-12-08  BoWang transplant from RTPPP (Network Precise Position,NPP)
@endverbatim
* Copyright (c) 2020, Wuhan University. All rights reserved.
*
* @file		    gDelaunay.h
* @brief		Delaunay class to triangluate points set in 2D.
*               The procedure uses Double List for holding data, it can be optimized 
*               by using another data structure such as DAG, Quad-edge, etc.
*
* @author       BoWang Wuhan University
* @version		1.0.0
* @date		    2020-12-08
*
*/

#ifndef GDELAUNAY_H
#define GDELAUNAY_H

#include "gutils/ExportINS.h"
#include "gutils/gtriple.h"

using namespace gnut;

namespace great
{
	struct point_float
	{
		double x;
		double y;
		//char station[4];
		int indx;
	};

	struct  points3_ABC
	{
		point_float A;
		point_float B;
		point_float C;
	};

	// TYPES ////////////
	typedef struct VERTEX2D_TYP
	{
		double x;
		double y;
		string station;
		//char station[4];
		int indx;

	} VERTEX2D, * VERTEX2D_PTR;
	//顶点类型


	typedef struct EDGE_TYP
	{
		VERTEX2D v1;
		VERTEX2D v2;

	} EDGE, * EDGE_PTR;
	//角度类型

	typedef struct TRIANGLE_TYP
	{
		int i1; // vertex index
		int i2;
		int i3;

		TRIANGLE_TYP* pNext;
		TRIANGLE_TYP* pPrev;

	} TRIANGLE, * TRIANGLE_PTR;
	//三角形类型


	typedef struct MESH_TYP
	{
		int vertex_num;
		int triangle_num;

		VERTEX2D_PTR pVerArr; // point to outer vertices arrary
		TRIANGLE_PTR pTriArr; // point to outer triangles arrary
	} MESH, * MESH_PTR;
	//包含顶点类型，三角类型

	/**
	* @brief virtual class for delaunay method
	*/
	class INS_LIBRARY_EXPORT t_gDelaunay
	{
	public:
		t_gDelaunay();
		virtual ~t_gDelaunay();

		int getIntSite(t_gtriple& x_r, vector<t_gtriple>& posList, int outIndx[3]);

	private:

		void InitMesh(MESH_PTR pMesh, int ver_num);
		void UnInitMesh(MESH_PTR pMesh);

		void AddBoundingBox(MESH_PTR pMesh);
		void RemoveBoundingBox(MESH_PTR pMesh);
		double CounterClockWise(VERTEX2D_PTR pa, VERTEX2D_PTR pb, VERTEX2D_PTR pc);
		;
		void IncrementalDelaunay(MESH_PTR pMesh);

		void Insert(MESH_PTR pMesh, int ver_index);
		bool FlipTest(MESH_PTR pMesh, TRIANGLE_PTR pTestTri);

		double InCircle(VERTEX2D_PTR pa, VERTEX2D_PTR pb, VERTEX2D_PTR pp, VERTEX2D_PTR  pd);
		double InTriangle(MESH_PTR pMesh, VERTEX2D_PTR pVer, TRIANGLE_PTR pTri);

		void InsertInTriangle(MESH_PTR pMesh, TRIANGLE_PTR pTargetTri, int ver_index);
		void InsertOnEdge(MESH_PTR pMesh, TRIANGLE_PTR pTargetTri, int ver_index);

		// Helper functions
		void RemoveTriangleNode(MESH_PTR pMesh, TRIANGLE_PTR pTri);
		TRIANGLE_PTR AddTriangleNode(MESH_PTR pMesh, TRIANGLE_PTR pPrevTri, int i1, int i2, int i3);

		// I/O functions
		void Input(char* pFile, MESH_PTR pMesh);
		void Input_pos(vector<t_gtriple>& posList, MESH_PTR pMesh);
		points3_ABC Output(MESH_PTR pMesh, point_float Point);


		int IsInTriangle(point_float A, point_float B, point_float C, point_float D);
		double GetTriangleSquar(point_float pf0, point_float pf1, point_float pf2);
	};

}
#endif