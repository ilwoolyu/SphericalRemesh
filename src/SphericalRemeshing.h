/*************************************************
*	SphericalRemeshing.cpp
*
*	Release: Sep 2016
*	Update: June 2021
*
*	Ulsan National Institute of Science and Technology
*	Department of Computer Science and Engineering
*	
*	Ilwoo Lyu, ilwoolyu@unist.ac.kr
*************************************************/

#pragma once

#include <stdio.h>
#include <vector>
#include <algorithm>
#include "Mesh.h"
#include "AABB_Sphere.h"

using namespace std;

class SphericalRemeshing
{
public:
	SphericalRemeshing(void);
	~SphericalRemeshing(void);
	SphericalRemeshing(const char *subject, const char *sphere, const char *dfield, bool keepColor, const char *sphere_t = NULL, const char *colormap = NULL, vector<string> property = vector<string>(), bool interpolation = true, int deg = -1, bool verbose = true, const char *pbtype = NULL, bool backward = false);
	int degree(void);
	void deform(int degree);
	void saveDeformedSurface(const char *filename);
	void saveDeformedSphere(const char *filename);
	void saveDeformedProperty(const char *filename);
	void saveBary(const char *filename);
	
private:
	void reconsCoord(const float *v0, float *v1, float *Y, float *coeff, int degree, float *pole, float *tan1, float *tan2);
	float dataInterpolation(float *refMap, int index, float *coeff, Mesh *mesh);
	int dataInterpolation(vector<int *> refMap, int index, float *coeff, Mesh *mesh, int channel);
	float dataMedian(float *refMap, int index, Mesh *mesh);
	float dataNN(float *refMap, int index, float *coeff, Mesh *mesh);
	void deformSurface(void);
	void deformData(void);

private:
	int m_degree;
	float *m_coeff;
	float *m_x, *m_y, *m_z;
	float m_pole[2];
	float m_tan1[3], m_tan2[3];
	AABB_Sphere *m_tree;
	bool m_interpolation;
	bool m_backward;
	bool m_keepColor;
	bool m_verbose;
	Mesh *m_sphere, *m_sphere_subj, *m_subj, *m_remesh;
	vector<float *> m_refMap;
	vector<float *> m_deData;
	vector<string> m_property;
	vector<int *> m_color;
	vector<int *> m_color_base;
	const char *m_pbtype;
	const char *m_subject;
	float *m_Y;
	float *m_v0;
};


