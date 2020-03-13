/**
* This file is part of DSO.
*
* Copyright 2016 Technical University of Munich and Intel.
* Developed by Jakob Engel <engelj at in dot tum dot de>,
* for more information see <http://vision.in.tum.de/dso>.
* If you use this code, please cite the respective publications as
* listed on the above website.
*
* DSO is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* DSO is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with DSO. If not, see <http://www.gnu.org/licenses/>.
*/


#pragma once

#include "Eigen/Core"
#include "sophus/sim3.hpp"
#include "sophus/se3.hpp"


namespace dso
{

// CAMERA MODEL TO USE


#define SSEE(val,idx) (*(((float*)&val)+idx))


#define MAX_RES_PER_POINT 8
#define NUM_THREADS 6


#define todouble(x) (x).cast<double>()


using SE3R = Sophus::SE3d;
using Sim3 = Sophus::Sim3d;
using SO3 = Sophus::SO3d;



#define CPARS 4


using MatXX = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>;
using MatCC = Eigen::Matrix<double,CPARS,CPARS>;
#define MatToDynamic(x) MatXX(x)



using MatC10 = Eigen::Matrix<double,CPARS,10>;
using Mat1010 = Eigen::Matrix<double,10,10>;
using Mat1313 = Eigen::Matrix<double,13,13>;

using Mat810 = Eigen::Matrix<double,8,10>;
using Mat83 = Eigen::Matrix<double,8,3>;
using Mat66 = Eigen::Matrix<double,6,6>;
using Mat53 = Eigen::Matrix<double,5,3>;
using Mat43 = Eigen::Matrix<double,4,3>;
using Mat42 = Eigen::Matrix<double,4,2>;
using Mat33 = Eigen::Matrix<double,3,3>;
using Mat22 = Eigen::Matrix<double,2,2>;
using Mat8C = Eigen::Matrix<double,8,CPARS>;
using MatC8 = Eigen::Matrix<double,CPARS,8>;
using Mat8Cf = Eigen::Matrix<float,8,CPARS>;
using MatC8f = Eigen::Matrix<float,CPARS,8>;

using Mat88 = Eigen::Matrix<double,8,8>;
using Mat77 = Eigen::Matrix<double,7,7>;

using VecC = Eigen::Matrix<double,CPARS,1>;
using VecCf = Eigen::Matrix<float,CPARS,1>;
using Vec13 = Eigen::Matrix<double,13,1>;
using Vec10 = Eigen::Matrix<double,10,1>;
using Vec9 = Eigen::Matrix<double,9,1>;
using Vec8 = Eigen::Matrix<double,8,1>;
using Vec7 = Eigen::Matrix<double,7,1>;
using Vec6 = Eigen::Matrix<double,6,1>;
using Vec5 = Eigen::Matrix<double,5,1>;
using Vec4 = Eigen::Matrix<double,4,1>;
using Vec3 = Eigen::Matrix<double,3,1>;
using Vec2 = Eigen::Matrix<double,2,1>;
using VecX = Eigen::Matrix<double,Eigen::Dynamic,1>;

using Mat33f = Eigen::Matrix<float,3,3>;
using Mat103f = Eigen::Matrix<float,10,3>;
using Mat22f = Eigen::Matrix<float,2,2>;
using Vec3f = Eigen::Matrix<float,3,1>;
using Vec2f = Eigen::Matrix<float,2,1>;
using Vec6f = Eigen::Matrix<float,6,1>;



using Mat49 = Eigen::Matrix<double,4,9>;
using Mat89 = Eigen::Matrix<double,8,9>;

using Mat94 = Eigen::Matrix<double,9,4>;
using Mat98 = Eigen::Matrix<double,9,8>;

using Mat81 = Eigen::Matrix<double,8,1>;
using Mat18 = Eigen::Matrix<double,1,8>;
using Mat91 = Eigen::Matrix<double,9,1>;
using Mat19 = Eigen::Matrix<double,1,9>;


using Mat84 = Eigen::Matrix<double,8,4>;
using Mat48 = Eigen::Matrix<double,4,8>;
using Mat44 = Eigen::Matrix<double,4,4>;


using VecNRf = Eigen::Matrix<float,MAX_RES_PER_POINT,1>;
using Vec12f = Eigen::Matrix<float,12,1>;
using Mat18f = Eigen::Matrix<float,1,8>;
using Mat66f = Eigen::Matrix<float,6,6>;
using Mat88f = Eigen::Matrix<float,8,8>;
using Mat84f = Eigen::Matrix<float,8,4>;
using Vec8f = Eigen::Matrix<float,8,1>;
using Vec10f = Eigen::Matrix<float,10,1>;
using Mat66f = Eigen::Matrix<float,6,6>;
using Vec4f = Eigen::Matrix<float,4,1>;
using Mat44f = Eigen::Matrix<float,4,4>;
using Mat1212f = Eigen::Matrix<float,12,12>;
using Vec12f = Eigen::Matrix<float,12,1>;
using Mat1313f = Eigen::Matrix<float,13,13>;
using Mat1010f = Eigen::Matrix<float,10,10>;
using Vec13f = Eigen::Matrix<float,13,1>;
using Mat99f = Eigen::Matrix<float,9,9>;
using Vec9f = Eigen::Matrix<float,9,1>;

using Mat42f = Eigen::Matrix<float,4,2>;
using Mat62f = Eigen::Matrix<float,6,2>;
using Mat12f = Eigen::Matrix<float,1,2>;

using VecXf = Eigen::Matrix<float,Eigen::Dynamic,1>;
using MatXXf = Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic>;


using MatPCPC = Eigen::Matrix<double,8+CPARS+1,8+CPARS+1>;
using MatPCPCf = Eigen::Matrix<float,8+CPARS+1,8+CPARS+1>;
using VecPC = Eigen::Matrix<double,8+CPARS+1,1>;
using VecPCf = Eigen::Matrix<float,8+CPARS+1,1>;

using Mat1414f = Eigen::Matrix<float,14,14>;
using Vec14f = Eigen::Matrix<float,14,1>;
using Mat1414 = Eigen::Matrix<double,14,14>;
using Vec14 = Eigen::Matrix<double,14,1>;






// transforms points from one frame to another.
struct AffLight
{
	AffLight(double a_, double b_) : a(a_), b(b_) {};
	AffLight() : a(0), b(0) {};

	// Affine Parameters:
	double a,b;	// I_frame = exp(a)*I_global + b. // I_global = exp(-a)*(I_frame - b).

	static Vec2 fromToVecExposure(float exposureF, float exposureT, AffLight g2F, AffLight g2T)
	{
		if(exposureF==0 || exposureT==0)
		{
			exposureT = exposureF = 1;
			//printf("got exposure value of 0! please choose the correct model.\n");
			//assert(setting_brightnessTransferFunc < 2);
		}

		double a = exp(g2T.a-g2F.a) * exposureT / exposureF;
		double b = g2T.b - a*g2F.b;
		return Vec2(a,b);
	}

	Vec2 vec()
	{
		return Vec2(a,b);
	}
};

}
