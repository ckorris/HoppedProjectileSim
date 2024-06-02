#pragma once
#include "MathStructs.h"
#include "Stats.h"
#include <vector>
#include <functional>


namespace hps
{
	using namespace std;

	extern "C" {
		__declspec(dllimport) bool Simulate(float startspeedmps, float distbetweensamples, bool applyphysics, float3 camPosOffset, float3 camRotOffset,
			float spinrpm, float bbDiameterMM, float bbMassGrams, float tempcelsius, float relhumidity01, float pressureHPa,
			float3 gravityVector, int2& collisionpoint, float& collisiondepth, float& totaltime, std::vector<float3>& linePoints,
			std::function<bool(const float3&, const float3&)> collisionFunc);
	}

	class Simulation
	{
	public:

		using CollisionDetectionFunc = std::function<bool(const float3&, const float3&)>;

		//static bool Simulate(sl::Mat depthmat, float speedmps, float distbetweensamples, bool applyphysics, sl::SensorsData sensordata,
		//	int2& collisionpoint, float& collisiondepth, float& totaltime, bool drawline, cv::Mat& drawlinetomat, cv::Scalar linecolor);

		static bool Simulate(float startspeedmps, float distbetweensamples, bool applyphysics, float3 camPosOffset, float3 camRotOffset,
			float spinrpm, float bbDiameterMM, float bbMassGrams, float tempcelsius, float relhumidity01, float pressureHPa,
			float bbToAirFrictionCoef, float3 gravityVector, CollisionDetectionFunc collisionDetectionFunc,
			float& collisiondepth, float& totaltime, vector<float3>& linePoints, Stats& stats);

		static float CalculateAirDensity(float totalairpressurehpa, float tempcelsius, float relhumidity01);

		static float CalculateAirViscosity(float tempcelsius);

		static float CalculateDragCoefficient(float spinrpm, float linearspeedmps, float bbdiametermeters,
			float airdensity, float airviscosity);

		static float CalculateLiftCoefficient(float spinrpm, float linearspeedmps, float bbdiametermeters);

		static float CalculateDragForce(float dragcoefficient, float airdensitykgm3, float surfaream2, float speedms);

		static float CalculateLiftForce(float liftcoefficient, float airdensitykgm3, float surfaream2, float speedms);

		static float CalculateMomentOfInertia(float bbmasskg, float bbdiameterm);

		static float CalculateSpinDragTorque(float airdensitykgm3, float airviscosity, float bbdiameterm, float spinrpm, float viscousfrictioncoefficient);

	private:
		static float3 CrossProduct(float3 v1, float3 v2);
		static float3 RotateVectorAroundAxis(float3 vectortorot, float3 axis, float angledegrees);
	};
}