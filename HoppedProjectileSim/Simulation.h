#pragma once
#include "MathStructs.h"
#include "Stats.h"
#include "PhysicsArgs.h"
#include <vector>
#include <functional>


namespace hps
{
	using namespace std;

#ifdef BUILD_DLL
#define DLL_API __declspec(dllexport)
#else
#define DLL_API __declspec(dllimport)
#endif

	extern "C" {

		typedef bool(__cdecl* CollisionDetectionFunc)(const float3&, const float3&);

		DLL_API bool RunSimulation(float distbetweensamples, bool applyphysics, float3 camPosOffset, float3 camRotOffset,
			PhysicsArgs physicsArgs, float3 gravityVector, CollisionDetectionFunc collisionDetectionFunc,
			float& collisiondepth, float& totaltime, float3** linePoints, int& linePointsCount, Stats& stats);
	}

	class Simulation
	{
	public:

		using CollisionDetectionFunc = std::function<bool(const float3&, const float3&)>;

		static bool Simulate(float distbetweensamples, bool applyphysics, float3 camPosOffset, float3 camRotOffset,
			PhysicsArgs physicsArgs, float3 gravityVector, CollisionDetectionFunc collisionDetectionFunc,
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