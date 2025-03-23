#include <math.h>
#include "Simulation.h"
#include "BasicMat.h"
#include "PhysicsArgs.h"
#include "SampleStats.h"

using namespace std;

namespace hps
{

	const float GRAVITY_ACCEL = 9.81;
	const float MAX_DISTANCE = 60.0;

	const float GAS_CONSTANT_DRY_AIR = 287.058; //J / (kg*K).
	const float GAS_CONSTANT_WATER_VAPOR = 461.495; //J / (kg*K).

	const float PI = 3.14159265358979323846f;


	DLL_API bool RunSimulation(float sampleTime, int maxSamples, float3 camPosOffset, float3 camRotOffset, 
		PhysicsArgs physicsArgs, float3 gravityVector, CollisionDetectionFunc collisionDetectionFunc, 
		float& collisiondepth, float& totaltime, float3** linePoints, int& linePointsCount, SampleStats** sampleStats,
		SimStats& stats)
	{
		std::vector<float3> linePointsVec;
		std::vector<SampleStats> sampleStatsVec;

		bool result = Simulation::Simulate(sampleTime, maxSamples, camPosOffset, camRotOffset,
			physicsArgs, gravityVector, collisionDetectionFunc, collisiondepth, totaltime, linePointsVec, sampleStatsVec, stats);

		//Update linePointsCount with the actual number of points.
		linePointsCount = static_cast<int>(linePointsVec.size());

		//Allocate memory for linePoints on the C++ side.
		*linePoints = new float3[linePointsVec.size()];
		*sampleStats = new SampleStats[linePointsVec.size()];

		//Copy data to the allocated array.
		std::copy(linePointsVec.begin(), linePointsVec.end(), *linePoints);
		std::copy(sampleStatsVec.begin(), sampleStatsVec.end(), *sampleStats);

		return result;
	}

	bool Simulation::Simulate(float sampleTime, int maxSamples, float3 camPosOffset, float3 camRotOffset, PhysicsArgs 
		physicsArgs, float3 gravityVector, CollisionDetectionFunc collisionDetectionFunc, float& collisiondepth, float& totaltime, 
		vector<float3>& linePoints, vector<SampleStats>& sampleStats, SimStats& stats)
	{
		float downspeed = 0;
		totaltime = 0;

		float3 lastvalidpoint = camPosOffset;
		float3 currentpoint = camPosOffset;

		float3 forwardnormal(0, 0, 1); //TODO: Account for rotation using Config

		//I'm handling the rotations without full-on 3D math because it's relative to the forward direction of the barrel,
		//where Z rot does not (yet) matter (because the IMU should account for that when we add gravity). 
		//Y rot first. 
		float3 yrotnormal = forwardnormal;
		float yanglerad = camRotOffset.y * PI / 180.0;
		yrotnormal.z = forwardnormal.z * cos(yanglerad) - forwardnormal.x * sin(yanglerad);
		yrotnormal.x = forwardnormal.z * sin(yanglerad) + forwardnormal.x * cos(yanglerad);

		float3 finalrotnormal = yrotnormal;
		float xanglerad = camRotOffset.x * PI / 180.0;
		finalrotnormal.z = yrotnormal.z * cos(xanglerad) - yrotnormal.y * sin(xanglerad);
		finalrotnormal.y = yrotnormal.z * sin(xanglerad) - yrotnormal.y * cos(xanglerad);


		//Calculate the gun upwards angle, which is the vector pointing upward from the gun.
		//This is used to calculate hop up, since the backspin of the bb will cause it to rise in this direction. 

		float3 vectorup(0, 1, 0); //Could just put 0 and 1 in the math directly, but this makes code easier to read. 
		float3 gunupnormal = vectorup;
		float zanglerad = camRotOffset.z * PI / 180.0;
		gunupnormal.x = vectorup.x * cos(zanglerad) - vectorup.y * sin(zanglerad);
		gunupnormal.y = vectorup.x * sin(zanglerad) + vectorup.y * cos(zanglerad);

		float3 velocity = finalrotnormal * physicsArgs.StartSpeedMPS; //Startingvelocity
		//float downspeedaddpersample = GRAVITY_ACCEL * timebetweendots * timebetweendots; //GRAVITY_ACCEL * timebetweendots is downward accel added per dot but STILL IN MPS. Multiply again to get how much it should change. 

		float currentSpinRPM = physicsArgs.SpinRPM;

		//Numbers needed for drag and hop-up/Magnus that don't change sample to sample. 
		float bbdiameterm = physicsArgs.BBDiameterMM / 1000;
		float crosssectionalarea = PI * pow(bbdiameterm * 0.5, 2);
		float airdensity = Simulation::CalculateAirDensity(physicsArgs.PressureHPa, physicsArgs.TemperatureCelsius, physicsArgs.RelativeHumidity01);
		float airviscosity = Simulation::CalculateAirViscosity(physicsArgs.TemperatureCelsius);
		float bbmasskg = physicsArgs.BBMassGrams / 1000.0;
		float momentofinertia = Simulation::CalculateMomentOfInertia(bbmasskg, bbdiameterm);
		//Update stats with air values, which won't change between samples. 
		stats.AirPressureHPa = physicsArgs.PressureHPa;
		stats.AirDensityKGM3 = airdensity;
		stats.AirViscosityPaS = airviscosity;

		int samplecount = 0; //Mostly for debug but also used here and there. 
		while (currentpoint.z < MAX_DISTANCE && samplecount < maxSamples) //Sliiiightly concerned this will be infinite. 
		{
			//Catch potential infinite loop. 
			if (velocity.z <= 0)
			{
				//cout << "Velocity z factor is less than or equal to zero. Aborting simulation." << endl;
				break;
			}

			//Calculate the time that passes in order to travel the distance per sample. 
			float speed = sqrt(pow(velocity.x, 2) + pow(velocity.y, 2) + pow(velocity.z, 2));

			//Declare physics values outside if statement so we can use them to update stats later if needed.
			//This does rely us on checking applyphysics again at that point to avoid passing inappropriate zeroes. 
			float dragcoefficient = 0;
			float dragforcenewtons = 0;
			float liftcoefficient = 0;
			float liftforcenewtons = 0;


			float3 velocitynorm = velocity / speed;

			//Gravity. Simplest first. 
			float3 changeDueToGravity = gravityVector * GRAVITY_ACCEL * sampleTime;

			SampleStats sampleStat;

			//Drag. 
			//float dragcoefficient = CalculateDragCoefficient(spinrpm, speed / timebetweendots, bbdiameterm, airdensity, airviscosity); //Using per-sample speed. 
			dragcoefficient = Simulation::CalculateDragCoefficient(currentSpinRPM, speed, bbdiameterm, airdensity, airviscosity, &sampleStat.ReynoldsDragCoef);
			dragforcenewtons = Simulation::CalculateDragForce(dragcoefficient, airdensity, crosssectionalarea, speed);
			//cout <<"Air Density: " << airdensity <<  " Newtons: " << dragforcenewtons << endl;

			float dragspeedchange = dragforcenewtons * sampleTime / bbmasskg;

			float3 changeDueToDrag = velocitynorm * -dragspeedchange;


			float3 changeDueToBackspin;


			//Hop-up/Magnus.
			if (currentSpinRPM != 0) //If it's not spinning, don't spend the cycles doing all this. 
			{
				//Calculate hop-up normal. It's orthogonal to the velocity. 
				//We get that by rotating around the cross product of the gun up normal and velocity by 90 degrees. 
				//Recall that the gun up normal is the upward direction of the gun, ie the direction the backspin will push the bb. 
				float3 hopupcross = Simulation::CrossProduct(velocity, gunupnormal);
				//Normalize it. 
				float hopupcrossmagnitude = sqrt(pow(hopupcross.x, 2) + pow(hopupcross.y, 2) + pow(hopupcross.z, 2));
				hopupcross = hopupcross / hopupcrossmagnitude;
				sampleStat.HopUpCross = hopupcross;

				//TODO: The nor

				float3 hopupnormal = Simulation::RotateVectorAroundAxis(velocitynorm, hopupcross, 90);
				sampleStat.HopUpNormal = hopupnormal;

				liftcoefficient = Simulation::CalculateLiftCoefficient(currentSpinRPM, speed, bbdiameterm); //If change to speed works, apply to drag. 

				//TEST
				liftcoefficient = -liftcoefficient;

				liftforcenewtons = Simulation::CalculateLiftForce(liftcoefficient, airdensity, crosssectionalarea, speed); //Force across a second. 
				//That force is how much force will occur over a second. 
				float liftspeedchange = liftforcenewtons * sampleTime / bbmasskg; //Speed change. 
				changeDueToBackspin = hopupnormal * liftspeedchange;

				//float spindragtorque = CalculateSpinDragTorque(speed, airdensity, airviscosity, bbdiameterm, currentSpinRPM,
				//	physicsArgs.BBToAirFrictionCoefficient, &sampleStat.ReynoldsSpinDragTorque);
				float spindragtorque = CalculateSpinDragTorqueSimple(airdensity, bbdiameterm / 2, currentSpinRPM);
				
				float dragangularaccel = spindragtorque / momentofinertia; //This is in radians per second.
				float sampledecayradpersec = dragangularaccel * sampleTime;
				float sampledecayRPM = sampledecayradpersec / (2 * PI) * 60;
				currentSpinRPM -= (int)round(sampledecayRPM);

				sampleStat.SpinDragTorque = spindragtorque;
				sampleStat.AngularDragAccel = spindragtorque;
				sampleStat.MomentOfInertia = momentofinertia;
				sampleStat.SampleDecayRPM = sampledecayRPM;

				if (samplecount < 1)
				{
					cout << "Sample " << samplecount << " torque: " << spindragtorque << ": dragangularaccel = " << dragangularaccel <<
						"Sample time: " << sampleTime << ". Decay this sample: " << sampledecayradpersec << "rad/sec. RPM change = " << sampledecayRPM << endl;
				}


				//if (currentSpinRPM < 0) currentSpinRPM = 0;
			}
			else
			{
				changeDueToBackspin = float3(0, 0, 0);
			}

			if (samplecount == 0)
			{
				stats.DragCoefficientBarrel = dragcoefficient;
				stats.DragForceBarrel = dragforcenewtons;
				stats.MagnusCoefficientBarrel = liftcoefficient;
				stats.MagnusForceBarrel = liftforcenewtons;
			}

			//Add the changes in velocity.
			velocity += changeDueToGravity;
			velocity += changeDueToDrag;
			velocity += changeDueToBackspin;

			currentpoint += velocity * sampleTime;

			float pointdepth = currentpoint.z; //Shorthand. 
			totaltime += sampleTime;

			//Add the position to the vector.

			linePoints.push_back(currentpoint);

			//Add stats.

			sampleStat.Velocity = velocity;
			sampleStat.ChangeDueToGravity = changeDueToGravity;
			sampleStat.ChangeDueToDrag = changeDueToDrag;
			sampleStat.ChangeDueToBackspin = changeDueToBackspin;
			sampleStat.CurrentSpinRPM = currentSpinRPM;
			sampleStat.DragCoefficient = dragcoefficient;
			sampleStat.DragForceNewtons = dragforcenewtons;
			sampleStat.LiftCoefficient = liftcoefficient;
			sampleStat.LiftForceNewtons = liftforcenewtons;

			sampleStats.push_back(sampleStat);

			bool hit = collisionDetectionFunc(lastvalidpoint, currentpoint);

			if (hit)
			{
				collisiondepth = lastvalidpoint.z;

				//Update stats for the impact point. 
				stats.SpeedImpact = speed;

				stats.DragCoefficientImpact = dragcoefficient;
				stats.DragForceImpact = dragforcenewtons;
				stats.MagnusCoefficientImpact = liftcoefficient;
				stats.MagnusForceImpact = liftforcenewtons;

				return true;
			}

			lastvalidpoint = currentpoint;
			samplecount++;
		}

		collisiondepth = lastvalidpoint.z;
		stats.DragCoefficientImpact = 0;
		stats.DragForceImpact = 0;
		stats.MagnusCoefficientImpact = 0;
		stats.MagnusForceImpact = 0;

		return false; //We didn't hit anything. 
	}



	//Returns the air density in kg/m^3.
	float Simulation::CalculateAirDensity(float totalairpressurehpa, float tempcelsius, float relhumidity01)
	{
		//Calculated with: https://www.omnicalculator.com/physics/air-density#how-to-calculate-the-air-density
		float tempkelvins = tempcelsius + 273.15;
		float saturationvaporpressure = 6.1078 * pow(10, 7.5 * tempcelsius / (tempcelsius + 237.3));
		float vaporpressure = saturationvaporpressure * relhumidity01;
		float dryairpressure = totalairpressurehpa - vaporpressure;

		float dryairpressure_pa = dryairpressure * 100; //Hectopascals to Pascals.
		float vaporpressure_pa = vaporpressure * 100;

		//ρ = (pd / (Rd * T)) + (pv / (Rv * T))

		float airdensity = (dryairpressure_pa / (GAS_CONSTANT_DRY_AIR * tempkelvins)) + (vaporpressure_pa / (GAS_CONSTANT_WATER_VAPOR * tempkelvins));
		return airdensity;
	}

	float Simulation::CalculateAirViscosity(float tempcelsius) //Pascal seconds
	{
		float tempkelvins = tempcelsius + 273.15;

		//We'll use Sutherland's formula for calculating this, which requires a reference temperature and viscosity. 
		//We'll use 0°C (273.15°L) where viscosity is 0.00001716.
		//Note we use the number 383.55, which is reference number 273.15 + the Sutherland temp of 110.4°K.
		//Used as reference: https://www.cfd-online.com/Wiki/Sutherland's_law and https://www.lmnoeng.com/Flow/GasViscosity.php

		//return 0.00001716 * pow(tempkelvins / 273.15, 1.5) * (383.55 / (tempkelvins + 110.4));
		//float viscositycp = (0.00001458 * pow(tempkelvins, 1.5)) / (tempkelvins + 110.4); //Poise. 
		//return viscositycp * 0.1; //To Pascal seconds. 

		double temperatureKelvin = tempcelsius + 273.15; // Convert temperature to Kelvin
		double referenceViscosity = 1.81e-5; // Reference viscosity at 273.15 K, in Pa*s
		double referenceTemperature = 273.15; // Reference temperature, in Kelvin
		double sutherlandsConstant = 120.0; // Sutherland's constant, in Kelvin

		return referenceViscosity * ((referenceTemperature + sutherlandsConstant) / (temperatureKelvin + sutherlandsConstant))
			* pow(temperatureKelvin / referenceTemperature, 1.5);
	}


	float Simulation::CalculateDragCoefficient(float spinrpm, float linearspeedmps, float bbdiametermeters,
		float airdensity, float airviscosity, float* reynolds) //Note speed in mps, not a smaller fraction. 
	{
		//This equation is made to fit experimental data, from Dr. Dyrkacz's "The Physics of Paintball" which he based on 
		//Achenbach, E., J. Fluid Mech. 54,565 (1972). See his page here:
		//https://web.archive.org/web/20040617080904/http://home.comcast.net/~dyrgcmn/pball/pballIntro.html

		//Calculate the Reynolds number, which we'll need in a few places. 
		//double reynolds = bbdiametermeters * airdensity * linearspeedmps * airviscosity;
		//double reynolds = bbdiametermeters * airdensity * linearspeedmps / airviscosity;
		float angvelocityradpersec = spinrpm * 2 * PI / 60.0;
		//double reynolds = pow(bbdiametermeters / 2, 2) * angvelocityradpersec * airdensity * airviscosity;
		//*reynolds = pow(bbdiametermeters / 2, 2) * linearspeedmps * airdensity * airviscosity;
		*reynolds = bbdiametermeters * airdensity * linearspeedmps / airviscosity;

		//cout << "Reynolds: " << reynolds << endl;

		//First calculate without spin. We're gonna start using doubles from now on. 
		double spinlessCD = (0.4274794 + 0.000001146254 * *reynolds - 7.559635 * pow(10, -12) * pow(*reynolds, 2) - 3.817309 * pow(10, -18) *
			pow(*reynolds, 3) + 2.389417 * pow(10, -23) * pow(*reynolds, 4)) / (1 - 0.000002120623 * *reynolds + 2.952772 * pow(10, -11) * pow(*reynolds, 2) -
				1.914687 * pow(10, -16) * pow(*reynolds, 3) + 3.125996 * pow(10, -22) * pow(*reynolds, 4));

		//If we have any spin, we then take that and apply a new formula. 
		if (spinrpm == 0) return (float)spinlessCD; //If no spin, we can quit here. 

		//Calculate the surface speed, which to my understanding is the speed at which a point on the outside of the circle is moving. 
		float surfacespeed = bbdiametermeters * PI * spinrpm / 60;

		//Calculate the ratio of spin velocity to linear velocity. Those two are represented as V and U mathematically, hence calling it vu.
		float vu = surfacespeed / linearspeedmps;

		//Aaaaand another crazy formula created via lots of experimentation by someone way better at the mathses than yours truly. 
		double spinCD = (spinlessCD + 2.2132291 * vu - 10.345178 * pow(vu, 2) + 16.15703 * pow(vu, 3) + -5.27306480 * pow(vu, 4)) /
			(1 + 3.1077276 * vu - 13.6598678 * pow(vu, 2) + 24.00539887 * pow(vu, 3) - 8.340493152 * pow(vu, 4));

		return (float)spinCD;
	}

	//Returns the magnitude of the drag force vector in Newtons. 
	float Simulation::CalculateDragForce(float dragcoefficient, float airdensitykgm3, float surfaream2, float speedms)
	{
		//Drag coefficient is estimated to be 0.47 with no spin, 0.43 with high spin.
		//Area for a 6mm BB is 0.000028274m^2.
		return dragcoefficient * 0.5 * airdensitykgm3 * surfaream2 * pow(speedms, 2);
	}

	float Simulation::CalculateLiftCoefficient(float spinrpm, float linearspeedmps, float bbdiametermeters)
	{
		//This equation is made to fit experimental data, from Dr. Dyrkacz's "The Physics of Paintball" which he based on 
		//Achenbach, E., J. Fluid Mech. 54,565 (1972). See his page here:
		//https://web.archive.org/web/20040617080904/http://home.comcast.net/~dyrgcmn/pball/pballIntro.html

		//Calculate the surface speed, which to my understanding is the speed at which a point on the outside of the circle is moving. 
		float surfacespeed = bbdiametermeters * PI * spinrpm / 60;

		//Calculate the ratio of spin velocity to linear velocity. Those two are represented as V and U mathematically, hence calling it vu.
		float vu = surfacespeed / linearspeedmps;

		//Aaaaand another crazy formula created via lots of experimentation by someone way better at the mathses than yours truly. 
		double liftcoef = (-0.0020907 - 0.208056226 * vu + 0.768791456 * pow(vu, 2) - 0.84865215 * pow(vu, 3) + 0.75365982 * pow(vu, 4)) /
			(1 - 4.82629033 * vu + 9.95459464 * pow(vu, 2) - 7.85649742 * pow(vu, 3) + 3.273765328 * pow(vu, 4));

		//cout << " VU: " << vu << " CL: " << liftcoef << endl;

		return (float)liftcoef;
	}

	//Returns the magnitute of the lift force from the Magnus effect (from hop-up) in Newtons.
	float Simulation::CalculateLiftForce(float liftcoefficient, float airdensitykgm3, float surfaream2, float speedms)
	{
		return liftcoefficient * airdensitykgm3 * pow(speedms, 2) * surfaream2;
	}


	float3 Simulation::CrossProduct(float3 v1, float3 v2)
	{
		float3 cross;
		cross.x = v1.y * v2.z - v1.z * v2.y;
		cross.y = v1.z * v2.x - v1.x * v2.z;
		cross.z = v1.x * v2.y - v1.y * v2.x;

		return cross;
	}

	/*
	float3 Simulation::RotateVectorAroundAxis(float3 srcvec, float3 axis, float angledegrees)
	{
		//Convert to radians.
		float anglerad = angledegrees * PI / 180;

		float sinang = sin(anglerad); //Shorthand.
		float cosang = cos(anglerad); //Shorthand.

		//Build the rotation matrix that does the job.
		//Referencing: https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
		cv::Mat rotmatrix = cv::Mat(3, 3, CV_32FC1);

		rotmatrix.at<float>(0, 0) = cosang + pow(axis.x, 2) * (1 - cosang);
		rotmatrix.at<float>(0, 1) = axis.x * axis.y * (1 - cosang) - axis.z * cosang;
		rotmatrix.at<float>(0, 2) = axis.x * axis.z * (1 - cosang) + axis.y * sinang;

		rotmatrix.at<float>(1, 0) = axis.y * axis.x * (1 - cosang) + axis.z * sinang;
		rotmatrix.at<float>(1, 1) = cosang + pow(axis.y, 2) * (1 - cosang);
		rotmatrix.at<float>(1, 2) = axis.y * axis.z * (1 - cosang) - axis.x * sinang;

		rotmatrix.at<float>(2, 0) = axis.z * axis.x * (1 - cosang) - axis.y * sinang;
		rotmatrix.at<float>(2, 1) = axis.z * axis.y * (1 - cosang) + axis.x * sinang;
		rotmatrix.at<float>(2, 2) = cosang + pow(axis.z, 2) * (1 - cosang);

		//Multiply.
		cv::Mat srcvecmat = cv::Mat(3, 1, CV_32FC1);
		srcvecmat.at<float>(0) = srcvec.x;
		srcvecmat.at<float>(1) = srcvec.y;
		srcvecmat.at<float>(2) = srcvec.z;

		cv::Mat multmat = rotmatrix * srcvecmat;

		return sl::float3(multmat.at<float>(0), multmat.at<float>(1), multmat.at<float>(2));
	}
	*/

	float3 Simulation::RotateVectorAroundAxis(float3 srcvec, float3 axis, float angledegrees) {
		// Convert to radians
		float anglerad = angledegrees * PI / 180.0f;

		float sinang = std::sin(anglerad); // Shorthand
		float cosang = std::cos(anglerad); // Shorthand

		// Build the rotation matrix
		// Referencing: https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
		BasicMat rotmatrix(3, 3);

		rotmatrix.at(0, 0) = cosang + std::pow(axis.x, 2) * (1 - cosang);
		rotmatrix.at(0, 1) = axis.x * axis.y * (1 - cosang) - axis.z * sinang;
		rotmatrix.at(0, 2) = axis.x * axis.z * (1 - cosang) + axis.y * sinang;

		rotmatrix.at(1, 0) = axis.y * axis.x * (1 - cosang) + axis.z * sinang;
		rotmatrix.at(1, 1) = cosang + std::pow(axis.y, 2) * (1 - cosang);
		rotmatrix.at(1, 2) = axis.y * axis.z * (1 - cosang) - axis.x * sinang;

		rotmatrix.at(2, 0) = axis.z * axis.x * (1 - cosang) - axis.y * sinang;
		rotmatrix.at(2, 1) = axis.z * axis.y * (1 - cosang) + axis.x * sinang;
		rotmatrix.at(2, 2) = cosang + std::pow(axis.z, 2) * (1 - cosang);

		// Convert srcvec to matrix form
		BasicMat srcvecmat(3, 1);
		srcvecmat.at(0, 0) = srcvec.x;
		srcvecmat.at(1, 0) = srcvec.y;
		srcvecmat.at(2, 0) = srcvec.z;

		// Multiply
		BasicMat multmat = rotmatrix * srcvecmat;

		return float3(multmat.at(0, 0), multmat.at(1, 0), multmat.at(2, 0));
	}

	float Simulation::CalculateMomentOfInertia(float bbmasskg, float bbdiamaterm)
	{
		//Moment of inertia for a sphere is known as 2/5mr². 
		//I found this on Wikipedia and a textbook I own. Wikipedia link:  
		//https://en.wikipedia.org/wiki/List_of_moments_of_inertia
		return .4 * bbmasskg * pow(bbdiamaterm * 0.5, 2);
	}

	float Simulation::CalculateSpinDragTorqueSimple(float airDensitykgm3, float bbradiusmeters, float spinRPM)
	{		const float dragCoefficient = 0.47f;  // Typical drag coefficient for a sphere

		//Convert RPM to radians per second.
		float angularVelocity = spinRPM * 2.0f * PI / 60.0f;

		//Calculate the cross-sectional area of the BB.
		float area = PI * bbradiusmeters * bbradiusmeters;

		//Calculate the relative speed due to rotation at the surface of the BB.
		float relativeSpeed = angularVelocity * bbradiusmeters;

		//Calculate the drag force based on rotation.
		float dragForce = 0.5f * dragCoefficient * airDensitykgm3 * area * (relativeSpeed * relativeSpeed);

		//Calculate and return the drag torque.
		float dragTorque = dragForce * bbradiusmeters;

		return dragTorque;
	}

	//Should return in Newton-meters. 
	float Simulation::CalculateSpinDragTorque(float speedMPS, float airdensitykgm3, float airviscosity, float bbdiameterm, float spinrpm,
		float viscousfriccoef, float* reynolds)
	{
		//This is based on, again, the Airsoft Trajectory Project, however this section is incomplete, namely
		//in that it doesn't provide info for providing the viscous friction coefficient. For this reason, 
		//I'm putting that coefficient as a configurable setting and will tweak it based on experimental data. 
		//I found only one site that provides friction coeffients for various materials against air, but it does
		//not include plastics: https://www.tribonet.org/wiki/friction-coefficients-in-atmosphere-and-vacuum/
		//However, I'm guessing it's relatively small, <0.1, as the listed coefficient for 
		//Teflon (polytetrafluoroethylene) is the closest material to BB plastic (non-biodegradable airsoft
		//BBs are made of acrylonitrile butadiene styrene, or ABS). And literally my only basis for comparison
		//is that both are polymers. So yeah... gonna tweak the setting based on testing. 

		//Also trying equations from: https://www.deepdyve.com/lp/wiley/a-high-reynolds-number-rotating-disk-rheometer-E3TzStyy9q

		float bbradius = bbdiameterm / 2.0; //Shorthand. 

		//TEST: Set density and viscosity to those in example 6.2 here: 
		//https://www.sciencedirect.com/topics/engineering/rotational-reynolds-number
		//airdensitykgm3 = 5.2;
		//airviscosity = 0.00003;


		//NEXT STEP: Take the math for the Cm (moment coefficient) there, which applies to cylinders, 
		//and figure out how to apply to spheres. See about finding about area relationships with radius or something. 
		//That number should give us the power, which we should be able to convert to torque rather simply. 
		//Also crap, we might have to do the same thing for the reynolds number. Yeah, Re = 233 for 3mm, 933 for 6mm. 

		//Calculate angular velocity from RPM into radians per sec. I *THINK* this is proper units. 
		//float angvelocityradpersec = spinrpm / 60.0 * 2 * PI;
		float angvelocityradpersec = spinrpm * 2 * PI / 60.0;

		//Calculate the Reynolds number for centerline rotation. 
		//2024 change: First line was commented, second was taken, even though the formula follows the first. Huh?
		//float reynoldsrot = airdensitykgm3 * pow(bbradius, 2) * angvelocityradpersec / viscousfriccoef;
		//float reynoldsrot = pow(bbradius, 2) * angvelocityradpersec * airdensitykgm3 / airviscosity; 
		float reynoldsrot = (airdensitykgm3 * angvelocityradpersec * pow(bbradius, 2)) / airviscosity;
		//float reynoldsrot = pow(bbradius, 2) * speedMPS * airdensitykgm3 * airviscosity;
		*reynolds = reynoldsrot;


		//float reynoldsrot = 6220.0;

		//Calculate the torque coefficient due to drag. 
		//This is taken straight from ATP, but I can't find that formula elsewhere online,  (especially the constants)
		//so this is an act of faith. I totally understand the author's hesitation to publish details, as there really 
		//are lots of blanks and mixed signals in calculating this stuff. 
		float torquecoef = 6.45 / pow(reynoldsrot, 0.5) + 32.1 / reynoldsrot;

		//cout << "AngRadPerSec: " << angvelocityradpersec << " Air density: " << airdensitykgm3 << " Viscosity: " << airviscosity << 
		//	" Reynolds: " << reynoldsrot << " torque coef: " << torquecoef << endl;

		//Calculate the actual torque due to drag. 
		float dragtorque = 0.5 * torquecoef * airdensitykgm3 * pow(bbradius, 3) * pow(angvelocityradpersec, 2);
		//float dragtorque = 0.5 * torquecoef * airviscosity * pow(bbradius, 3) * pow(angvelocityradpersec, 2);
		return dragtorque;
	}


}