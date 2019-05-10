#pragma once
#ifndef NATIVE_CAPI
#ifdef NATIVE_EXPORTS
#define NATIVE_CAPI extern "C" __declspec(dllexport)
#else
#define NATIVE_CAPI extern "C" __declspec(dllimport)
#endif
#endif

#ifndef NATIVE_API
#ifdef NATIVE_EXPORTS
#define NATIVE_API __declspec(dllexport)
#else
#define NATIVE_API __declspec(dllimport)
#endif
#endif

#ifndef NATIVE_TYPES
#define NATIVE_TYPES
typedef float datatype;
typedef __int64 sizetype;
typedef unsigned __int64 usizetype;
#endif

#include <memory>
#include <vector>

namespace GPC
{
	class Cycle;
	class Dataset;
	class API;

	class NATIVE_API Pairs
	{
	public:
		void* warpper;

		Pairs();
		~Pairs();
	};

	typedef std::shared_ptr<Pairs> PairsPtr;

	class NATIVE_API Cycle
	{
	public:
		void* warpper;

		Cycle();
		~Cycle();

		bool AddFrame(float time, int phase, const std::vector<float>& vals);
		bool AddFrame1(float time, int phase, const float* vals, unsigned int valcount);
	};

	typedef std::shared_ptr<Cycle> CyclePtr;


	class NATIVE_API Dataset
	{
	public:
		void* warpper;

		Dataset();
		~Dataset();

		CyclePtr AllocateCycle(bool complete);
		void AddCycle(CyclePtr ptr);

		Cycle* AllocateCycle1(bool complete);
		void AddCycle1(Cycle* ptr);
	};

	typedef std::shared_ptr<Dataset> DatasetPtr;

	class NATIVE_API API
	{
	public:
		static DatasetPtr CreateDataset(unsigned int channelcount, unsigned int phasecount);

		static bool Clustering(DatasetPtr ptr, float edgeratio=0.1);

		static PairsPtr ExtractFeaturePairs(DatasetPtr ptr, int desiredfeaturecount, int duplicate_neigborsize=6);

		static bool ComputeFeatures(DatasetPtr ptr, PairsPtr ptr1, const wchar_t* outpath, unsigned int standardizecycleframecount=100);
	
		static bool SaveDataset(DatasetPtr ptr, const wchar_t* path);
	
		static DatasetPtr LoadDataset(const wchar_t* path, unsigned int channelcount, unsigned int phasecount);
	};
}