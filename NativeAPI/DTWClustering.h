#pragma once
#include "Config.h"
#include "Data.h"
#include "Optional.h"
#include <tuple>
#include <unordered_set>
#include <unordered_map>
#include <concurrent_vector.h>
#include <concurrent_unordered_map.h>

class DTWCluster
{
public:
	static void MakeDTWLabelDistances(const std::vector<sizetype>& labels, std::map<sizetype, std::map<sizetype, float>>& dists);
	static void Clustering(SequenceDataSet* data, const std::vector<sizetype>& channels, const std::vector<sizetype>& labels, float lamda, const std::map<sizetype, std::map<sizetype, float>>& labeldistances, float labeldistancesmoothterm);

	class DTWCache : noncopyable
	{
	public:
		sizetype Length;
		Table<datatype>* Distances;
		Table<datatype>* PathCosts;
		Table<sizetype>* Step1;
		Table<sizetype>* Step2;

		DTWCache(sizetype len);
		~DTWCache();
		void Reset();
	};

	static datatype DTWCost(DTWCache* cache, const SequenceChannel* sch1, const SequenceChannel* sch2, const std::map<sizetype, std::map<sizetype, float>>& labeldistances);

	static datatype DTWPath(DTWCache* cache, const SequenceChannel* sch1, const SequenceChannel* sch2, const std::map<sizetype, std::map<sizetype, float>>& labeldistances, std::vector<std::tuple<sizetype, sizetype>>& paths);

	static datatype DTWPath2(DTWCache* cache, const SequenceChannel* sch1, const std::vector<datatype>& sch2, const std::vector<std::map<sizetype, float>>& schlab2, const std::map<sizetype, std::map<sizetype, float>>& labeldistances, float lamda, std::vector<std::tuple<sizetype, sizetype>>& paths);

	static void DBScan(DTWCache* cache, const concurrency::concurrent_vector<SequenceChannel*>& samples, sizetype order, const std::map<sizetype, std::map<sizetype, float>>& labeldistances, std::unordered_map<sizetype, std::vector<sizetype>>& clusters);

	static void DBA(DTWCache* cache, const concurrency::concurrent_vector<SequenceChannel*>& samples, const std::vector<sizetype>& indexs, const std::map<sizetype, std::map<sizetype, float>>& labeldistances, float lamda, size_t maxiterations, SequenceTemplate* dtw);

};