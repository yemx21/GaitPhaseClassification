#include "DTWClustering.h"
#include <map>
#include <ppl.h>
#include <numeric>
#include <string>
using namespace concurrency;

void DTWCluster::MakeDTWLabelDistances(const std::vector<sizetype>& labels, std::map<sizetype, std::map<sizetype, float>>& dists)
{
	int labelcount = (int)labels.size();

	double K = (double)labels[0];
	for (sizetype l : labels)
	{
		if (l > K)
			K = l;
	}

	float LK = ceil((K - 1.0f) / 2.0f);
	float UK = ceil((K + 1.0f) / 2.0f);


	for (int i = 0; i < labelcount; i++)
	{
		dists.insert(std::make_pair(labels[i], std::map<sizetype, float>()));
		for (int j = 0; j < labelcount; j++)
		{
			
			float d = LK - abs(((((int)floor(2.0f*abs(labels[i] - labels[j]) / (K + 1.0f))) % 2 == 0) ? 1.0 : 0.0)*LK -
				(((int)(abs(labels[i] - labels[j]))) % ((int)UK)));

			dists[labels[i]].insert(std::make_pair(labels[j], d));
		}
	}


	//if (labelcount == 9)
	//{
	//	dists[1][2] = 1.0f;
	//	dists[1][3] = 2.0f;
	//	dists[1][4] = 3.0f;
	//	dists[1][5] = 4.0f;
	//	dists[1][6] = 4.0f;
	//	dists[1][7] = 3.0f;
	//	dists[1][8] = 2.0f;
	//	dists[1][9] = 1.0f;

	//	dists[2][1] = 1.0f;
	//	dists[2][3] = 1.0f;
	//	dists[2][4] = 2.0f;
	//	dists[2][5] = 3.0f;
	//	dists[2][6] = 4.0f;
	//	dists[2][7] = 4.0f;
	//	dists[2][8] = 3.0f;
	//	dists[2][9] = 2.0f;

	//	dists[3][1] = 2.0f;
	//	dists[3][2] = 1.0f;
	//	dists[3][4] = 1.0f;
	//	dists[3][5] = 2.0f;
	//	dists[3][6] = 3.0f;
	//	dists[3][7] = 4.0f;
	//	dists[3][8] = 4.0f;
	//	dists[3][9] = 3.0f;

	//	dists[4][1] = 3.0f;
	//	dists[4][2] = 2.0f;
	//	dists[4][3] = 1.0f;
	//	dists[4][5] = 1.0f;
	//	dists[4][6] = 2.0f;
	//	dists[4][7] = 3.0f;
	//	dists[4][8] = 4.0f;
	//	dists[4][9] = 4.0f;

	//	dists[5][1] = 4.0f;
	//	dists[5][2] = 3.0f;
	//	dists[5][3] = 2.0f;
	//	dists[5][4] = 1.0f;
	//	dists[5][6] = 1.0f;
	//	dists[5][7] = 2.0f;
	//	dists[5][8] = 3.0f;
	//	dists[5][9] = 4.0f;

	//	dists[6][1] = 4.0f;
	//	dists[6][2] = 4.0f;
	//	dists[6][3] = 3.0f;
	//	dists[6][4] = 2.0f;
	//	dists[6][5] = 1.0f;
	//	dists[6][7] = 1.0f;
	//	dists[6][8] = 2.0f;
	//	dists[6][9] = 3.0f;

	//	dists[7][1] = 3.0f;
	//	dists[7][2] = 4.0f;
	//	dists[7][3] = 4.0f;
	//	dists[7][4] = 3.0f;
	//	dists[7][5] = 2.0f;
	//	dists[7][6] = 1.0f;
	//	dists[7][8] = 1.0f;
	//	dists[7][9] = 2.0f;

	//	dists[8][1] = 2.0f;
	//	dists[8][2] = 3.0f;
	//	dists[8][3] = 4.0f;
	//	dists[8][4] = 4.0f;
	//	dists[8][5] = 3.0f;
	//	dists[8][6] = 2.0f;
	//	dists[8][7] = 1.0f;
	//	dists[8][9] = 1.0f;

	//	dists[9][1] = 1.0f;
	//	dists[9][2] = 2.0f;
	//	dists[9][3] = 3.0f;
	//	dists[9][4] = 4.0f;
	//	dists[9][5] = 4.0f;
	//	dists[9][6] = 3.0f;
	//	dists[9][7] = 2.0f;
	//	dists[9][8] = 1.0f;
	//}

	
}


constexpr float SQRT2 = 1.414213562373095048801688724209698078;

DTWCluster::DTWCache::~DTWCache()
{
	if (Distances) { delete Distances; Distances = nullptr; }
	if (PathCosts) { delete PathCosts; PathCosts = nullptr; }
	if (Step1) { delete Step1; Step1 = nullptr; }
	if (Step2) { delete Step2; Step2 = nullptr; }
}

DTWCluster::DTWCache::DTWCache(sizetype len)
{
	Length = len;
	sizetype upper = Length + 1;

	Distances = new Table<datatype>(upper);
	PathCosts = new Table<datatype>(upper);
	Step1 = new Table<sizetype>(upper);
	Step2 = new Table<sizetype>(upper);
}

void DTWCluster::DTWCache::Reset()
{
	datatype* address_PathCosts = PathCosts->At(Length);
	for (size_t i = 0; i <= Length; i++)
	{
		address_PathCosts[i] = std::numeric_limits<datatype>::infinity();
	}

	for (size_t i = 0; i <= Length; i++)
	{
		PathCosts->At(i, Length) = std::numeric_limits<datatype>::infinity();
	}

	Distances->Clear();
	Step1->Clear();
	Step2->Clear();
}

datatype DTWCluster::DTWCost(DTWCache* cache, const SequenceChannel* sch1, const SequenceChannel* sch2, const std::map<sizetype, std::map<sizetype, float>>& labeldistances)
{
	sizetype len = cache->Length;
	cache->Reset();

	for (size_t i = 0; i < len; i++)
	{
		datatype* currentDistances = cache->Distances->At(i);
		auto xVal = sch1->Values[i];
		auto xLabel = sch1->Labels[i];
		for (size_t j = 0; j < len; j++)
		{
			datatype dist = (xVal - sch2->Values[j]) * labeldistances.at(xLabel).at(sch2->Labels[j]);
			currentDistances[j] += dist * dist;
		}
	}

	for (int i = 0; i < len; i++)
	{
		datatype* currentDistances = cache->Distances->At(i);
		for (int j = 0; j < len; j++)
			currentDistances[j] = sqrt(currentDistances[j]);
	}

	for (sizetype i = len - 1; i >= 0; i--)
	{
		datatype* currentRowDistances = cache->Distances->At(i);
		datatype* currentRowPathCost = cache->PathCosts->At(i);
		datatype* previousRowPathCost = cache->PathCosts->At(i + 1);

		sizetype* currentRowPredecessorStepX = cache->Step1->At(i);
		sizetype* currentRowPredecessorStepY = cache->Step2->At(i);

		for (int j = len - 1; j >= 0; j--)
		{
			datatype diagonalNeighbourCost = previousRowPathCost[j + 1];
			datatype xNeighbourCost = previousRowPathCost[j];
			datatype yNeighbourCost = currentRowPathCost[j + 1];

			if (std::isinf(diagonalNeighbourCost) && (i == j))
				currentRowPathCost[j] = currentRowDistances[j];
			else
			{
				if (diagonalNeighbourCost <= xNeighbourCost && diagonalNeighbourCost <= yNeighbourCost)
				{
					currentRowPathCost[j] = diagonalNeighbourCost + currentRowDistances[j];
					currentRowPredecessorStepX[j] = 1;
					currentRowPredecessorStepY[j] = 1;
				}
				else if (xNeighbourCost <= yNeighbourCost)
				{
					currentRowPathCost[j] = xNeighbourCost + currentRowDistances[j];
					currentRowPredecessorStepX[j] = 1;
					currentRowPredecessorStepY[j] = 0;
				}
				else
				{
					currentRowPathCost[j] = yNeighbourCost + currentRowDistances[j];
					currentRowPredecessorStepX[j] = 0;
					currentRowPredecessorStepY[j] = 1;
				}
			}
		}
	}

	datatype cost = *cache->PathCosts->AddressOf();
	if (std::isnan(cost)) return 0.0;
	return cost / len / SQRT2;
}

datatype DTWCluster::DTWPath(DTWCache* cache, const SequenceChannel* sch1, const SequenceChannel* sch2, const std::map<sizetype, std::map<sizetype, float>>& labeldistances, std::vector<std::tuple<sizetype, sizetype>>& paths)
{
	paths.clear();
	sizetype len = cache->Length;
	cache->Reset();

	for (size_t i = 0; i < len; i++)
	{
		datatype* currentDistances = cache->Distances->At(i);
		auto xVal = sch1->Values[i];
		auto xLabel = sch1->Labels[i];
		for (size_t j = 0; j < len; j++)
		{
			datatype dist = (xVal - sch2->Values[j]) * labeldistances.at(xLabel).at(sch2->Labels[j]);
			currentDistances[j] += dist * dist;
		}
	}

	for (int i = 0; i < len; i++)
	{
		datatype* currentDistances = cache->Distances->At(i);
		for (int j = 0; j < len; j++)
			currentDistances[j] = sqrtf(currentDistances[j]);
	}

	for (sizetype i = len - 1; i >= 0; i--)
	{
		datatype* currentRowDistances = cache->Distances->At(i);
		datatype* currentRowPathCost = cache->PathCosts->At(i);
		datatype* previousRowPathCost = cache->PathCosts->At(i + 1);

		sizetype* currentRowPredecessorStepX = cache->Step1->At(i);
		sizetype* currentRowPredecessorStepY = cache->Step2->At(i);

		for (int j = len - 1; j >= 0; j--)
		{
			datatype diagonalNeighbourCost = previousRowPathCost[j + 1];
			datatype xNeighbourCost = previousRowPathCost[j];
			datatype yNeighbourCost = currentRowPathCost[j + 1];

			if (std::isinf(diagonalNeighbourCost) && (i == j))
				currentRowPathCost[j] = currentRowDistances[j];
			else
			{
				if (diagonalNeighbourCost <= xNeighbourCost && diagonalNeighbourCost <= yNeighbourCost)
				{
					currentRowPathCost[j] = diagonalNeighbourCost + currentRowDistances[j];
					currentRowPredecessorStepX[j] = 1;
					currentRowPredecessorStepY[j] = 1;
				}
				else if (xNeighbourCost <= yNeighbourCost)
				{
					currentRowPathCost[j] = xNeighbourCost + currentRowDistances[j];
					currentRowPredecessorStepX[j] = 1;
					currentRowPredecessorStepY[j] = 0;
				}
				else
				{
					currentRowPathCost[j] = yNeighbourCost + currentRowDistances[j];
					currentRowPredecessorStepX[j] = 0;
					currentRowPredecessorStepY[j] = 1;
				}
			}
		}
	}

	sizetype indexX = 0;
	sizetype indexY = 0;

	paths.push_back(std::make_tuple(indexX, indexY));
	while (indexX < len - 1 || indexY < len - 1)
	{
		sizetype stepX = cache->Step1->At(indexX, indexY);
		sizetype stepY = cache->Step2->At(indexX, indexY);
		indexX += stepX;
		indexY += stepY;
		paths.push_back(std::make_tuple(indexX, indexY));
	}

	datatype cost = *cache->PathCosts->AddressOf();
	if (std::isnan(cost)) return 0.0;
	return cost / len / SQRT2;
}

datatype DTWCluster::DTWPath2(DTWCache* cache, const SequenceChannel* sch1, const std::vector<datatype>& sch2, const std::vector<std::map<sizetype, float>>& schlab2, const std::map<sizetype, std::map<sizetype, float>>& labeldistances, float lamda, std::vector<std::tuple<sizetype, sizetype>>& paths)
{
	paths.clear();
	sizetype len = cache->Length;
	cache->Reset();

	if (lamda >= 0.0)
	{
		for (size_t i = 0; i < len; i++)
		{
			datatype* currentDistances = cache->Distances->At(i);
			auto xVal = sch1->Values[i];
			auto xLab = sch1->Labels[i];
			for (size_t j = 0; j < len; j++)
			{
				if (schlab2[j].empty())
				{
					datatype dist = exp(-(xVal - sch2[j]) / lamda);
					datatype edist = dist * dist;
					currentDistances[j] += edist;
				}
				else
				{
					datatype labfactor = 0.0f;
					for (auto lf : schlab2[j])
					{
						labfactor += labeldistances.at(xLab).at(lf.first) * lf.second;
					}

					datatype dist = exp(-(xVal - sch2[j]) * labfactor / lamda);
					datatype edist = dist * dist;
					currentDistances[j] += edist;
				}
			}
		}
	}
	else
	{
		for (size_t i = 0; i < len; i++)
		{
			datatype* currentDistances = cache->Distances->At(i);
			auto xVal = sch1->Values[i];
			for (size_t j = 0; j < len; j++)
			{
				datatype dist = xVal - sch2[j];
				datatype edist = dist * dist;
				currentDistances[j] += edist;
			}
		}
	}

	for (int i = 0; i < len; i++)
	{
		datatype* currentDistances = cache->Distances->At(i);
		for (int j = 0; j < len; j++)
			currentDistances[j] = sqrtf(currentDistances[j]);
	}

	for (sizetype i = len - 1; i >= 0; i--)
	{
		datatype* currentRowDistances = cache->Distances->At(i);
		datatype* currentRowPathCost = cache->PathCosts->At(i);
		datatype* previousRowPathCost = cache->PathCosts->At(i + 1);

		sizetype* currentRowPredecessorStepX = cache->Step1->At(i);
		sizetype* currentRowPredecessorStepY = cache->Step2->At(i);

		for (int j = len - 1; j >= 0; j--)
		{
			datatype diagonalNeighbourCost = previousRowPathCost[j + 1];
			datatype xNeighbourCost = previousRowPathCost[j];
			datatype yNeighbourCost = currentRowPathCost[j + 1];

			if (std::isinf(diagonalNeighbourCost) && (i == j))
				currentRowPathCost[j] = currentRowDistances[j];
			else
			{
				if (diagonalNeighbourCost <= xNeighbourCost && diagonalNeighbourCost <= yNeighbourCost)
				{
					currentRowPathCost[j] = diagonalNeighbourCost + currentRowDistances[j];
					currentRowPredecessorStepX[j] = 1;
					currentRowPredecessorStepY[j] = 1;
				}
				else if (xNeighbourCost <= yNeighbourCost)
				{
					currentRowPathCost[j] = xNeighbourCost + currentRowDistances[j];
					currentRowPredecessorStepX[j] = 1;
					currentRowPredecessorStepY[j] = 0;
				}
				else
				{
					currentRowPathCost[j] = yNeighbourCost + currentRowDistances[j];
					currentRowPredecessorStepX[j] = 0;
					currentRowPredecessorStepY[j] = 1;
				}
			}
		}
	}

	sizetype indexX = 0;
	sizetype indexY = 0;

	paths.push_back(std::make_tuple(indexX, indexY));
	while (indexX < len - 1 || indexY < len - 1)
	{
		sizetype stepX = cache->Step1->At(indexX, indexY);
		sizetype stepY = cache->Step2->At(indexX, indexY);
		indexX += stepX;
		indexY += stepY;
		paths.push_back(std::make_tuple(indexX, indexY));
	}

	datatype cost = *cache->PathCosts->AddressOf();
	if (std::isnan(cost)) return 0.0;
	return cost / len / SQRT2;
}


#pragma optimize("", off)
__declspec(noinline)
void Spin(unsigned int work)
{
	for (unsigned long long i = 0; i < work*work*work; ++i);
}
#pragma optimize("",on)

class SortFunc
{
public:
	explicit SortFunc(unsigned int work = 0)
		:_work(work)
	{
	}

	bool operator()(datatype left, datatype right) const
	{
		Spin(_work);
		return left < right;
	}

private:
	const unsigned int _work;
};

bool RegionNearestAverageDistance(size_t count, sizetype index, const std::vector<std::vector<datatype>>& distmap, sizetype num, datatype& averdist)
{
	averdist = std::numeric_limits<datatype>::quiet_NaN();
	concurrency::concurrent_vector<datatype> sampledists(count - 1);
	const std::vector<datatype>& distmap_selected = distmap[index];
	for (size_t i = 0, n = 0; i < count; i++)
	{
		if (i == index) continue;
		sampledists[n++] = distmap_selected[i];
	}

	concurrency::parallel_sort(sampledists.begin(), sampledists.end(), SortFunc(10));

	sizetype sumcount = 0;
	double sum = 0.0;
	for (size_t i = 0; i < count - 1; i++)
	{
		sum += sampledists[i];
		sumcount++;
		if (sumcount >= num) break;
	}
	if (count == 0) return false;
	averdist = sum / sumcount;
	return true;
}

bool ComputeEpsMinPtS(size_t count, const std::vector<std::vector<datatype>>& distmap, int order, double& epsilon, sizetype& minPts)
{
	std::vector<datatype> averdists;
	averdists.reserve(order);
	for (int o = 2; o <= order; o++)
	{
		std::vector<datatype> taverdists;
		taverdists.reserve(count);
		for (size_t i = 0; i < count; i++)
		{
			datatype averdist = std::numeric_limits<datatype>::quiet_NaN();
			if (RegionNearestAverageDistance(count, i, distmap, o, averdist))
				taverdists.push_back(averdist);
		}
		averdists.push_back(taverdists.size() ? std::accumulate(taverdists.begin(), taverdists.end(), 0.0f) / taverdists.size() : 0.0);
	}

	std::vector<sizetype> curve_order;
	std::vector<datatype> curve_eps;
	size_t nonzeroeps = 0;
	for (int o = 2, n = 0; o <= order; o++, n++)
	{
		datatype eps = averdists[n];
		if (std::isnan(eps)) continue;
		if (eps > 0.0) nonzeroeps++;
		curve_order.push_back(o);
		curve_eps.push_back(eps);
	}

	if (nonzeroeps > 1)
	{
		try
		{
			epsilon = std::numeric_limits<datatype>::max();
			minPts = 2;
			for (size_t i = 0; i<curve_eps.size(); i++)
			{
				if (curve_eps[i] < epsilon)
				{
					epsilon = curve_eps[i];
					minPts = curve_order[i];
				}
			}
			return true;
		}
		catch (...)
		{

		}
	}
	epsilon = 0.1;
	minPts = 2;
	return false;
}

void ExpandCluster(size_t count, sizetype index, const std::vector<std::vector<datatype>>& distmap, std::vector<bool>& isvisited, std::vector<sizetype>& clusterids, std::vector<sizetype>& neighborPts, int clusterId, double epsilon, int minPts)
{
	clusterids[index] = clusterId;
	for (sizetype i = 0; i < neighborPts.size(); i++)
	{
		sizetype in = neighborPts[i];
		if (!isvisited[in])
		{
			isvisited[in] = true;
			std::unordered_set<sizetype> neighborPts2;
			for (size_t n = 0; n < count; n++)
			{
				if (distmap[in][n] <= epsilon)
					neighborPts2.insert(n);
			}
			if (neighborPts2.size() >= minPts)
			{
				for (sizetype n : neighborPts2)
				{
					if (std::find(neighborPts.begin(), neighborPts.end(), n) != neighborPts.end()) continue;
					neighborPts.push_back(n);
				}
			}
		}
		if (clusterids[in] == 0)
			clusterids[in] = clusterId;
	}
}

void DTWCluster::DBScan(DTWCache* cache, const concurrency::concurrent_vector<SequenceChannel*>& samples, sizetype order, const std::map<sizetype, std::map<sizetype, float>>& labeldistances, std::unordered_map<sizetype, std::vector<sizetype>>& clusters)
{
	size_t len = samples.size();
	std::vector<std::vector<datatype>> distmap;
	distmap.resize(len);
	for (size_t i = 0; i < len; i++) distmap[i].resize(len, -1.0f);

	for (size_t i = 0; i < len; i++)
	{
		SequenceChannel* sample1 = samples[i];
		for (size_t j = 0; j < len; j++)
		{
			if (i == j)
			{
				distmap[i][j] = 0.0f;
				distmap[j][i] = 0.0f;
				continue;
			}
			datatype cost = DTWCost(cache, sample1, samples[j], labeldistances);
			if (distmap[i][j] < 0.0f) distmap[i][j] = cost;
			if (distmap[j][i] < 0.0f) distmap[j][i] = cost;
		}
	}

	double eps = 0.1;
	sizetype minpts = 2;
	if (!ComputeEpsMinPtS(len, distmap, 5, eps, minpts))
	{
		std::vector<sizetype> sampleindexs0(len);;
		for (size_t i = 0; i< len; i++) sampleindexs0[i] = i;
		clusters.insert(std::make_pair(2, sampleindexs0));
		return;
	}

	std::vector<bool> isvisited;
	std::vector<sizetype> clusterids;
	isvisited.resize(len, false);
	clusterids.resize(len, 0);

	int clusterId = 0;
	for (size_t i = 0; i < len; i++)
	{
		if (isvisited[i]) continue;
		isvisited[i] = true;

		std::vector<sizetype> neighborPts;
		for (size_t n = 0; n < len; n++)
		{
			if (distmap[i][n] <= eps) neighborPts.push_back(n);
		}

		if (neighborPts.size() < minpts)
			clusterids[i] = -1;
		else
		{
			clusterId++;
			ExpandCluster(len, i, distmap, isvisited, clusterids, neighborPts, clusterId, eps, minpts);
		}
	}

	for (sizetype i = 0; i < len; i++)
	{
		sizetype cid = clusterids[i];
		if (cid > 0)
		{
			if (clusters.find(cid) == clusters.end()) clusters.insert(std::make_pair(cid, std::vector<sizetype>{}));
			clusters[cid].push_back(i);
		}
	}

	for (auto cliter = clusters.begin(); cliter != clusters.end();)
	{
		if (cliter->second.size() <= minpts)
			cliter = clusters.erase(cliter);
		else
			cliter++;
	}

	if (clusters.size() == 1)
	{
		if (clusters.begin()->second.size() <= minpts)
		{
			std::vector<sizetype> sampleindexs0(len);
			for (size_t i = 0; i< len; i++) sampleindexs0[i] = i;
			clusters.insert(std::make_pair(2, sampleindexs0));
		}
	}
	else if (clusters.size() == 0)
	{
		std::vector<sizetype> sampleindexs0(len);
		for (size_t i = 0; i< len; i++) sampleindexs0[i] = i;
		clusters.insert(std::make_pair(2, sampleindexs0));
	}
}

void Detrend(SequenceChannel* input, std::vector<datatype>& output)
{
	int len = input->Length;
	double step = (input->Values[len - 1] - input->Values[0]) / len;
	output.resize(len);
	for (int i = 0; i < len; i++)
	{
		output[i] = input->Values[i] - step * i;
	}
}

size_t IndexOfMax(const std::vector<datatype>& input)
{
	size_t maxIndex = -1;
	datatype maxValue = input[0];

	for (size_t i = 0; i < input.size(); i++)
	{
		if (input[i] >= maxValue)
		{
			maxIndex = i;
			maxValue = input[i];
		}
	}

	return maxIndex;
}

size_t IndexOfMin(const std::vector<datatype>& input)
{
	size_t minIndex = -1;
	datatype minValue = input[0];

	for (size_t i = 0; i < input.size(); i++)
	{
		if (input[i] <= minValue)
		{
			minIndex = i;
			minValue = input[i];
		}
	}

	return minIndex;
}

datatype Median(const std::vector<size_t>& input)
{
	int midIndex = input.size() / 2;

	std::vector<size_t> input_copy = input;
	std::sort(input_copy.begin(), input_copy.end());

	return input.size() % 2 == 0
		? (input_copy[midIndex] + input_copy[midIndex - 1]) / 2
		: input_copy[midIndex];
}

datatype Median(const std::vector<datatype>& input)
{
	int midIndex = input.size() / 2;

	std::vector<datatype> input_copy = input;
	std::sort(input_copy.begin(), input_copy.end());

	return input.size() % 2 == 0
		? (input_copy[midIndex] + input_copy[midIndex - 1]) / 2
		: input_copy[midIndex];
}

void DTWCluster::DBA(DTWCache* cache, const concurrency::concurrent_vector<SequenceChannel*>& samples, const std::vector<sizetype>& indexs, const std::map<sizetype, std::map<sizetype, float>>& labeldistances, float lamda, size_t maxiterations, SequenceTemplate* dtw)
{
	if (indexs.size() == 1)
	{
		dtw->Values = samples[indexs[0]]->Values;
		dtw->Probabilities.resize(samples[indexs[0]]->Values.size());
		const auto& labelrefs = samples[indexs[0]]->Labels;
		for (sizetype i = 0; i<dtw->Values.size(); i++)
		{
			dtw->Probabilities[i].insert(std::make_pair(labelrefs[i], 1.0f));
		}
		dtw->GlobalWeight = 0.0;
		dtw->LocalWeight = 0.0;
	}

	size_t length = samples[indexs[0]]->Values.size();
	std::vector<size_t> minindexes;
	std::vector<size_t> maxindexes;
	std::vector<datatype> tempsample;

	for (sizetype idx : indexs)
	{
		Detrend(samples[idx], tempsample);
		minindexes.push_back(IndexOfMin(tempsample));
		maxindexes.push_back(IndexOfMax(tempsample));
	}

	datatype medianMaxIndex = Median(maxindexes);
	datatype medianMinIndex = Median(minindexes);

	std::vector<datatype> distances;
	for (size_t i = 0; i < minindexes.size(); i++)
	{
		distances.push_back(pow(maxindexes[i] - medianMaxIndex, 2.0) + pow(minindexes[i] - medianMinIndex, 2.0));
	}
	size_t selectedindex = IndexOfMin(distances);

	std::vector<datatype> average = samples[indexs[selectedindex]]->Values;
	std::vector<std::map<sizetype, float>> averagelabels(length);
	const std::vector<sizetype>& averagelabels0 = samples[indexs[selectedindex]]->Labels;
	for (size_t i = 0; i < length; i++) averagelabels[i].insert(std::make_pair(averagelabels0[i], 1.0f));

	std::vector<std::vector<datatype>> points(length);
	std::vector<std::vector<sizetype>> paths(length);

	std::vector<std::tuple<sizetype, sizetype, datatype>> ref_paths;

	double prevTotalDist = -1;
	double totalDist = -2;
	int count = 0;

	while (totalDist != prevTotalDist && count < maxiterations)
	{
		prevTotalDist = totalDist;

		for (size_t i = 0; i < length; i++)
		{
			points[i].clear();
			paths[i].clear();
		}

		std::vector<std::tuple<sizetype, sizetype>> tmp_paths;
		for (sizetype idx : indexs)
		{
			SequenceChannel* sc = samples[idx];
			const std::vector<datatype>& svals = sc->Values;
			const std::vector<sizetype>& slabs = sc->Labels;
			datatype ncost = DTWPath2(cache, sc, average, averagelabels, labeldistances, lamda, tmp_paths);

			for (auto it : tmp_paths)
			{
				points[std::get<1>(it)].push_back(svals[std::get<0>(it)]);
				paths[std::get<1>(it)].push_back(slabs[std::get<0>(it)]);
				ref_paths.push_back(std::make_tuple(std::get<1>(it), sc->Labels[std::get<0>(it)], ncost));
			}
		}

		for (size_t i = 0; i < length; i++)
		{
			const std::vector<datatype>& ps = points[i];
			average[i] = std::accumulate(ps.begin(), ps.end(), 0.0f) / ps.size();

			averagelabels[i].clear();
			const std::vector<sizetype>& ls = paths[i];
			float lscount = (float)ls.size();
			for (sizetype ils : ls)
			{
				if (averagelabels[i].find(ils) == averagelabels[i].end())
					averagelabels[i].insert(std::make_pair(ils, 1.0f));
				else
					averagelabels[i][ils] += 1.0f;
			}

			for (auto& xls : averagelabels[i])
			{
				xls.second /= lscount;
			}
		}

		totalDist = 0.0;

		for (sizetype idx : indexs)
		{
			const std::vector<datatype>& svals = samples[idx]->Values;
			const std::vector<sizetype>& slabs = samples[idx]->Labels;

			for (size_t i = 0; i < length; i++)
			{
				auto xval = svals[i];
				auto xlab = slabs[i];
				auto yval = average[i];
				const std::map<sizetype, float> ylabs = averagelabels[i];
				if (averagelabels[i].empty())
				{
					datatype dist = (xval - yval);
					datatype edist = dist * dist;
					totalDist += edist;
				}
				else
				{
					datatype labfactor = 0.0f;
					for (const auto& yl : ylabs)
					{
						labfactor += labeldistances.at(xlab).at(yl.first) * yl.second;
					}

					datatype dist = (xval - yval) * labfactor;
					datatype edist = dist * dist;
					totalDist += edist;
				}


				//totalDist += pow(svals[i] - average[i], 2.0);
			}
		}
		count++;
	}

	dtw->Values = average;
	size_t ref_pathcount = ref_paths.size();
	datatype ref_path_mincost = std::numeric_limits<datatype>::max();
	datatype ref_path_maxcost = std::numeric_limits<datatype>::min();
	for (size_t i = 0; i<ref_pathcount; i++)
	{
		datatype pcost = std::get<2>(ref_paths[i]);
		if (pcost < ref_path_mincost) ref_path_mincost = pcost;
		if (pcost > ref_path_maxcost) ref_path_maxcost = pcost;
	}
	datatype ref_path_dcost = ref_path_maxcost - ref_path_mincost;

	dtw->GlobalWeight = 0.0f;
	dtw->LocalWeight = 0.0f;
	dtw->Probabilities.resize(length);
	dtw->Labels.resize(length, -1);
	for (int nn = 0; nn < length; nn++)
	{
		double total_ncost = 0.0;
		std::unordered_map<sizetype, std::vector<double>> votes;
		for (auto np : ref_paths)
		{
			sizetype n2 = std::get<1>(np);
			if (std::get<0>(np) == nn && n2 != -1)
			{
				datatype npcost = std::get<2>(np);
				dtw->LocalWeight += npcost;
				if (votes.find(n2) == votes.end())
				{
					std::vector<double> vote_p;
					datatype ncost = (ref_path_maxcost - npcost) / ref_path_dcost;
					if (ncost < 0.0) ncost = 0.0;
					else if (ncost > 1.0) ncost = 1.0;
					vote_p.push_back(ncost);
					votes.insert(std::make_pair(n2, vote_p));
					total_ncost += ncost;
				}
				else
				{
					datatype ncost = (ref_path_maxcost - npcost) / ref_path_dcost;
					if (ncost < 0.0) ncost = 0.0;
					else if (ncost > 1.0) ncost = 1.0;
					votes[n2].push_back(ncost);
					total_ncost += ncost;
				}
			}
		}

		dtw->Probabilities[nn].clear();
		for (auto vv : votes)
		{
			if (total_ncost == 0.0)
				dtw->Probabilities[nn].insert(std::make_pair(vv.first, 0.0));
			else
			{
				float sum = std::accumulate(vv.second.begin(), vv.second.end(), 0.0);
				dtw->Probabilities[nn].insert(std::make_pair(vv.first, sum));
			}
		}

		float maxprob = 0.0f;
		dtw->Labels[nn] = -1;

		for (const auto& pp : dtw->Probabilities[nn])
		{
			if (pp.second > maxprob)
			{
				dtw->Labels[nn] = pp.first;
				maxprob = pp.second;
			}
		}
		dtw->GlobalWeight += total_ncost;
	}
	dtw->GlobalWeight /= length;
	dtw->LocalWeight /= ref_paths.size();
}

void DTWCluster::Clustering(SequenceDataSet* data, const std::vector<sizetype>& channels, const std::vector<sizetype>& labels, float lamda, const std::map<sizetype, std::map<sizetype, float>>& labeldistances, float labeldistanceexpterm)
{
	const std::vector<std::map<sizetype, SequenceChannel*>>& samples = data->SampleRefs;
	sizetype samplelength = data->Length;

	concurrent_unordered_map<sizetype, concurrent_vector<SequenceChannel*>> concurrent_samples;
	for (sizetype ch : channels) concurrent_samples.insert(std::make_pair(ch, concurrent_vector<SequenceChannel*>{}));

	std::unordered_map<sizetype, sizetype> majorhits;
	for (sizetype i = 0; i < samplelength; i++)
	{
		auto& sc = samples[i];
		for (auto cc : sc)
		{
			if (cc.second->Completed)
			{
				auto l = cc.second->Length;
				if (majorhits.find(l) != majorhits.end())
				{
					majorhits[l]++;
				}
				else
				{
					majorhits.insert(std::make_pair(l, 0));
				}
			}
		}
	}
	sizetype majorhitcount = majorhits.begin()->second;
	sizetype majorcompletedcount = majorhits.begin()->first;
	for (auto mp : majorhits)
	{
		if (mp.second > majorhitcount)
		{
			majorhitcount = mp.second;
			majorcompletedcount = mp.first;
		}
	}

	for (sizetype i = 0; i < samplelength; i++)
	{
		auto& sc = samples[i];
		for (auto cc : sc)
		{
			if (cc.second->Completed && cc.second->ValidLength == majorcompletedcount)
			{
				concurrent_samples[cc.first].push_back(cc.second);
			}
		}
	}

	std::map<sizetype, std::map<sizetype, float>> labeldistances1 = labeldistances;

	for (sizetype i :labels)
	{
		for (sizetype j : labels)
		{
			labeldistances1.at(i).at(j) = pow(labeldistanceexpterm, labeldistances1.at(i).at(j));
		}
	}

	concurrent_unordered_map<sizetype, concurrent_vector<SequenceTemplate>> concurrent_templates;
	for (sizetype ch : channels) concurrent_templates.insert(std::make_pair(ch, concurrent_vector<SequenceTemplate>{}));

	parallel_for_each(std::begin(channels), std::end(channels), [&concurrent_samples, &majorcompletedcount, &labeldistances1, &lamda, &concurrent_templates](sizetype ch)
	{
		concurrent_vector<SequenceChannel*>& pr = concurrent_samples[ch];
		if (pr.size() > 1)
		{
			DTWCache* dtwcache = new DTWCache(majorcompletedcount);
			std::unordered_map<sizetype, std::vector<sizetype>> clusters;
			DTWCluster::DBScan(dtwcache, pr, 5, labeldistances1, clusters);

			std::vector<float> weis;
			for (auto hs : clusters)
			{
				SequenceTemplate tdtw;
				DBA(dtwcache, pr, hs.second, labeldistances1, lamda, 100, &tdtw);
				size_t tvalcount = tdtw.Values.size();
				std::vector<datatype> ttimes(tvalcount, 0.0f);
				for (sizetype hsidx : hs.second)
				{
					SequenceChannel* psc = pr[hsidx];
					const std::vector<datatype>& psctimes = psc->Times;
					for (size_t tvi = 0; tvi < tvalcount; tvi++)
					{
						ttimes[tvi] += psctimes[tvi];
					}
					tdtw.Sources.push_back(psc);
				}

				size_t timecount = hs.second.size();
				tdtw.Times.resize(tvalcount);
				for (size_t tvi = 0; tvi < tvalcount; tvi++)
				{
					tdtw.Times[tvi] = ttimes[tvi] / timecount;
				}
				weis.push_back(tdtw.LocalWeight);
				concurrent_templates[ch].push_back(tdtw);
			}

			float weisum = std::accumulate(weis.begin(), weis.end(), 0.0f);
			for (size_t wi = 0; wi < weis.size(); wi++)
			{
				concurrent_templates[ch][wi].LocalWeight = weis[wi] / weisum;
			}

			delete dtwcache;
			dtwcache = nullptr;
		}
		else
		{
			SequenceTemplate tdtw;
			tdtw.GlobalWeight = 1.0f;
			tdtw.LocalWeight = 1.0f;
			SequenceChannel* sc0 = pr[0];
			tdtw.Values = sc0->Values;
			tdtw.Probabilities.resize(tdtw.Values.size());
			tdtw.Times = sc0->Times;
			for (size_t i = 0; i < tdtw.Values.size(); i++)
			{
				tdtw.Probabilities[i].insert(std::make_pair(sc0->Labels[i], 1.0f));
			}
			concurrent_templates[ch].push_back(tdtw);
		}
	});

	datatype ch_sumcost = 0.0;
	for (auto t : concurrent_templates)
	{
		for (const auto& tt : t.second)
		{
			ch_sumcost += tt.GlobalWeight;
		}
	}
	ch_sumcost = 1.0 / ch_sumcost;

	data->_templates.clear();

	for (auto t : concurrent_templates)
	{
		data->_templates.insert(std::make_pair(t.first, std::vector<SequenceTemplate*>{}));
		for (auto& tt : t.second)
		{
			data->_templates[t.first].push_back(SequenceTemplate::Create(tt, ch_sumcost));
		}
	}

}

